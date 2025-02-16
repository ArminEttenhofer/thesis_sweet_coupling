/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_HELPERGEOSTROPICBALANCE_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_HELPERGEOSTROPICBALANCE_HPP


#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "../Shack.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {

class HelperGeostropicBalance
{

	sweet::Shacks::Dictionary *shackDict;
	sweet::Data::Sphere2D::Operators *ops;

	PDE_SWESphere2D::Shack *shackPDESWESphere2D;

	sweet::Data::Sphere2D::DataGrid fg;

public:
	HelperGeostropicBalance()	:
		shackDict(nullptr),
		ops(nullptr),
		shackPDESWESphere2D(nullptr)
	{
	}

	void shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPDESWESphere2D = shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	}

	void setup(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;

		if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
			fg = ops->getFG_fSphere2D(shackPDESWESphere2D->sphere2d_fsphere2d_f0);
		else
			fg = ops->getFG_rotatingSphere2D(shackPDESWESphere2D->sphere2d_rotating_coriolis_omega);
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_nonlinear(
			sweet::Data::Sphere2D::DataSpectral &i_vort,
			sweet::Data::Sphere2D::DataSpectral &i_div,
			sweet::Data::Sphere2D::DataSpectral &o_phi
	)
	{
		/*
		 * Compute vorticity and divergence from velocities
		 */
		sweet::Data::Sphere2D::DataGrid ug(o_phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid vg(o_phi.sphere2DDataConfig);

		ops->vrtdiv_2_uv(i_vort, i_div, ug, vg);

		sweet::Data::Sphere2D::DataGrid vrtg = i_vort.toGrid();

		sweet::Data::Sphere2D::DataGrid tmpg1 = ug*(vrtg+fg);
		sweet::Data::Sphere2D::DataGrid tmpg2 = vg*(vrtg+fg);

		sweet::Data::Sphere2D::DataSpectral tmpspec1(o_phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral tmpspec2(o_phi.sphere2DDataConfig);

		ops->uv_2_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		sweet::Data::Sphere2D::DataSpectral phispec = ops->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);

		o_phi = phispec;
	}



	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_linear(
			sweet::Data::Sphere2D::DataSpectral &i_vort,
			sweet::Data::Sphere2D::DataSpectral &i_div,
			sweet::Data::Sphere2D::DataSpectral &o_phi
	)
	{
		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = o_phi.sphere2DDataConfig;
		/*
		 * Setup Coriolis effect
		 */
		sweet::Data::Sphere2D::DataGrid f(sphere2DDataConfig);
		f.grid_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega*mu;
			}
		);

		/*
		 * Compute vorticity and divergence from velocities
		 */
		sweet::Data::Sphere2D::DataGrid u(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid v(sphere2DDataConfig);

		ops->vrtdiv_2_uv(i_vort, i_div, u, v);

		sweet::Data::Sphere2D::DataGrid vrtg = i_vort.toGrid();

		sweet::Data::Sphere2D::DataGrid tmpg1 = u*f;
		sweet::Data::Sphere2D::DataGrid tmpg2 = v*f;

		sweet::Data::Sphere2D::DataSpectral tmpspec1(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral tmpspec2(sphere2DDataConfig);

		ops->uv_2_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		sweet::Data::Sphere2D::DataSpectral phispec = ops->inv_laplace(tmpspec1);

		o_phi = shackPDESWESphere2D->h0*shackPDESWESphere2D->gravitation + phispec;
	}


	void clear()
	{
		fg.clear();
	}

};

}}

#endif
