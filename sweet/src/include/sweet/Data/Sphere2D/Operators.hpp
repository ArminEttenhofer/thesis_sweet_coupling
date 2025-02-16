/*
 *  Created on: 12 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_OPERATORS_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_OPERATORS_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Helpers/Helpers_SPHIdentities.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>

#if SWEET_THREADING
#include <omp.h>
#endif

namespace sweet {
namespace Data {
namespace Sphere2D {


class Operators	:
	public Helpers_SPHIdentities
{
	friend Config;
	typedef std::complex<double> Tcomplex;
	typedef double Treal;

public:
	Error::Base error;

	const Sphere2D::Config *sphere2DDataConfig;

private:
	double r;		///! radius
	double r2;		///! radius^2

	double ir;		///! 1/radius
	double ir2;		///! 1/radius^2


public:
	Operators(
		const Sphere2D::Config *i_sphere2DDataConfig,
		const Sphere2D::Shack *i_shackSphere2DDataOps
	);


public:
	Operators();


	Sphere2D::DataGrid getFG_fSphere2D(
			double i_fsphere2d_f0
	)	const;

	void getFG_fSphere2D(
			double i_fsphere2d_f0,
			Sphere2D::DataGrid o_fg
	)	const;


	Sphere2D::DataGrid getFG_rotatingSphere2D(
			double i_sphere2d_rotating_coriolis_omega
	)	const;


	void getFG_rotatingSphere2D(
			double i_sphere2d_rotating_coriolis_omega,
			Sphere2D::DataGrid &o_fg
	)	const;


public:
	void setup(
		const Sphere2D::Config *i_sphere2DDataConfig,
		const Sphere2D::Shack *i_shackSphere2DDataOps
	);

	void clear();


	/*!
	 * Compute gradient
	 */
	void grad_2_vec(
			const Sphere2D::DataSpectral &i_phi,
			Sphere2D::DataGrid &o_u,
			Sphere2D::DataGrid &o_v,
			double i_radius

	)	const;



	/*!
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void vrtdiv_2_uv(
			const Sphere2D::DataSpectral &i_vrt,
			const Sphere2D::DataSpectral &i_div,
			Sphere2D::DataGrid &o_u,
			Sphere2D::DataGrid &o_v

	)	const;



	/*!
	 * Convert spectral scalar field to physical one
	 */
	void scalar_spectral_2_physical(
			const Sphere2D::DataSpectral &i_spectral,
			Sphere2D::DataGrid &o_physical
	)	const;



	/*!
	 * Convert spectral scalar field to physical one
	 */
	Sphere2D::DataGrid scalar_spectral_2_physical(
			const Sphere2D::DataSpectral &i_spectral

	)	const;



	/*!
	 * Convert spectral scalar field to physical one
	 */
	void scalar_grid_2_spectral(
			const Sphere2D::DataGrid &i_physical,
			Sphere2D::DataSpectral &o_spectral
	)	const;



	/*!
	 * Convert spectral scalar field to physical one
	 */
	DataSpectral scalar_grid_2_spectral(
			const Sphere2D::DataGrid &i_physical
	)	const;



	DataSpectral uv_2_vort(
			const Sphere2D::DataGrid &i_u,
			const Sphere2D::DataGrid &i_v

	)	const;



	/*!
	 * Compute Nonlinear advection terms
	 *
	 * U \cdot \grad phi = \div \cdot (V*phi) - \nabla
	 */
	DataSpectral V_dot_grad_scalar(
			const Sphere2D::DataGrid &i_u_phys,		//!< u velocity
			const Sphere2D::DataGrid &i_v_phys,		//!< v velocity
			const Sphere2D::DataGrid &i_div_phys,		//!< divergence in physical space to avoid transformation
			const Sphere2D::DataGrid &i_scalar_phys	//!< scalar
	)	const;


	DataSpectral uv_2_div(
			const Sphere2D::DataGrid &i_u,
			const Sphere2D::DataGrid &i_v
	)	const;



	void uv_2_vrtdiv(
			const Sphere2D::DataGrid &i_u,
			const Sphere2D::DataGrid &i_v,
			Sphere2D::DataSpectral &o_vrt,
			Sphere2D::DataSpectral &o_div

	)	const;


	DataSpectral spectral_one_minus_sinphi_squared_diff_lat_mu(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;


	DataSpectral spectral_cosphi2_diff_lat_mu(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;



	/*!
	 * (1-mu^2) d/dmu ()
	 */
	DataSpectral spectral_one_minus_mu_squared_diff_lat_mu(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;



	/*!
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	DataSpectral mu(
			const Sphere2D::DataSpectral &i_sphere2d_data
	)	const;


	/*!
	 * Solve a Helmholtz problem given by
	 *
	 * (1 + b D^2) x = rhs
	 */
	DataSpectral implicit_helmholtz(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			const double &i_b,
			double i_sphere2d_radius
	)	const;



	/*!
	 * Compute multiplication with "J" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	DataSpectral implicit_J(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			double i_dt_two_omega
	)	const;


	/*!
	 * Compute multiplication with "J^-1" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	DataSpectral implicit_Jinv(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			double i_dt_two_omega
	)	const;


	/*!
	 * Compute multiplication with "F" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	DataSpectral implicit_F(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			double i_dt_two_omega
	)	const;


	/*!
	 * Compute multiplication with "F" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	DataSpectral implicit_FJinv(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			double i_dt_two_omega
	)	const;

	/*!
	 * Compute multiplication with "L" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 *
	 * return -dt * D^2 * sphere2d_data
	 */
	DataSpectral implicit_L(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			double i_dt
	)	const;


	/*!
	 * Compute multiplication with "Linv" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	DataSpectral implicit_Linv(
			const Sphere2D::DataSpectral &i_sphere2d_data,
			double i_dt
	)	const;



	/*!
	 * Compute
	 * mu*mu*F(\lambda,\mu)
	 */
	DataSpectral mu2(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;


	/*!
	 * Laplace operator
	 */
	DataSpectral laplace(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

	/*!
	 * root Laplace operator
	 */
	DataSpectral root_laplace(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

	/*!
	 * Laplace operator
	 */
	DataSpectral inv_laplace(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;


	/*!
	 * inverse root Laplace operator
	 */
	DataSpectral inv_root_laplace(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;


	/*!
	 * Calculates implicit diffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	__attribute__((deprecated))
	DataSpectral implicit_diffusion(
			const Sphere2D::DataSpectral &i_data,
			double i_coef,
			/*int i_order*/
			double i_r
	);


	/*!
	 * Calculates implicit hyperdiffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	 DataSpectral implicit_hyperdiffusion(
			const Sphere2D::DataSpectral &i_data,
			double i_coef,
			int i_order,
			double i_r
	);
};

}}}

#endif
