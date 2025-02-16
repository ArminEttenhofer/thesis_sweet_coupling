/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_2_GEOSTROPHIC_BALANCE_LINEAR_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_2_GEOSTROPHIC_BALANCE_LINEAR_HPP


#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {


class williamson_2_geostrophic_balance_linear	:
		public BaseInterface
{
	sweet::Data::Sphere2D::DataGrid fg;

public:
	williamson_2_geostrophic_balance_linear()
	{

	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson2_linear"			||
				benchmark_name == "geostrophic_balance_linear"	||
				benchmark_name == "geostrophic_balance_linear_2"	||
				benchmark_name == "geostrophic_balance_linear_4"	||
				benchmark_name == "geostrophic_balance_linear_8"	||
				benchmark_name == "geostrophic_balance_linear_16"	||
				benchmark_name == "geostrophic_balance_linear_32"	||
				benchmark_name == "geostrophic_balance_linear_64"	||
				benchmark_name == "geostrophic_balance_linear_128"	||
				benchmark_name == "geostrophic_balance_linear_256"	||
				benchmark_name == "geostrophic_balance_linear_512"	||
				false
		;
	}



	void setup_1_shackData()
	{
		if (shackParallelization->isMPIRoot)
		{
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
		}

		shackPDESWESphere2D->sphere2d_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere2D->gravitation = 9.80616;
		shackSphere2DDataOps->sphere2d_radius = 6.37122e6;
		shackPDESWESphere2D->h0 = 29400.0/shackPDESWESphere2D->gravitation;
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;

		if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
			fg = ops->getFG_fSphere2D(shackPDESWESphere2D->sphere2d_fsphere2d_f0);
		else
			fg = ops->getFG_rotatingSphere2D(shackPDESWESphere2D->sphere2d_rotating_coriolis_omega);

	}


	void clear()
	{
	}


	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  WILLIAMSON #2 (variant):" << std::endl;
		stream << "     'williamson2_linear'" << std::endl;
		stream << "     'geostrophic_balance_linear': Geostrophic balance for linear SWE" << std::endl;
		return stream.str();
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_linear(
			const sweet::Data::Sphere2D::Operators *i_ops,
			const sweet::Data::Sphere2D::DataSpectral &i_vort,
			const sweet::Data::Sphere2D::DataSpectral &i_div,
			sweet::Data::Sphere2D::DataSpectral &o_phi
	)
	{
		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = o_phi.sphere2DDataConfig;

		/*
		 * Compute vorticity and divergence from velocities
		 */
		sweet::Data::Sphere2D::DataGrid u(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid v(sphere2DDataConfig);

		i_ops->vrtdiv_2_uv(i_vort, i_div, u, v);

		sweet::Data::Sphere2D::DataGrid vrtg = i_vort.toGrid();

		using namespace sweet;
		sweet::Data::Sphere2D::DataGrid tmpg1 = u*fg;
		sweet::Data::Sphere2D::DataGrid tmpg2 = v*fg;

		sweet::Data::Sphere2D::DataSpectral tmpspec1(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral tmpspec2(sphere2DDataConfig);

		i_ops->uv_2_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		o_phi = i_ops->inv_laplace(tmpspec1);
	}



	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{

		/*
		 * Linear-only!!!
		 *
		 * The original version is for the non-linear case
		 */

		double a = shackSphere2DDataOps->sphere2d_radius;
		//double omega = shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;

		double freq_multiplier = 1.0;

		if (benchmark_name == "geostrophic_balance_linear_2")
			freq_multiplier = 2.0;
		else if (benchmark_name == "geostrophic_balance_linear_4")
			freq_multiplier = 4.0;
		else if (benchmark_name == "geostrophic_balance_linear_8")
			freq_multiplier = 8.0;
		else if (benchmark_name == "geostrophic_balance_linear_16")
			freq_multiplier = 16.0;
		else if (benchmark_name == "geostrophic_balance_linear_32")
			freq_multiplier = 32.0;
		else if (benchmark_name == "geostrophic_balance_linear_64")
			freq_multiplier = 64.0;
		else if (benchmark_name == "geostrophic_balance_linear_128")
			freq_multiplier = 128.0;
		else if (benchmark_name == "geostrophic_balance_linear_256")
			freq_multiplier = 256.0;
		else if (benchmark_name == "geostrophic_balance_linear_512")
			freq_multiplier = 512.0;


		/*
		 * Setup U
		 */
		sweet::Data::Sphere2D::DataGrid ug(ops->sphere2DDataConfig);
		ug.grid_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				// Eq. 90, Williamson TC paper
				o_data = u0*(std::cos(phi*freq_multiplier)*std::cos(alpha) + std::cos(lon)*std::sin(phi*freq_multiplier)*std::sin(alpha));
			}
		);

		/*
		 * Setup V
		 */
		sweet::Data::Sphere2D::DataGrid vg(ops->sphere2DDataConfig);
		vg.grid_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				// Eq. 91, Williamson TC paper
				o_data = -u0*std::sin(lon*freq_multiplier)*std::sin(alpha);
			}
		);

		ops->uv_2_vrtdiv(ug, vg, o_vrt, o_div);



		/**
		 * TEST for non-divergent test case
		 */
		double div_zero = o_div.toGrid().grid_reduce_max_abs();
		if (div_zero > 1e-12)
		{
			if (shackParallelization->isMPIRoot)
			{
				std::cout << "Divergence: " << div_zero << std::endl;
			}
			SWEETErrorFatal("Divergence should be close to 0, maybe there are some numerical round-off errors?");
		}

		if (freq_multiplier == 1.0)
		{
			/**
			 * TEST for correct vorticity
			 */
			/*
			 * Setup relative vorticity
			 */
			sweet::Data::Sphere2D::DataGrid vortg(ops->sphere2DDataConfig);
			vortg.grid_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 94, Williamson TC paper

					// absolute vorticity, but we like the relative one
					//o_data = (2.0*u0/a + 2.0*omega)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));

					// relative vorticity
					o_data = (2.0*u0/a)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));
				}
			);

			double vort_diff = (o_vrt.toGrid() - vortg).grid_reduce_max_abs();
			if (vort_diff > 1e-12)
			{
				if (shackParallelization->isMPIRoot)
				{
					std::cout << "Vorticity difference: " << vort_diff << std::endl;
				}
				SWEETErrorFatal("Vorticity fields differ (should be close to 0), maybe there are some numerical round-off errors?");
			}
		}


		if (shackParallelization->isMPIRoot)
			std::cout << "[MULE] use_analytical_geostrophic_setup: 1" << std::endl;

		computeGeostrophicBalance_linear(
				ops,
				o_vrt,
				o_div,
				o_phi_pert
		);

		/*
		 * TODO: Maybe we need to add
		 */
	}


	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}

};

}}

#endif
