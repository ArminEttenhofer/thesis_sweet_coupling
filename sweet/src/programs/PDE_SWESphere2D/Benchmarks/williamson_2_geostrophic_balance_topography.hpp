/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_2_TOPOGRAPHY_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_2_TOPOGRAPHY_HPP


#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {


class williamson_2_geostrophic_balance_topography	:
		public BaseInterface
{
	HelperGeostropicBalance helperGeostropicBalance;

public:
	williamson_2_geostrophic_balance_topography()
	{
	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		BaseInterface::shackRegistration(io_shackDict);
		helperGeostropicBalance.shackRegistration(io_shackDict);
		return true;
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson2_topography"			||
				benchmark_name == "geostrophic_balance_topography"	||
				benchmark_name == "geostrophic_balance_1_topography"	||
				benchmark_name == "geostrophic_balance_2_topography"	||
				benchmark_name == "geostrophic_balance_4_topography"	||
				benchmark_name == "geostrophic_balance_8_topography"	||
				benchmark_name == "geostrophic_balance_16_topography"	||
				benchmark_name == "geostrophic_balance_32_topography"	||
				benchmark_name == "geostrophic_balance_64_topography"	||
				benchmark_name == "geostrophic_balance_128_topography"	||
				benchmark_name == "geostrophic_balance_256_topography"	||
				benchmark_name == "geostrophic_balance_512_topography"	||
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
		//shackPDESWESphere2D->h0 = 29400.0/shackPDESWESphere2D->gravitation;
		///shackPDESWESphere2D->h0 = 100.;
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;

		helperGeostropicBalance.setup(ops);

	}



	void clear()
	{
	}


	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  WILLIAMSON #2 with topography:" << std::endl;
		stream << "     'williamson2_topography'" << std::endl;
		stream << "     'geostrophic_balance': Geostrophic balance, one wave (standard)" << std::endl;
		stream << "     'geostrophic_balance_[N]': Geostrophic balance, with [N] waves" << std::endl;
		return stream.str();
	}


	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{

		/*
		 * geostrophic_balance / geostrophic_balance_1:
		 * Williamson test case 2 for geostrophic balance.
		 *
		 * WARNING: This test uses a balanced solution for the full non-linear equations
		 * See Williamson paper for accurate setup
		 *
		 * "geostrophic_balance_N" means that N is the multiplier for the frequency
		 * in the direction of the Latitude
		 */
		double a = shackSphere2DDataOps->sphere2d_radius;
		double omega = shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;


		double freq_multiplier = 1.0;

		if (benchmark_name == "geostrophic_balance_2_topography")
			freq_multiplier = 2.0;
		else if (benchmark_name == "geostrophic_balance_4_topography")
			freq_multiplier = 4.0;
		else if (benchmark_name == "geostrophic_balance_8_topography")
			freq_multiplier = 8.0;
		else if (benchmark_name == "geostrophic_balance_16_topography")
			freq_multiplier = 16.0;
		else if (benchmark_name == "geostrophic_balance_32_topography")
			freq_multiplier = 32.0;
		else if (benchmark_name == "geostrophic_balance_64_topography")
			freq_multiplier = 64.0;
		else if (benchmark_name == "geostrophic_balance_128_topography")
			freq_multiplier = 128.0;
		else if (benchmark_name == "geostrophic_balance_256_topography")
			freq_multiplier = 256.0;
		else if (benchmark_name == "geostrophic_balance_512_topography")
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

		/*
		 * Setup Phi
		 */
		sweet::Data::Sphere2D::DataGrid phig(ops->sphere2DDataConfig);
		phig.grid_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				// Eq. 91, Williamson TC paper
				o_data = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0 * 0.;
			}
		);


		///o_phi_pert.loadSphere2DDataGrid(phig + shackPDESWESphere2D->gravitation*shackPDESWEBenchmarks->topography.topography_grid);
		o_phi_pert.loadSphere2DDataGrid(phig);


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

	}
public:
	void setup_topography() override
	{
		if (shackDict == nullptr)
			SWEETErrorFatal("Benchmarks are not yet initialized");

		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);

		double a = shackSphere2DDataOps->sphere2d_radius;
		double g = shackPDESWESphere2D->gravitation;
		double omega = shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;

		double freq_multiplier = 1.0;

		if (benchmark_name == "geostrophic_balance_2_topography")
			freq_multiplier = 2.0;
		else if (benchmark_name == "geostrophic_balance_4_topography")
			freq_multiplier = 4.0;
		else if (benchmark_name == "geostrophic_balance_8_topography")
			freq_multiplier = 8.0;
		else if (benchmark_name == "geostrophic_balance_16_topography")
			freq_multiplier = 16.0;
		else if (benchmark_name == "geostrophic_balance_32_topography")
			freq_multiplier = 32.0;
		else if (benchmark_name == "geostrophic_balance_64_topography")
			freq_multiplier = 64.0;
		else if (benchmark_name == "geostrophic_balance_128_topography")
			freq_multiplier = 128.0;
		else if (benchmark_name == "geostrophic_balance_256_topography")
			freq_multiplier = 256.0;
		else if (benchmark_name == "geostrophic_balance_512_topography")
			freq_multiplier = 512.0;

		/*
		 * Setup topography
		 */
		sweet::Data::Sphere2D::DataGrid ug(ops->sphere2DDataConfig);
		shackPDESWEBenchmarks->topography.topography_grid.grid_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				// Eq. 90, Williamson TC paper
				////o_data = 1. / g * (a * omega * u0 + .5 * u0 * u0) * std::sin(phi*freq_multiplier)*std::sin(phi*freq_multiplier)*std::cos(alpha);
				o_data = 2.94e4 / g - 1. / g * (a * omega * u0 + .5 * u0 * u0) * std::sin(phi*freq_multiplier)*std::sin(phi*freq_multiplier)*std::cos(alpha);
				///o_data = u0*(std::cos(phi*freq_multiplier)*std::cos(alpha) + std::cos(lon)*std::sin(phi*freq_multiplier)*std::sin(alpha));
			}
		);

		shackPDESWEBenchmarks->topography.topography.loadSphere2DDataGrid(shackPDESWEBenchmarks->topography.topography_grid);

		shackPDESWEBenchmarks->use_topography = true;

	}





};

}}

#endif
