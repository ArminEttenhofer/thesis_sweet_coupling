/*
 * Author: Francois Hamon, Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_5_FLOW_OVER_MOUNTAIN_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_5_FLOW_OVER_MOUNTAIN_HPP


#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {


class williamson_5_flow_over_mountain	:
		public BaseInterface
{
	////sweet::Shacks::Dictionary *shackDict = nullptr;
	////sweet::Data::Sphere2D::Operators *ops = nullptr;

	HelperGeostropicBalance helperGeostropicBalance;

public:
	williamson_5_flow_over_mountain()
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
				benchmark_name == "williamson5"	||
				benchmark_name == "flow_over_mountain" ||
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

		//! Setup Williamson's parameters
		shackPDESWESphere2D->sphere2d_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere2D->gravitation = 9.80616;
		shackSphere2DDataOps->sphere2d_radius = 6.37122e6;
		shackPDESWESphere2D->h0 = 5960;
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
		shackPDESWEBenchmarks->topography.topography_grid.clear();
	}


	std::string getHelp()
	{
		std::ostringstream stream;

		if (shackParallelization->isMPIRoot)
		{
			stream << "  WILLIAMSON #5:" << std::endl;
			stream << "     'williamson5'/'flow_over_mountain': Flow over mountain benchmark" << std::endl;
		}
		return stream.str();
	}



	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{

		double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;

		const double u0 = 20.0;

		/*
		 * Setup V=0
		 */
		sweet::Data::Sphere2D::DataGrid vg(o_phi_pert.sphere2DDataConfig);
		vg.grid_setZero();

		/*
		 * Setup U=...
		 * initial velocity along longitude
		 */
		sweet::Data::Sphere2D::DataGrid ug(o_phi_pert.sphere2DDataConfig);
		ug.grid_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				o_data = u0 * std::cos(phi);
			}
		);

		ops->uv_2_vrtdiv(ug, vg, o_vrt, o_div);

		sweet::Data::Sphere2D::DataGrid hg(o_phi_pert.sphere2DDataConfig);

		bool use_analytical_geostrophic_setup = false;
		if (use_analytical_geostrophic_setup)
		{
			helperGeostropicBalance.computeGeostrophicBalance_nonlinear(
					o_vrt,
					o_div,
					o_phi_pert
			);

			o_phi_pert -= gh0;
		}
		else
		{

			double a = shackSphere2DDataOps->sphere2d_radius;
			double omega = shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
			double alpha = 0;
			// Squared term in Eq. 95, Williamson TC paper
			sweet::Data::Sphere2D::DataGrid r2(ops->sphere2DDataConfig);
			r2.grid_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					o_data = -std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha);
					o_data = o_data*o_data;
				}
			);

			// Eq. 95, Williamson TC paper
			sweet::Data::Sphere2D::DataGrid phig = -(a*omega*u0 + u0*u0/2.0)*r2;

			////o_phi_pert.loadSphere2DDataGrid(phig);
			o_phi_pert.loadSphere2DDataGrid(phig - shackPDESWESphere2D->gravitation*shackPDESWEBenchmarks->topography.topography_grid);
		}


	}


public:
	static
	void setup_topography_(
			sweet::Data::Sphere2D::DataGrid &o_h_topo,
			double i_R            = M_PI/9.,
			double i_h_topo_0     = 2000.,
			double i_center_lon   = 3.*M_PI/2.,
			double i_center_lat   = M_PI/6.
	)
	{
		const double center_lat = i_center_lat;
		const double center_lon = i_center_lon;

		o_h_topo.grid_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{

				const double phi1    = asin(mu);
				const double phi2    = center_lat;
				const double lambda1 = lon;
				const double lambda2 = center_lon;

				const double r_squared = std::min( i_R*i_R, (phi1-phi2)*(phi1-phi2) + (lambda1-lambda2)*(lambda1-lambda2) );

				o_data = i_h_topo_0 * ( 1. - sqrt(r_squared) / i_R );
			}
		);
	}

public:
	void setup_topography() override
	{
		if (shackDict == nullptr)
			SWEETErrorFatal("Benchmarks are not yet initialized");

		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);

		if (benchmark_name == "flow_over_mountain" || benchmark_name == "williamson5")
		{
			// set the topography flag to true
			shackPDESWEBenchmarks->use_topography = true;

			// setup the parameters for the flow-over-mountain test case
			const double R			= M_PI/9.;
			const double h_topo_0	 = 2000.;
			const double i_center_lon = 3.*M_PI/2.;
			const double i_center_lat = M_PI/6.;

			shackPDESWEBenchmarks->topography.topography_grid.grid_setZero();

			// setup the topography vector
			setup_topography_(
					shackPDESWEBenchmarks->topography.topography_grid,
					R,
					h_topo_0,
					i_center_lon,
					i_center_lat
			);

			shackPDESWEBenchmarks->topography.topography.loadSphere2DDataGrid(shackPDESWEBenchmarks->topography.topography_grid);
		}
		else
		{
			// set the topography flag to false
			shackPDESWEBenchmarks->use_topography = false;
		}
	}



};

}}

#endif
