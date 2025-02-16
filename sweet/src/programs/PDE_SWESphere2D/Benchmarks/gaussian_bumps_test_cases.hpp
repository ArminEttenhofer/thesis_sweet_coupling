/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_GAUSSIAN_BUMPS_TEST_CASES_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_GAUSSIAN_BUMPS_TEST_CASES_HPP


#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {


class gaussian_bumps_test_cases	:
		public BaseInterface
{
public:
	gaussian_bumps_test_cases()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "gaussian_bumps_test_cases"	||
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

	}


	void clear()
	{
	}

	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  'gaussian_bumps_test_cases':	Gaussian bumps, special test case" << std::endl;
		return stream.str();
	}


public:
	void setup_gaussian_bump(
			sweet::Data::Sphere2D::DataGrid &o_data,
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		o_data.grid_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				// https://en.wikipedia.org/wiki/Great-circle_distance
				// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
				// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

				double phi1 = asin(mu);
				double phi2 = i_center_lat;
				double lambda1 = lon;
				double lambda2 = i_center_lon;

				double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

				o_data = std::exp(-d*d*i_exp_fac);
			}
		);
	}


	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		sweet::Data::Sphere2D::DataGrid tmp(o_phi_pert.sphere2DDataConfig);

		setup_gaussian_bump(tmp, 2.0*M_PI*0.1, M_PI/3, 1.0);
		o_phi_pert.loadSphere2DDataGrid(tmp);
		o_phi_pert *= 0.1;
		o_phi_pert += shackPDESWESphere2D->h0*shackPDESWESphere2D->gravitation;

		setup_gaussian_bump(tmp, 2.0*M_PI*0.1, M_PI/3, 1.0);
		o_vrt.loadSphere2DDataGrid(tmp);
		o_vrt *= -1e-8;
		//o_vort *= 0;
		setup_gaussian_bump(tmp, 2.0*M_PI*0.1, M_PI/3, 1.0);
		o_div.loadSphere2DDataGrid(tmp);
		o_div *= 1e-8;

		/*
		 * Convert forward/backward to velocity space to apply a certain truncation
		 */
		sweet::Data::Sphere2D::DataGrid ug(o_phi_pert.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid vg(o_phi_pert.sphere2DDataConfig);
		ops->vrtdiv_2_uv(o_vrt, o_div, ug, vg);
		ops->uv_2_vrtdiv(ug, vg, o_vrt, o_div);
	}


	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}

};

}}

#endif
