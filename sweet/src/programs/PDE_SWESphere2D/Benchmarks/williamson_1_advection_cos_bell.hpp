/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_1_ADVECTION_COS_BELL_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_1_ADVECTION_COS_BELL_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "BaseInterface.hpp"
#include <ostream>

namespace PDE_SWESphere2D {
namespace Benchmarks {


class williamson_1_advection_cos_bell	:
		public BaseInterface
{
public:
	williamson_1_advection_cos_bell()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson1"		||
				benchmark_name == "adv_cosine_bell"	||
				benchmark_name == "advection_cosine_bell"	||
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
		shackPDESWESphere2D->h0 = 1000.0;
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
		stream << "  WILLIAMSON #1:" << std::endl;
		stream << "     'williamson1'" << std::endl;
		stream << "     'advection_cosine_bell'" << std::endl;
		stream << "     'adv_cosine_bell': Advection test case of cosine bell" << std::endl;
		stream << "         OPTION:" << std::endl;
		stream << "         --benchmark-advection-rotation-angle=[angle]" << std::endl;
		return stream.str();
	}


	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

		/*
		 * Advection test case
		 * See Williamson test case, eq. (77), (78), (79)
		 */


		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;
		double a = shackSphere2DDataOps->sphere2d_radius;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		sweet::Data::Sphere2D::DataGrid phi_phys(o_phi_pert.sphere2DDataConfig);

		phi_phys.grid_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double r = a * std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				if (r < R)
					io_data = shackPDESWESphere2D->h0/2.0*(1.0+std::cos(M_PI*r/R));
				else
					io_data = 0;

				io_data *= shackPDESWESphere2D->gravitation;
			}
		);

		o_phi_pert.loadSphere2DDataGrid(phi_phys);

		sweet::Data::Sphere2D::DataGrid stream_function(o_phi_pert.sphere2DDataConfig);

		stream_function.grid_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = shackPDESWEBenchmarks->benchmark_sphere2d_advection_rotation_angle;

				io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
			}
		);

		o_vrt = ops->laplace(stream_function);
		o_div.spectral_setZero();

		o_phi_pert -= gh0;
	}



	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}

};

}}

#endif
