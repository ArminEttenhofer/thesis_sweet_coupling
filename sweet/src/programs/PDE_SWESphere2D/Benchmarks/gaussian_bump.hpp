/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_GAUSSIAN_BUMP_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_GAUSSIAN_BUMP_HPP


#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {

class gaussian_bump	:
		public BaseInterface
{
public:
	gaussian_bump()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				i_benchmark_name == "gaussian_bump"			||
				i_benchmark_name == "gaussian_bump_phi"			||
				i_benchmark_name == "gaussian_bump_phi_pint"		||
				i_benchmark_name == "sharp_gaussian_bump" ||
				i_benchmark_name == "gaussian_bump_vrt"			||
				i_benchmark_name == "gaussian_bump_div"			||
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

		if (benchmark_name == "gaussian_bump_phi_pint")
		{
			shackPDESWESphere2D->h0 = 29400.0;
		}
		else
		{
			shackPDESWESphere2D->h0 = 29400.0/shackPDESWESphere2D->gravitation;
			// Scale geopotential to make NL influencing the stiffness stronger
			//shackPDESWESphere2D->h0 *= 0.2;
			//shackPDESWESphere2D->gravitation *= 0.2;
		}
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
		stream << "  'gaussian_bumps_pvd':	Gaussian bump slightly delocated on pot, vrt and div field" << std::endl;
		return stream.str();
	}

public:
	sweet::Data::Sphere2D::DataGrid get_gaussian_bump(
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		sweet::Data::Sphere2D::DataGrid o_h(ops->sphere2DDataConfig);

		o_h.grid_update_lambda_gaussian_grid(
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

		return o_h;
	}


	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{


		/*
		 * We scale the variables so that cancellation errors are reduced.
		 *
		 * Div is related to geopotential via the laplace operator which requires a scaling of 1/r^2
		 *
		 * Vrt and div need to be of the same order of magnitude due to the stream formulation
		 */
		double phi_scale;
		if (benchmark_name == "gaussian_bump_phi_pint")
			phi_scale = 6000 * shackPDESWESphere2D->gravitation;
		else
			phi_scale = 0.1*shackPDESWESphere2D->h0;
		double vrt_scale = phi_scale/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);
		double div_scale = phi_scale/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);

		if (benchmark_name == "gaussian_bump" || benchmark_name == "gaussian_bump_phi")
			o_phi_pert = get_gaussian_bump(M_PI, M_PI/3, 20.0)*phi_scale;
		else if (benchmark_name == "gaussian_bump_phi_pint")
			o_phi_pert = get_gaussian_bump(M_PI, M_PI/4., 20.)*phi_scale;
		else if (benchmark_name == "sharp_gaussian_bump")
			o_phi_pert = get_gaussian_bump(M_PI, M_PI/4, 40.0)*6000;
		else
			o_phi_pert.spectral_setZero();

		if (benchmark_name == "gaussian_bump_vrt")
			o_vrt = get_gaussian_bump(M_PI + M_PI*0.1, M_PI/3 + M_PI*0.05, 20.0)*vrt_scale;
		else
			o_vrt.spectral_setZero();

		if (benchmark_name == "gaussian_bump_div")
			o_div = get_gaussian_bump(M_PI + M_PI*0.05, M_PI/3 + M_PI*0.1, 20.0)*div_scale;
		else
			o_div.spectral_setZero();
	}


	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}


};

}}

#endif
