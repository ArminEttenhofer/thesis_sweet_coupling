/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP


#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {


class three_gaussian_bumps	:
		public BaseInterface
{
public:
	three_gaussian_bumps()
	{
	}


	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "three_gaussian_bumps"	||
			i_benchmark_name == "three_gaussian_bumps_phi_pint"	||
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
		if (benchmark_name == "three_gaussian_bumps_phi_pint")
			shackPDESWESphere2D->h0 = 29400.0;
		else
			shackPDESWESphere2D->h0 = 29400.0/shackPDESWESphere2D->gravitation;

#if 0
		// Scale geopotential to make NL influencing the stiffness stronger
		shackPDESWESphere2D->h0 *= 0.2;
		shackPDESWESphere2D->gravitation *= 0.2;
#endif
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
		stream << "  'three_gaussian_bumps':	Three Gaussian bumps on the geopotential field." << std::endl;
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

		double phi_scale;
		if (benchmark_name == "three_gaussian_bumps_phi_pint")
			phi_scale = 6000. * shackPDESWESphere2D->gravitation;
		else
			phi_scale = .1 * shackPDESWESphere2D->h0;


		o_phi_pert.spectral_setZero();
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.1, M_PI/3, 20.0)*0.1*phi_scale;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.6, M_PI/5.0, 80.0)*0.1*phi_scale;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.8, -M_PI/4, 360.0)*0.1*phi_scale;

		o_vrt.spectral_setZero();
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
