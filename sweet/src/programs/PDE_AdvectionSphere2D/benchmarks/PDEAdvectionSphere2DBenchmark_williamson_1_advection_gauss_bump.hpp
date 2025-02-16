/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_WILLIAMSON_1_ADVECTION_GAUSS_BUMP_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_WILLIAMSON_1_ADVECTION_GAUSS_BUMP_HPP


#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmarks_BaseInterface.hpp>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>



class PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump	:
	public PDEAdvectionSphere2DBenchmarks_BaseInterface
{
	double _default_sphere2d_radius = 6.37122e6;

public:
	PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson1b"		||
				benchmark_name == "adv_gauss_bump"		||
				benchmark_name == "advection_gauss_bump"	||
				false
		;
	}

	void setup_1_shackData()
	{
		std::cout << "!!! WARNING !!!" << std::endl;
		std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
		std::cout << "!!! WARNING !!!" << std::endl;

		shackSphere2DDataOps->sphere2d_radius = _default_sphere2d_radius;
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;
	}


	std::string printHelp()
	{
		std::ostringstream stream;
		stream << "  WILLIAMSON #1 (variant):" << std::endl;
		stream << "     'williamson1b'/" << std::endl;
		stream << "     'adv_gauss_bump'/" << std::endl;
		stream << "     'adv_gauss_bump': Advection test case of gaussian bump" << std::endl;
		stream << "         OPTION:" << std::endl;
		stream << "         --benchmark-advection-rotation-angle=[angle]" << std::endl;
		return stream.str();
	}

	void getInitialState(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &o_prognostic_fields,
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
	)
	{
		sweet::Data::Sphere2D::DataSpectral vrt(ops->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div(ops->sphere2DDataConfig);

		getInitialState_Spectral(o_prognostic_fields[0], vrt, div);
		ops->vrtdiv_2_uv(vrt, div, o_u, o_v);
	}

	void getInitialState_Spectral(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		/*
		 * Alternative to original Williamson #1 advection test case which is based on a Gaussian bell instead of a cosine bell.
		 * This allows to test for L_inf convergence.
		 */

		/*
		 * Coefficients from benchmark
		 */
		double gravitation = 9.80616;
		double h0 = 1000.0;
		double gh0 = gravitation * h0;
		

		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;
		//theta_c = M_PI*0.5*0.8;
		double a = _default_sphere2d_radius;

		//double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		sweet::Data::Sphere2D::DataGrid phi_phys(o_phi_pert.sphere2DDataConfig);

		phi_phys.grid_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				// Same as Williamson test case
				double d = std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				// Variation from test case
				double i_exp_fac = 20.0;
				io_data = std::exp(-d*d*i_exp_fac)*0.1*gh0;
			}
		);

		o_phi_pert.loadSphere2DDataGrid(phi_phys);

		/*
		 * Both versions are working
		 */
#if 1
		sweet::Data::Sphere2D::DataGrid stream_function(o_phi_pert.sphere2DDataConfig);

		stream_function.grid_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = shackPDEAdvBenchmark->sphere2d_advection_rotation_angle;

				io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
			}
		);

		o_vrt = ops->laplace(stream_function);

#else

		o_vrt.grid_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = shackDict->benchmark.sphere2d_advection_rotation_angle;

				io_data = 2.0*u0/a*(-std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha) + std::sin(i_theta)*std::cos(alpha));
			}
		);
#endif
		o_div.spectral_setZero();
	}
};

#endif
