/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 */


#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DTimeSteppers.hpp>
#include <programs/PDE_AdvectionSphere2D/ProgramPDEAdvectionSphere2D.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include "../programs/PDE_AdvectionSphere2D/time/ShackPDEAdvectionSphere2DTimeDiscretization.hpp"




class InterpolationTests
{
public:
	sweet::Data::Sphere2D::Config *sphere2DDataConfig;
	sweet::Data::Sphere2D::Config *sphere2DDataConfigOversampled;

	sweet::Data::Sphere2D::DataSpectral prog_h;

	int interpolation_order = 3;

	bool use_limiter = false;
	bool use_poles_pseudo_points = false;

	sweet::Data::Vector::Vector<double> posx_a, posy_a;

	#define MAX_GAUSSIANS 9
	double gaussian_center_array[MAX_GAUSSIANS][2] = {
			{0.0, M_PI*0.5},	// top
			{0.0, 0.0},		// equator
			{0.0, -M_PI*0.5},	// bottom

			{1.1, M_PI*0.5*0.8},	// a-kind top
			{3.0, 0.1},		// a-kind equator
			{2.0, -M_PI*0.5*0.75},	// a-kind bottom

			{0.1, M_PI*0.3},	// misc
			{1.0, M_PI*0.4},	// misc
			{2.3, -M_PI*0.34},	// misc
	};

	int gaussian_id = 0;
	//double center_lon = 0;
	//double center_lat = 0;
	double exp_fac = 20.0;

	sweet::Data::Sphere2D::Operators_Sampler_DataGrid sphere2DDataSampler;


public:
	InterpolationTests()	:
		sphere2DDataConfig(nullptr),
		sphere2DDataConfigOversampled(nullptr)
	{
	}

	void setup(
			sweet::Data::Sphere2D::Config &i_sphere2DDataConfig,
			sweet::Data::Sphere2D::Config &i_sphere2DDataConfigOversampled
	)
	{
		sphere2DDataConfig = &i_sphere2DDataConfig;
		sphere2DDataConfigOversampled = &i_sphere2DDataConfigOversampled;

		prog_h.setup(sphere2DDataConfig);

		sweet::Data::Sphere2D::DataSpectral tmp_vort(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral tmp_div(sphere2DDataConfig);

		sweet::Data::Sphere2D::DataGrid prog_h_phys(sphere2DDataConfig);

		prog_h_phys.grid_update_lambda(
			[&](double i_lon, double i_lat, double &o_data)
			{
				o_data = gaussianValue(i_lon, i_lat, exp_fac);
			}
		);
		prog_h.loadSphere2DDataGrid(prog_h_phys);

		posx_a.setup(sphere2DDataConfigOversampled->grid_number_elements);
		posy_a.setup(sphere2DDataConfigOversampled->grid_number_elements);

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position
		posx_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % sphere2DDataConfigOversampled->grid_num_lon;

				io_data = 2.0*M_PI*(double)i/(double)sphere2DDataConfigOversampled->grid_num_lon;
				SWEET_ASSERT(io_data >= 0);
				SWEET_ASSERT(io_data < 2.0*M_PI);
			}
		);
		posy_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				//int i = idx % sphere2DDataConfig->grid_data_size[0];
				int j = idx / sphere2DDataConfigOversampled->grid_num_lon;

				io_data = sphere2DDataConfigOversampled->lat[j];

				SWEET_ASSERT(io_data >= -M_PI*0.5);
				SWEET_ASSERT(io_data <= M_PI*0.5);
			}
		);

		sphere2DDataSampler.setup(sphere2DDataConfig);
	}


	double gaussianValue(
			double i_lon, double i_lat,
			double i_exp_fac
	)
	{
		if (gaussian_id >= 0)
		{
			return gaussianValue_(gaussian_center_array[gaussian_id][0], gaussian_center_array[gaussian_id][1], i_lon, i_lat, i_exp_fac);
		}
		else
		{
			double o_data = 0;

			for (int i_gaussians = 0; i_gaussians < MAX_GAUSSIANS; i_gaussians++)
			{
				double scalar = std::cos(i_lat*2)*std::cos(i_lon*4);
				scalar = 1.0;
				o_data += scalar*gaussianValue_(gaussian_center_array[i_gaussians][0], gaussian_center_array[i_gaussians][1], i_lon, i_lat, i_exp_fac);
				//o_data = scalar;
			}

			o_data /= MAX_GAUSSIANS;
			return o_data;
		}
	}

	double gaussianValue_(
			double i_center_lon, double i_center_lat,
			double i_lon, double i_lat,
			double i_exp_fac
	)
	{
#if 1

		double x0[3];
		sweet::LibMath::VectorMath::point_latlon_2_cartesian__scalar(i_center_lon, i_center_lat, x0[0], x0[1], x0[2]);

		double x[3];
		sweet::LibMath::VectorMath::point_latlon_2_cartesian__scalar(i_lon, i_lat, x[0], x[1], x[2]);

#if 0
		double d =	(x[0] - x0[0])*(x[0] - x0[0])*(2.0+i_lon*0.1) +
					(x[1] - x0[1])*(x[1] - x0[1])*(2.0+i_lat*0.1) +
					(x[2] - x0[2])*(x[2] - x0[2])*(2.0+i_lon*i_lat*0.1);
#else

		double d =	(x[0] - x0[0])*(x[0] - x0[0])*(2.0) +
					(x[1] - x0[1])*(x[1] - x0[1])*(3.0) +
					(x[2] - x0[2])*(x[2] - x0[2])*(4.0);
#endif

		return std::exp(-20*d);

#else
		double center_lon = i_center_lon;
		double center_lat = i_center_lat;

		double mu = std::sin(i_lat);
		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = i_lon;
		double lambda2 = center_lon;

		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		double d1 = acos(sin(phi1)*sin(phi2));
		double d2 = acos(cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		//return std::exp(-d*d*i_exp_fac)*(1.0-(d1*d1)/(M_PI*M_PI));//*0.1*simVars.sim.h0;// + simVars.sim.h0;
		return std::exp(-d*d*i_exp_fac);//*0.1*simVars.sim.h0;// + simVars.sim.h0;
#endif
	}


	double runTests()
	{
		sweet::Data::Vector::Vector<double> out_data(posx_a.numberOfElements);

		if (interpolation_order == 2)
		{
			sphere2DDataSampler.bilinear_scalar(
					prog_h.toGrid(),
					posx_a,
					posy_a,
					out_data.data,
					false,
					use_poles_pseudo_points
			);
		}
		else if (interpolation_order == 3)
		{
			sphere2DDataSampler.bicubic_scalar(
					prog_h.toGrid(),
					posx_a,
					posy_a,
					out_data.data,
					false,
					use_poles_pseudo_points,
					use_limiter
			);
		}
		else
		{
			SWEETErrorFatal("Interpolation order not available");
		}

		double max_error = 0;
		SWEET_ASSERT(posx_a.numberOfElements != 0);
		for (std::size_t i = 0; i < posx_a.numberOfElements; i++)
		{
			double value = gaussianValue(posx_a.data[i], posy_a.data[i], exp_fac);
			double err = std::abs(value - out_data.data[i]);
			max_error = std::max(max_error, err);
		}
		return max_error;
	}


	bool should_quit()
	{
		return false;
	}
};


int main(int i_argc, char *i_argv[])
{
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	int initial_spectral_modes = shackSphere2DDataOps->space_res_spectral[0];

	for (
			int use_poles_pseudo_points = 0;
			use_poles_pseudo_points < 2;
			use_poles_pseudo_points++
	)
	{
		std::cout << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << "* Running studies with or without pseudo pole points " << use_poles_pseudo_points << std::endl;
		std::cout << "*********************************************************" << std::endl;

		//int i_gaussians_start = -1;

		for (int gaussian_id = -1; gaussian_id < 9; gaussian_id++)
		//for (int i_gaussians = 8; i_gaussians >= -1; i_gaussians--)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* Running studies for Gaussian type " << gaussian_id << std::endl;
			std::cout << "*********************************************************" << std::endl;

			for (int interpolation_order = 2; interpolation_order <= 3; interpolation_order++)
			{
				std::cout << std::endl;
				std::cout << "*********************************************************" << std::endl;
				std::cout << "* Running studies for interpolation of order " << interpolation_order << std::endl;
				std::cout << "*********************************************************" << std::endl;

				int oversampling = 5;
				std::cout << "Using oversampling of " << oversampling << std::endl;

				double prev_max_error = -1;
				for (int i = initial_spectral_modes; i <= 256*2; i *= 2)
				{

					shackSphere2DDataOps->space_res_physical[0] = 2*i;
					shackSphere2DDataOps->space_res_physical[1] = i;

					shackSphere2DDataOps->space_res_spectral[0] = i;
					shackSphere2DDataOps->space_res_spectral[1] = i;

					sweet::Data::Sphere2D::Config sphere2DDataConfig;
					sphere2DDataConfig.setupAuto(shackSphere2DDataOps);

					/*
					 * Generate higher resolution
					 */
					sweet::Data::Sphere2D::Shack shackSphere2DDataOpsOversampled;
					shackSphere2DDataOpsOversampled = *shackSphere2DDataOps;

					shackSphere2DDataOpsOversampled.space_res_physical[0] *= oversampling;
					shackSphere2DDataOpsOversampled.space_res_physical[1] *= oversampling;

					shackSphere2DDataOpsOversampled.space_res_spectral[0] *= oversampling;
					shackSphere2DDataOpsOversampled.space_res_spectral[1] *= oversampling;

					sweet::Data::Sphere2D::Config sphere2DDataConfigOversampled;
					sphere2DDataConfigOversampled.setupAuto(shackSphere2DDataOpsOversampled);


					{
						InterpolationTests interpolationTests;

						// Update interpolation order
						interpolationTests.interpolation_order = interpolation_order;

						// center of Gaussian bump
						interpolationTests.gaussian_id = gaussian_id;

						interpolationTests.use_poles_pseudo_points = use_poles_pseudo_points;

						interpolationTests.setup(sphere2DDataConfig, sphere2DDataConfigOversampled);

						double max_error = interpolationTests.runTests();
						{
							std::cout << "Lmax error: " << max_error << std::endl;

							if (prev_max_error >= 0)
							{
								//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
								double conv = prev_max_error / max_error;
								std::cout << "Convergence: " << conv << std::endl;

								if (use_poles_pseudo_points == 1 && gaussian_id == -1)
								{
									if (conv*1.2 < std::pow(2.0, interpolation_order))
										SWEETErrorFatal("Convergence not given!");
								}
								else
								{
									if (conv*1.1 < std::pow(2.0, interpolation_order))
										SWEETErrorFatal("Convergence not given!");
								}
							}
							prev_max_error = max_error;
						}
					}
				}
			}
		}
	}

	return 0;
}
