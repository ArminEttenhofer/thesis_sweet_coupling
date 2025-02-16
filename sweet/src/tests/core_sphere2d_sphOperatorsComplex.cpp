/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Data/Sphere2DComplex/Sphere2DComplex.hpp>
#include <sweet/Error/Fatal.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>

#include "core_sphere2d_sphOperators/Sphere2DTestSolutions_Gaussian.hpp"


void test_header(const std::string &i_str)
{
	std::cout << "**********************************************" << std::endl;
	std::cout << i_str << std::endl;
	//std::cout << "**********************************************" << std::endl;
}

void run_tests(
		sweet::Data::Sphere2D::Config *i_sphere2DDataConfig,
		sweet::Data::Sphere2D::Shack *i_shackSphere2DDataOps
)
{
	//double eps = 1e-10;
	double eps = 1e-8;
	eps *= std::sqrt(i_sphere2DDataConfig->spectral_modes_n_max)*std::sqrt(i_sphere2DDataConfig->spectral_modes_m_max);
	std::cout << "Using max allowed error of " << eps << std::endl;

	// Use earth radius of 1
	i_shackSphere2DDataOps->sphere2d_radius = 1.0;
	sweet::Data::Sphere2DComplex::Operators ops(i_sphere2DDataConfig, i_shackSphere2DDataOps);


	if (true)
	{
		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		double alpha[] = {0, M_PI/3, M_PI/2};

		{
			test_header("Testing divergence freeness (non-Robert)");

			for (int i = 0; i < sizeof(alpha)/sizeof(double); i++)
			{
				double advection_rotation_angle = alpha[i];

				std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

				sweet::Data::Sphere2DComplex::DataGrid u(i_sphere2DDataConfig);
				u.grid_update_lambda(
					[&](double i_lon, double i_lat, std::complex<double> &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
								u0*(
									std::cos(i_theta)*std::cos(advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(advection_rotation_angle)
							);
					}
				);

				sweet::Data::Sphere2DComplex::DataGrid v(i_sphere2DDataConfig);
				v.grid_update_lambda(
					[&](double i_lon, double i_lat, std::complex<double> &io_data)
					{
						//double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(advection_rotation_angle)
							);
					}
				);

				sweet::Data::Sphere2DComplex::DataSpectral vort(i_sphere2DDataConfig);
				sweet::Data::Sphere2DComplex::DataSpectral div(i_sphere2DDataConfig);
				ops.uv_2_vrtdiv(u, v, vort, div);

				double div_max_error = div.toGrid().grid_reduce_max_abs();
				std::cout << " + div_max_error: " << div_max_error << std::endl;

				if (div_max_error > eps)
					SWEETErrorFatal(" + ERROR! max error exceeds threshold");
			}
		}
	}


	{
		Sphere2DTestSolutions_Gaussian testSolutions;

		if (true)
		{
			test_header("Testing Multiplication (a*b) with b=123.0");

			sweet::Data::Sphere2DComplex::DataSpectral data(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid data_phys(i_sphere2DDataConfig);

			data_phys.grid_update_lambda(
				[&](double x, double y, std::complex<double> &io_data)
				{
					io_data = y;
				}
			);

			data.loadSphere2DDataGrid(data_phys);
			data = data*123.0;

			sweet::Data::Sphere2DComplex::DataSpectral data2(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid data2_phys(i_sphere2DDataConfig);
			data2_phys.grid_update_lambda(
				[&](double x, double y, std::complex<double> &io_data)
				{
					io_data = y*123.0;
				}
			);
			data2.loadSphere2DDataGrid(data2_phys);

			double div_max_error = (data-data2).toGrid().grid_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETErrorFatal(" + ERROR! max error exceeds threshold");
		}



		{
			Sphere2DTestSolutions_Gaussian testSolutions;

			if (true)
			{
				test_header("Testing Multiplication (a *= b) with b=123.0");

				sweet::Data::Sphere2DComplex::DataSpectral data(i_sphere2DDataConfig);
				sweet::Data::Sphere2DComplex::DataGrid data_phys(i_sphere2DDataConfig);
				data_phys.grid_update_lambda(
						[&](double x, double y, std::complex<double> &io_data)
						{
							io_data = y;
						}
				);
				data.loadSphere2DDataGrid(data_phys);

				data *= 123.0;

				sweet::Data::Sphere2DComplex::DataSpectral data2(i_sphere2DDataConfig);
				sweet::Data::Sphere2DComplex::DataGrid data2_phys(i_sphere2DDataConfig);
				data2_phys.grid_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y*123.0;
					}
				);
				data2.loadSphere2DDataGrid(data2_phys);

				double div_max_error = (data-data2).toGrid().grid_reduce_max_abs();
				std::cout << " + div_max_error: " << div_max_error << std::endl;

				if (div_max_error > eps)
					SWEETErrorFatal(" + ERROR! max error exceeds threshold");
			}
		}

		if (true)
		{
			test_header("Testing add (a+b) operation with 123.0");

			sweet::Data::Sphere2DComplex::DataSpectral data(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid data_phys(i_sphere2DDataConfig);
			data_phys.grid_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y;
					}
			);

			data.loadSphere2DDataGrid(data_phys);
			data = data + 123.0;

			sweet::Data::Sphere2DComplex::DataSpectral data2(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid data2_phys(i_sphere2DDataConfig);
			data2_phys.grid_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y+123.0;
					}
			);
			data2.loadSphere2DDataGrid(data2_phys);

			double div_max_error = (data-data2).toGrid().grid_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETErrorFatal(" + ERROR! max error exceeds threshold");
		}


		if (true)
		{
			test_header("Testing add (a+=b) operation with 123.0");

			sweet::Data::Sphere2DComplex::DataSpectral data(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid data_phys(i_sphere2DDataConfig);
			data_phys.grid_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y;
					}
			);
			data.loadSphere2DDataGrid(data_phys);

			data += 123.0;

			sweet::Data::Sphere2DComplex::DataSpectral data2(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid data2_phys(i_sphere2DDataConfig);
			data2_phys.grid_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y+123.0;
					}
			);
			data2.loadSphere2DDataGrid(data2_phys);

			double div_max_error = (data-data2).toGrid().grid_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETErrorFatal(" + ERROR! max error exceeds threshold");
		}

		if (true)
		{
			test_header("Testing Gaussian latitude coordinates");

			sweet::Data::Sphere2DComplex::DataSpectral h(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid h_phys(i_sphere2DDataConfig);
			h_phys.grid_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function__grid_gaussian(a,b,val);
						c = val;
					}
			);
			h.loadSphere2DDataGrid(h_phys);

			sweet::Data::Sphere2DComplex::DataSpectral hphi(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid hphi_phys(i_sphere2DDataConfig);
			hphi_phys.grid_update_lambda(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function_phi__grid_phi(a,b,val);
						c = val;
					}
			);
			hphi.loadSphere2DDataGrid(hphi_phys);

			double div_max_error = (h-hphi).toGrid().grid_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETErrorFatal(" + ERROR! max error exceeds threshold");
		}


		if (true)
		{
			test_header("Testing multiplication with Gaussian latitude");

			// mu*F(\lambda,\mu)
			sweet::Data::Sphere2DComplex::DataSpectral h(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid h_phys(i_sphere2DDataConfig);
			h_phys.grid_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function__grid_gaussian(a,b,val);
						c = val;
					}
			);
			h.loadSphere2DDataGrid(h_phys);
			h = ops.mu(h);

			sweet::Data::Sphere2DComplex::DataSpectral result(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid result_phys(i_sphere2DDataConfig);
			result_phys.grid_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.correct_result_mu__grid_gaussian(a,b,val);
						c = val;
					}
			);
			result.loadSphere2DDataGrid(result_phys);

			double div_max_error = (h-result).toGrid().grid_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETErrorFatal(" + ERROR! max error exceeds threshold");
		}

		if (true)
		{
			test_header("Testing multiplication with pow2 of Gaussian latitude");

			// mu*mu*F(\lambda,\mu)
			sweet::Data::Sphere2DComplex::DataSpectral h(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid h_phys(i_sphere2DDataConfig);
			h_phys.grid_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function__grid_gaussian(a,b,val);
						c = val;
					}
			);
			h.loadSphere2DDataGrid(h_phys);
			h = ops.mu2(h);

			sweet::Data::Sphere2DComplex::DataSpectral result(i_sphere2DDataConfig);
			sweet::Data::Sphere2DComplex::DataGrid result_phys(i_sphere2DDataConfig);
			result_phys.grid_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &i_data)
					{
						double val;
						testSolutions.test_function__grid_gaussian(lat, mu, val);
						i_data = val;
						i_data *= mu*mu;
					}
			);
			result.loadSphere2DDataGrid(result_phys);


			double div_max_error = (h-result).toGrid().grid_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETErrorFatal(" + ERROR! max error exceeds threshold");
		}
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{

	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	sweet::Data::Sphere2D::Config sphere2DData_Config;
	sphere2DData_Config.setupAuto(shackSphere2DDataOps);
	run_tests(&sphere2DData_Config, shackSphere2DDataOps);

	std::cout << "All test successful" << std::endl;
}

