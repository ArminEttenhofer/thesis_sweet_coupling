/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if !SWEET_USE_CART2D_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>

#include <ostream>
#include <cmath>


int main(
		int i_argc,
		char *i_argv[]
)
{
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Cart2D::Shack *shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = shackCart2DDataOps->space_res_physical[0];
	std::size_t res_y = shackCart2DDataOps->space_res_physical[1];

	std::cout << "*************************************************************" << std::endl;
	std::cout << "Testing aliasing pattern with resolution " << res_x << " x " << res_y << std::endl;
	std::cout << "*************************************************************" << std::endl;
	std::size_t res[2] = {res_x, res_y};

	shackCart2DDataOps->space_res_physical[0] = res[0];
	shackCart2DDataOps->space_res_physical[1] = res[1];


	sweet::Data::Cart2D::Config cart2DDataConfigInstance;
	cart2DDataConfigInstance.setupAuto(*shackCart2DDataOps);

	std::size_t test_max_freqx = cart2DDataConfigInstance.spectral_representable_modes[0];
	std::size_t test_max_freqy = cart2DDataConfigInstance.spectral_representable_modes[1];

	std::cout << "*************************************************************" << std::endl;
	cart2DDataConfigInstance.printInformation();
	std::cout << "*************************************************************" << std::endl;
	std::cout << std::endl;

	sweet::Data::Cart2D::DataGrid h1_x(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h2_x(cart2DDataConfigInstance);

	sweet::Data::Cart2D::DataGrid h1_y(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h2_y(cart2DDataConfigInstance);

	sweet::Data::Cart2D::DataGrid h1h2_numerical_x(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h1h2_numerical(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h1h2_numerical_y(cart2DDataConfigInstance);

	sweet::Data::Cart2D::DataGrid h1h2_analytical_x_high(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h1h2_analytical_x_low(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h1h2_analytical_x(cart2DDataConfigInstance);


	sweet::Data::Cart2D::DataGrid h1h2_analytical_y_high(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h1h2_analytical_y_low(cart2DDataConfigInstance);
	sweet::Data::Cart2D::DataGrid h1h2_analytical_y(cart2DDataConfigInstance);

	sweet::Data::Cart2D::DataGrid h1h2_analytical(cart2DDataConfigInstance);

	double epsilon = 1e-10;

	std::size_t max_modes_x = cart2DDataConfigInstance.spectral_representable_modes[0];
	std::size_t max_modes_y = cart2DDataConfigInstance.spectral_representable_modes[1];

	/*
	 * Use cos(x)*cos(x) instead of sin*sin?
	 */
//	int coscos_y = 1;
//	int coscos_x = 1;
	for (int coscos_y = 0; coscos_y < 2; coscos_y++)
	for (int coscos_x = 0; coscos_x < 2; coscos_x++)
	{
		//std::size_t freq1_y = 0;
		for (std::size_t freq1_y = 0; freq1_y < test_max_freqy; freq1_y++)
		{
			/*
			 * Setup frequency
			 * 	sin(fx1*pi*x)
			 */
			h1_y.grid_update_lambda_unit_coordinates_corner_centered(
				[&](double x, double y, double &io_data)
				{
					if (coscos_y)
						io_data = std::cos((double)freq1_y*2.0*M_PI*y);
					else
						io_data = std::sin((double)freq1_y*2.0*M_PI*y);
				}
			);


			//std::size_t freq2_y = 0;
			for (std::size_t freq2_y = 0; freq2_y < test_max_freqy; freq2_y++)
			{
				h2_y.grid_update_lambda_unit_coordinates_corner_centered(
					[&](double x, double y, double &io_data)
					{
						if (coscos_y)
							io_data = std::cos((double)freq2_y*2.0*M_PI*y);
						else
							io_data = std::sin((double)freq2_y*2.0*M_PI*y);
					}
				);

				h1h2_numerical_y = h1_y*h2_y;

				/*
				 * Compute analytical solution
				 *
				 * Setup lower frequency
				 * \frac{1}{2}\cos\left(\left(k_{1}-k_{2}\right)x\pi\right)
				 */
				h1h2_analytical_y_low.grid_update_lambda_unit_coordinates_corner_centered(
						[&](double x, double y, double &io_data)
						{
							io_data = 0.5*std::cos((double)((int)freq1_y-(int)freq2_y)*2.0*M_PI*y);
						}
					);

				if (freq1_y+freq2_y < max_modes_y)
				{
					/*
					 * Setup higher frequency
					 * -\frac{1}{2}\cos\left(\left(k_{1}+k_{2}\right)y\pi\right)
					 */
					h1h2_analytical_y_high.grid_update_lambda_unit_coordinates_corner_centered(
							[&](double x, double y, double &io_data)
							{
								io_data = -0.5*std::cos((double)((int)freq1_y+(int)freq2_y)*2.0*M_PI*y);

								if (coscos_y)
									io_data = -io_data;
							}
						);
				}
				else
				{
					h1h2_analytical_y_high.grid_setValue(0.);
					std::cout << "Higher mode truncated in this case!" << std::endl;
				}


				/*
				 * Merge solutions
				 */
				h1h2_analytical_y = h1h2_analytical_y_low + h1h2_analytical_y_high;

				for (std::size_t freq1_x = 0; freq1_x < test_max_freqx; freq1_x++)
				{
					/*
					 * Setup frequency
					 * 	sin(fx1*pi*x)
					 */
					h1_x.grid_update_lambda_unit_coordinates_corner_centered(
						[&](double x, double y, double &io_data)
						{
							if (coscos_x)
								io_data = std::cos((double)freq1_x*2.0*M_PI*x);
							else
								io_data = std::sin((double)freq1_x*2.0*M_PI*x);
						}
					);


					std::size_t freq2_x = freq1_x;
					//for (std::size_t freq_x2 = 0; freq_x2 < test_max_freqx; freq_x2++)
					{
						std::cout << "***********************************************************" << std::endl;
						std::cout << "Resolution physical space: (" << cart2DDataConfigInstance.grid_res[0] << ", " << cart2DDataConfigInstance.grid_res[1] << ")" << std::endl;
						std::cout << "           spectral space: (" << cart2DDataConfigInstance.spectral_representable_modes[0] << ", " << cart2DDataConfigInstance.spectral_representable_modes[1] << ")" << std::endl;
						std::cout << "Testing for frequency fx1=" << freq1_x << " and fx2=" << freq2_x << ", max frequency x=" << test_max_freqx << std::endl;
						std::cout << "                      fy1=" << freq1_y << " and fy2=" << freq2_y << ", max frequency y=" << test_max_freqy << std::endl;
						std::cout << "                      Test functions: " << (coscos_x ? "coscos_x" : "sinsin_x") << ", " << (coscos_y ? "coscos_y" : "sinsin_y") << std::endl;
						/*
						 * Setup frequency
						 * 	sin(fx1*pi*x)
						 */
						h2_x.grid_update_lambda_unit_coordinates_corner_centered(
							[&](double x, double y, double &io_data)
							{
								if (coscos_x)
									io_data = std::cos((double)freq2_x*2.0*M_PI*x);
								else
									io_data = std::sin((double)freq2_x*2.0*M_PI*x);
							}
						);


						/*
						 * NOTE
						 *
						 * See also
						 * ./doc/software_development_discussions/antialiasing/antialiasing_rule_and_tests.pdf
						 */
						h1h2_numerical_x = h1_x*h2_x;

						/*
						 * Setup lower frequency
						 * \frac{1}{2}\cos\left(\left(k_{1}-k_{2}\right)x\pi\right)
						 */
						h1h2_analytical_x_low.grid_update_lambda_unit_coordinates_cell_centered(
								[&](double x, double y, double &io_data)
								{
									io_data = 0.5*std::cos((double)((int)freq1_x-(int)freq2_x)*2.0*M_PI*x);
								}
							);


						if (freq1_x+freq2_x < max_modes_x)
						{
							/*
							 * Setup higher frequency
							 * -\frac{1}{2}\cos\left(\left(k_{1}+k_{2}\right)x\pi\right)
							 */
							h1h2_analytical_x_high.grid_update_lambda_unit_coordinates_corner_centered(
									[&](double x, double y, double &io_data)
									{
										io_data = -0.5*std::cos((double)((int)freq1_x+(int)freq2_x)*2.0*M_PI*x);

										if (coscos_x)
											io_data = -io_data;
									}
								);
						}
						else
						{
							h1h2_analytical_x_high.grid_setValue(0.);
							std::cout << "Higher mode truncated in this case!" << std::endl;
						}


						h1h2_analytical_x = h1h2_analytical_x_low + h1h2_analytical_x_high;

						{
							h1h2_numerical = h1h2_numerical_x + h1h2_numerical_y;
							h1h2_analytical = h1h2_analytical_x + h1h2_analytical_y;

							double max_error = (h1h2_analytical - h1h2_numerical).grid_reduce_max_abs();
							std::cout << "Max error: " << max_error << std::endl;

							if (max_error > epsilon)
							{
								sweet::Data::Cart2D::DataSpectral tmp_spec(h1_y.cart2DDataConfig);

								std::cout << "H1_y:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1_y);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "H2_y:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h2_y);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "H1_x:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1_x);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "H2_x:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h2_x);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Numerical spectrum:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1h2_numerical);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum x:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1h2_analytical_x);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum y (low):" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1h2_analytical_y_low);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum y (high):" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1h2_analytical_y_high);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum y:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1h2_analytical_y);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum:" << std::endl;
								tmp_spec.loadCart2DDataGrid(h1h2_analytical);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								if (max_error > epsilon)
									SWEETErrorFatal("Error too high!");
							}
						}
					}
				}
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
