/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/Data/Cart2DComplex/Cart2DComplex.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>




void setupDataFreq(
		sweet::Data::Cart2D::DataSpectral &io_data,
		int fx,	//!< frequency x
		int fy	//!< frequency y
)
{
	int res_x = io_data.cart2DDataConfig->grid_res[0];
	int res_y = io_data.cart2DDataConfig->grid_res[1];

	// shift by half a cell to generate exactly this mode in spectral space
	double phase_shift = 0.0;

	sweet::Data::Cart2D::DataGrid tmp(io_data.cart2DDataConfig);

	tmp.grid_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = 0.0;

				if (fx >= 0)
					o_data += std::cos(((double)x+phase_shift)*(double)fx*M_PI*2.0/(double)res_x);

				if (fy >= 0)
					o_data += std::cos(((double)y+phase_shift)*(double)fy*M_PI*2.0/(double)res_y);
			}
	);

	io_data.loadCart2DDataGrid(tmp);
}


void setupData123(
		sweet::Data::Cart2D::DataSpectral &io_data
)
{

	sweet::Data::Cart2D::DataGrid tmp(io_data.cart2DDataConfig);

	tmp.grid_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = (x+1.0)+(y+3.0)*y;
			}
	);

	io_data.loadCart2DDataGrid(tmp);
}

int main(int i_argc, char *i_argv[])
{
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Cart2D::Shack *shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	if (shackCart2DDataOps->space_res_spectral[0] <= 0)
	{
		shackCart2DDataOps->space_res_spectral[0] = shackCart2DDataOps->space_res_physical[0];
		shackCart2DDataOps->space_res_spectral[1] = shackCart2DDataOps->space_res_physical[1];
	}
	else
	{
		shackCart2DDataOps->space_res_physical[0] = shackCart2DDataOps->space_res_spectral[0];
		shackCart2DDataOps->space_res_physical[1] = shackCart2DDataOps->space_res_spectral[1];
	}


	int max_res = 64;

	double epsilon = 1e-12;

	for (int res = shackCart2DDataOps->space_res_spectral[0]; res <= max_res; res+=2)
	{
		shackCart2DDataOps->space_res_physical[0] = res;
		shackCart2DDataOps->space_res_physical[1] = res;
		shackCart2DDataOps->space_res_spectral[0] = 0;
		shackCart2DDataOps->space_res_spectral[1] = 0;

		/**
		 * Here we enforce the same physical and spectral resolution
		 */
		sweet::Data::Cart2D::Config cart2DDataConfig;
		cart2DDataConfig.setupAuto(shackCart2DDataOps);
		cart2DDataConfig.printInformation();
		std::cout << std::endl;

		if ((res & 1) == 1)
			SWEETErrorFatal("Only even resolutions supported");

		sweet::Data::Cart2D::DataSpectral a(cart2DDataConfig);

		// relative spectral resolution deltas to test restriction / interpolation
		for (int spec_delta = -4; spec_delta <= 4; spec_delta+=2)
		{
			int dst_res_physical[2] = {
					(int)shackCart2DDataOps->space_res_physical[0]+spec_delta,
					(int)shackCart2DDataOps->space_res_physical[1]+spec_delta
			};

			if (dst_res_physical[0] < 4 || dst_res_physical[1] < 4)
				continue;

			int dst_res_spectral[2] = {0, 0};

			sweet::Data::Cart2D::Config cart2DDataConfigDst;
			cart2DDataConfigDst.setupAuto(dst_res_physical, dst_res_spectral, shackCart2DDataOps->reuse_spectral_transformation_plans);

			/*
			 * Iterate over relative frequencies
			 * Only iterate up to real_modes/2 frequency since higher frequencies would be only
			 * representable as aliased ones.
			 */
			for (int freq_x = 0; freq_x < (int)cart2DDataConfig.spectral_representable_modes[0]; freq_x += 1)
			{
				for (int freq_y = 0; freq_y < (int)cart2DDataConfig.spectral_representable_modes[1]; freq_y += 1)
				{
					std::cout << std::endl;
					std::cout << std::endl;
					std::cout << "*************************************************************" << std::endl;
					std::cout << "Testing (" << res << ", " << res << ") -> (" << cart2DDataConfigDst.spectral_modes[0] << ", " << cart2DDataConfigDst.spectral_modes[1] << ")"  << std::endl;
					std::cout << " + testing frequency in source data: " << freq_x << ", " << freq_y << std::endl;

					int test_freq_x = freq_x;
					if (cart2DDataConfigDst.spectral_representable_modes[0] <= (std::size_t)freq_x)
						test_freq_x = -1;

					int test_freq_y = freq_y;
					if (cart2DDataConfigDst.spectral_representable_modes[1] <= (std::size_t)freq_y)
						test_freq_y = -1;

					std::cout << " + testing frequency in destination data: " << test_freq_x << ", " << test_freq_y << std::endl;
					std::cout << "*************************************************************" << std::endl;

					double error = 0;

					/*
					 * Setup data with highest possible frequency
					 */

					{
						setupDataFreq(a, freq_x, freq_y);

#if 0
						std::cout << "A (physical):" << std::endl;
						a.print_physicalArrayData();
						std::cout << std::endl;

						std::cout << "A (spectral):" << std::endl;
						a.print_spectralData_zeroNumZero();
						std::cout << std::endl;
#endif

						{
							/*
							 * Test for conserving high frequencies in Fourier transformations
							 */
							sweet::Data::Cart2D::DataSpectral tmp = a;

							error = (a-tmp).spectral_reduce_max_abs();
							if (error > epsilon)
							{
								std::cout << "Error: " << error << std::endl;
								SWEETErrorFatal("Test for conserving high frequencies failed! Results should be identical!");
							}
						}

						sweet::Data::Cart2D::DataSpectral b = a.spectral_returnWithDifferentModes(cart2DDataConfigDst);

						{
							sweet::Data::Cart2D::DataSpectral test(cart2DDataConfigDst);
							setupDataFreq(test, test_freq_x, test_freq_y);

#if 0
							std::cout << "Spectral data of b:" << std::endl;
							b.print_spectralData_zeroNumZero();
							std::cout << std::endl;

							std::cout << "Spectral data of test:" << std::endl;
							test.print_spectralData_zeroNumZero();
							std::cout << std::endl;
#endif

							error = (b-test).spectral_reduce_max_abs();

							if (error > epsilon)
							{
								a.print_spectralData_zeroNumZero();
								std::cout << "**************************************************" << std::endl;
								std::cout << "* ERROR" << std::endl;
								std::cout << "* a = freq(fx,fy)" << std::endl;
								std::cout << "* b = interp_restrict(a)" << std::endl;
								std::cout << "* test = expected_interp_restrict(fx,fy)" << std::endl;
								std::cout << "**************************************************" << std::endl;
								std::cout << "Spectral data of A:" << std::endl;
								a.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Grid data of A:" << std::endl;
								a.toGrid().print_physicalData_zeroNumZero();
								std::cout << std::endl;


								std::cout << "Spectral data of b:" << std::endl;
								b.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Grid data of b:" << std::endl;
								b.toGrid().print_physicalData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Spectral data of test:" << std::endl;
								test.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Grid data of test:" << std::endl;
								test.toGrid().print_physicalData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Error: " << error << std::endl;
								SWEETErrorFatal("No modes changed! Results should be identical!");
							}

							std::cout << "PASSED (freq test) with error of " << error << std::endl;
						}
					}

					if (spec_delta >= 0)
					{
						std::cout << "TESTING for conservation of modes" << std::endl;

						setupData123(a);

						sweet::Data::Cart2D::DataSpectral b = a.spectral_returnWithDifferentModes(cart2DDataConfigDst);
						sweet::Data::Cart2D::DataSpectral test = b.spectral_returnWithDifferentModes(cart2DDataConfig);

						error = (a-test).spectral_reduce_max_abs();
						if (error > epsilon)
						{
							std::cout << "**************************************************" << std::endl;
							std::cout << "* ERROR" << std::endl;
							std::cout << "**************************************************" << std::endl;
							std::cout << "Spectral data of A:" << std::endl;
							a.print_spectralData_zeroNumZero();

							std::cout << "Spectral data of b:" << std::endl;
							b.print_spectralData_zeroNumZero();

							std::cout << "Spectral data of test:" << std::endl;
							test.print_spectralData_zeroNumZero();

							std::cout << "Error: " << error << std::endl;
							SWEETErrorFatal("Mode extension requested, but interpolation followed by restriction does not return identical results!");
						}

						std::cout << "PASSED (setup 123) with error of " << error << std::endl;
					}
				}
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
