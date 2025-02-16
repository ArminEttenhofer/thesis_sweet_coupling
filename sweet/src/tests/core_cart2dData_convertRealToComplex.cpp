/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *      
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/Error/Base.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>


#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Cart2DComplex_DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataSpectral_2_Cart2D_DataSpectral.hpp>
#include <sweet/Data/Cart2D/Convert/DataGrid_2_Cart2DComlex_DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataGrid_2_Cart2D_DataGrid.hpp>


class TestCart2DDataModes
{
public:
	sweet::Error::Base error;


	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::Error::Base error;

		sweet::Data::Cart2D::Config cart2DDataConfig;
		sweet::Data::Cart2D::Operators ops;

		sweet::Data::Cart2D::DataSpectral prog_h;


		bool setup(sweet::Data::Cart2D::Shack *i_shackCart2DDataOps)
		{
			/*
			 * Setup Cart2D Data Config & Operators
			 */
			cart2DDataConfig.setupAuto(*i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);

			ops.setup(cart2DDataConfig, *i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			prog_h.setup(cart2DDataConfig);

			return true;
		}

		void clear()
		{
			prog_h.clear();

			ops.clear();
			cart2DDataConfig.clear();
		}
	};

	// Simulation data
	Data data;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;

public:
	TestCart2DDataModes(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackCart2DDataOps(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}


	bool setup()
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.setup();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.printShackData();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		data.setup(shackCart2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(data);

		return true;
	}


	void clear()
	{
		data.clear();

		shackCart2DDataOps = nullptr;
		shackProgArgDict.clear();
	}

	void test_cart2ddata_cart2ddatacomplex_physicalgrid_convert()
	{
		sweet::Data::Cart2D::DataGrid test(data.cart2DDataConfig);
		sweet::Data::Cart2DComplex::DataGrid testcplx(data.cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral test_spc(data.cart2DDataConfig);
		sweet::Data::Cart2DComplex::DataSpectral test_spc_cplx(data.cart2DDataConfig);

		for (std::size_t y = 0; y < data.cart2DDataConfig.grid_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.cart2DDataConfig.grid_data_size[0]; x++)
			{
				std::cout << "test_cart2ddata_cart2ddatacomplex_physicalgrid_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test.grid_setZero();
				test.grid_setValue(y, x, 1.0);

				testcplx = sweet::Data::Cart2D::Convert::DataGrid_2_Cart2DComplex_DataGrid::convert(test);

				test_spc_cplx.loadCart2DDataGrid(testcplx);

				test_spc_cplx.test_realphysical();
				sweet::Data::Cart2D::DataGrid tmp = sweet::Data::Cart2DComplex::Convert::DataGrid_2_Cart2D_DataGrid::convert_real(testcplx);

				double error = (test-tmp).grid_reduce_max_abs();

				if (error > 1e-8)
				{
					test_spc.loadCart2DDataGrid(test);
					sweet::Data::Cart2D::DataSpectral tmp_spc(tmp.cart2DDataConfig);
					tmp_spc.loadCart2DDataGrid(tmp);

					std::cout << std::endl;
					std::cout << "Original" << std::endl;
					test_spc.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original complex" << std::endl;
					test_spc_cplx.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original converted" << std::endl;
					tmp_spc.print_spectralData_zeroNumZero();

					SWEETErrorFatal("Inconsistency detected a");
				}
			}
		}


		for (std::size_t y = 0; y < data.cart2DDataConfig.spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.cart2DDataConfig.spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_physicalgrid_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_setZero();
					test_spc.spectral_set(y, x, 1.0);
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toGrid();

					testcplx = sweet::Data::Cart2D::Convert::DataGrid_2_Cart2DComplex_DataGrid::convert(test);
					test_spc_cplx.loadCart2DDataGrid(testcplx);
					test_spc_cplx.test_realphysical();

					sweet::Data::Cart2D::DataGrid tmp = sweet::Data::Cart2DComplex::Convert::DataGrid_2_Cart2D_DataGrid::convert_real(testcplx);

					double error = (test-tmp).grid_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadCart2DDataGrid(test);
						sweet::Data::Cart2D::DataSpectral tmp_spc(tmp.cart2DDataConfig);
						tmp_spc.loadCart2DDataGrid(tmp);

						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test_spc.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						test_spc_cplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp_spc.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected b");
					}
				}

				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_physicalgrid_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_setZero();
					test_spc.spectral_set(y, x, std::complex<double>(0.0, 1.0));
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toGrid();

					testcplx = sweet::Data::Cart2D::Convert::DataGrid_2_Cart2DComplex_DataGrid::convert(test);
					test_spc_cplx.loadCart2DDataGrid(testcplx);
					test_spc_cplx.test_realphysical();

					sweet::Data::Cart2D::DataGrid tmp = sweet::Data::Cart2DComplex::Convert::DataGrid_2_Cart2D_DataGrid::convert_real(testcplx);

					double error = (test-tmp).grid_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadCart2DDataGrid(test);
						sweet::Data::Cart2D::DataSpectral tmp_spc(tmp.cart2DDataConfig);
						tmp_spc.loadCart2DDataGrid(tmp);

						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test_spc.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						test_spc_cplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp_spc.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected c");
					}
				}

				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_physicalgrid_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_setZero();
					test_spc.spectral_set(y, x, std::complex<double>(1.0, 1.0));
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toGrid();

					testcplx = sweet::Data::Cart2D::Convert::DataGrid_2_Cart2DComplex_DataGrid::convert(test);
					test_spc_cplx.loadCart2DDataGrid(testcplx);
					test_spc_cplx.test_realphysical();

					sweet::Data::Cart2D::DataGrid tmp = sweet::Data::Cart2DComplex::Convert::DataGrid_2_Cart2D_DataGrid::convert_real(testcplx);

					double error = (test-tmp).grid_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadCart2DDataGrid(test);
						sweet::Data::Cart2D::DataSpectral tmp_spc(tmp.cart2DDataConfig);
						tmp_spc.loadCart2DDataGrid(tmp);

						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test_spc.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						test_spc_cplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp_spc.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected d");
					}
				}
			}
		}
	}


	void test_cart2ddata_cart2ddatacomplex_physicalspectral_convert()
	{
		sweet::Data::Cart2D::DataSpectral test(data.cart2DDataConfig);
		sweet::Data::Cart2DComplex::DataSpectral testcplx(data.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid test_phys(data.cart2DDataConfig);
		sweet::Data::Cart2DComplex::DataGrid test_phys_cplx(data.cart2DDataConfig);

		for (std::size_t y = 0; y < data.cart2DDataConfig.grid_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.cart2DDataConfig.grid_data_size[0]; x++)
			{
				std::cout << "test_cart2ddata_cart2ddatacomplex_physicalspectral_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test_phys.grid_setZero();
				test_phys.grid_setValue(y, x, 1.0);

				test.loadCart2DDataGrid(test_phys);

				//testcplx = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(test);
				testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(test);
				testcplx.test_realphysical();

				sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

				double error = (test-tmp).spectral_reduce_max_abs();

				if (error > 1e-8)
				{
					std::cout << std::endl;
					std::cout << "Original" << std::endl;
					test.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original complex" << std::endl;
					testcplx.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original converted" << std::endl;
					tmp.print_spectralData_zeroNumZero();

					SWEETErrorFatal("Inconsistency detected e");
				}
			}
		}


		for (std::size_t y = 0; y < data.cart2DDataConfig.spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.cart2DDataConfig.spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_physicalspectral_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_setZero();
					test.spectral_set(y, x, 1.0);
					test.spectral_zeroAliasingModes();

					test.loadCart2DDataGrid(test.toGrid()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(test);
					testcplx.test_realphysical();

					sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected f");
					}
				}

				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_physicalspectral_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_setZero();
					test.spectral_set(y, x, std::complex<double>(0.0, 1.0));
					test.spectral_zeroAliasingModes();

					test.loadCart2DDataGrid(test.toGrid()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(test);
					testcplx.test_realphysical();

					sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected g");
					}
				}

				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_physicalspectral_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_setZero();
					test.spectral_set(y, x, std::complex<double>(1.0, 1.0));
					test.spectral_zeroAliasingModes();

					test.loadCart2DDataGrid(test.toGrid()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(test);
					testcplx.test_realphysical();

					sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected h");
					}
				}
			}
		}
	}

	void test_cart2ddata_cart2ddatacomplex_spectralgrid_convert()
	{
		sweet::Data::Cart2D::DataSpectral test(data.cart2DDataConfig);
		sweet::Data::Cart2DComplex::DataSpectral testcplx(data.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid test_phys(data.cart2DDataConfig);
		sweet::Data::Cart2DComplex::DataGrid test_phys_cplx(data.cart2DDataConfig);


		for (std::size_t y = 0; y < data.cart2DDataConfig.grid_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.cart2DDataConfig.grid_data_size[0]; x++)
			{
				std::cout << "test_cart2ddata_cart2ddatacomplex_spectralgrid_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test_phys.grid_setZero();
				test_phys.grid_setValue(y, x, 1.0);
				test.loadCart2DDataGrid(test_phys);

				testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::spectral_convert(test);
				testcplx.test_realphysical();

				sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

				double error = (test-tmp).spectral_reduce_max_abs();

				if (error > 1e-8)
				{
					std::cout << std::endl;
					std::cout << "Original" << std::endl;
					test.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original complex" << std::endl;
					testcplx.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original converted" << std::endl;
					tmp.print_spectralData_zeroNumZero();

					SWEETErrorFatal("Inconsistency detected i");
				}
			}
		}


		for (std::size_t y = 0; y < data.cart2DDataConfig.spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.cart2DDataConfig.spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_spectralgrid_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_setZero();
					test.spectral_set(y, x, 1.0);
					test.spectral_zeroAliasingModes();

					testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::spectral_convert(test);
					testcplx.test_realphysical();

					sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected j");
					}
				}

				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_spectralgrid_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_setZero();
					test.spectral_set(y, x, std::complex<double>(0.0, 1.0));
					test.spectral_zeroAliasingModes();

					testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::spectral_convert(test);
					testcplx.test_realphysical();

					sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);
					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected k");
					}
				}

				{
					std::cout << "test_cart2ddata_cart2ddatacomplex_spectralgrid_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_setZero();
					test.spectral_set(y, x, std::complex<double>(1.0, 1.0));
					test.spectral_zeroAliasingModes();

					testcplx = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::spectral_convert(test);
					testcplx.test_realphysical();

					sweet::Data::Cart2D::DataSpectral tmp = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_grid_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETErrorFatal("Inconsistency detected l");
					}
				}
			}
		}
	}

	void run_tests()
	{
		test_cart2ddata_cart2ddatacomplex_physicalgrid_convert();
		test_cart2ddata_cart2ddatacomplex_physicalspectral_convert();
		test_cart2ddata_cart2ddatacomplex_spectralgrid_convert();
	}
};



int main(int i_argc, char *i_argv[])
{
	TestCart2DDataModes simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.run_tests();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.clear();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);



	std::cout << "FIN" << std::endl;
	return 0;
}
