/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

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

	std::size_t max_res = 64;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	for (; res_x <= max_res && res_y <= max_res; res_x+=2, res_y+=2)
	{
		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing aliasing pattern with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		shackCart2DDataOps->space_res_physical[0] = res[0];
		shackCart2DDataOps->space_res_physical[1] = res[1];

		sweet::Data::Cart2D::Config cart2DDataConfig;
		cart2DDataConfig.setupAuto(shackCart2DDataOps);

		cart2DDataConfig.printInformation();

		std::cout << "PHYS RES: " << cart2DDataConfig.grid_res[0] << " x " << cart2DDataConfig.grid_res[1] << std::endl;
		std::cout << "PHYS SIZE: " << cart2DDataConfig.grid_data_size[0] << " x " << cart2DDataConfig.grid_data_size[1] << std::endl;
		std::cout << "SPEC MODES: " << cart2DDataConfig.spectral_modes[0] << " x " << cart2DDataConfig.spectral_modes[1] << std::endl;
		std::cout << "SPEC SIZE: " << cart2DDataConfig.spectral_data_size[0] << " x " << cart2DDataConfig.spectral_data_size[1] << std::endl;
		std::cout << std::endl;

#define PRINT_SPECTRUM	1

#if PRINT_SPECTRUM
		std::size_t res_max = 32;
#endif


		sweet::Data::Cart2D::DataSpectral h(cart2DDataConfig);

		h.spectral_setZero();

		h = h.spectral_addScalarAll(1.0);

#if PRINT_SPECTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "All one spectrum:" << std::endl;
			h.print_spectralData_zeroNumZero();
			std::cout << std::endl;
		}
#endif

		for (std::size_t i = 0; i < cart2DDataConfig.spectral_array_data_number_of_elements; i++)
			h.spectral_space_data[i] = {1.0,0.0};


#if PRINT_SPECTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "All one spectrum:" << std::endl;
			h.print_spectralData();
			std::cout << std::endl;
		}
#endif

		h.spectral_zeroAliasingModes();

#if PRINT_SPECTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "Zero aliasing spectrum:" << std::endl;
			h.print_spectralData();
			std::cout << std::endl;
		}
#endif

		h = h.spectral_addScalarAll(1.0);

#if PRINT_SPECTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "Add scalar 1 to non-aliasing spectrum:" << std::endl;
			h.print_spectralData();
			std::cout << std::endl;
			std::cout << "(There should be only 2 and 0)" << std::endl;
		}
#endif

		for (std::size_t i = 0; i < cart2DDataConfig.spectral_array_data_number_of_elements; i++)
		{
			if (h.spectral_space_data[i] != 2.0 && h.spectral_space_data[i] != 0.0)
			{
				SWEETErrorFatal("INCONSISTENT ALIASING !!!");
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
