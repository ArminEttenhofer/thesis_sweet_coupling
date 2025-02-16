/*
 * Author: Pedro Peixoto <ppeixoto@usp.br>
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 * MULE_SCONS_OPTIONS: --eigen=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>

#include "PDE_SWECart2D_normalModeConversion/SWECart2DNormalModes.hpp"

#include <ostream>
#include <cmath>


typedef std::complex<double> complex;


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

	SWECart2DNormalModes sweCart2DNormalModes;
	sweCart2DNormalModes.shackRegistration(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	sweet::Data::Cart2D::Config cart2DDataConfig;
	cart2DDataConfig.setupAuto(shackCart2DDataOps);

	sweet::Data::Cart2D::Operators ops;
	ops.setup(cart2DDataConfig, shackCart2DDataOps);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(ops);



	double eps = 1e-9;

	std::cout << "*************************************************************" << std::endl;
	std::cout << "Testing normal mode conversion tool" << std::endl;
	std::cout << "*************************************************************" << std::endl;


	std::cout << "*************************************************************" << std::endl;
	cart2DDataConfig.printInformation();
	std::cout << "*************************************************************" << std::endl;
	std::cout << std::endl;

	sweet::Data::Cart2D::DataSpectral h(cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral u(cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral v(cart2DDataConfig);

	//Normal modes to be generated
	double geo_mode=1.0;
	double igwest_mode=1.0;
	double igeast_mode=0.0;
	
	h.spectral_setZero();
	u.spectral_setZero();
	v.spectral_setZero();

	std::cout << "Adding normal mode with coefficients:" << std::endl;
	std::cout << " Geost: " << geo_mode<<std::endl;
	std::cout << " IGWest: " << igwest_mode<<std::endl;
	std::cout << " IGEast: " << igeast_mode<<std::endl;
	std::cout << "To wavenumber  ";

	for (std::size_t ik1 = 0; ik1 < cart2DDataConfig.spectral_data_size[1]/4; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < cart2DDataConfig.spectral_data_size[0]/2; ik0++)
		{
			std::cout << " (" << ik0 << "," << ik1<< ") , ";

			sweCart2DNormalModes.add_normal_mode(
									ik0, ik1,
									geo_mode,
									igwest_mode,
									igeast_mode,
									h,
									u,
									v
							);
		}
	}
	
	std::cout << "\n\n Extracting Normal Modes:" << std::endl;

	sweet::Data::Cart2D::DataSpectral geo(cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral igwest(cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral igeast(cart2DDataConfig);
	complex geo_mode_c;
	complex igwest_mode_c;
	complex igeast_mode_c;

	sweCart2DNormalModes.convert_allspectralmodes_2_normalmodes(
		h,	
		u,
		v,
		geo,
		igwest,
		igeast
	);

	for (std::size_t ik1 = 0; ik1 < cart2DDataConfig.spectral_data_size[1]/4; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < cart2DDataConfig.spectral_data_size[0]/2; ik0++)
		{
			std::cout << "From wavenumber (" << ik0 << "," << ik1<< "):";
			geo_mode_c = geo.spectral_get(ik1, ik0);
			igwest_mode_c = igwest.spectral_get(ik1, ik0);
			igeast_mode_c = igeast.spectral_get(ik1, ik0);

			std::cout << " Geost: " << geo_mode_c;
			std::cout << " IGWest: " << igwest_mode_c;
			std::cout << " IGEast: " << igeast_mode_c;

			double error = abs(geo_mode_c.real()-geo_mode)+abs(igwest_mode_c.real()-igwest_mode)+abs(igeast_mode_c.real()-igeast_mode);
			std::cout << " Error: " << error ;
			if (error<10e-10)
			{
				std::cout << " PASSED " << std::endl;
			}
			else //mode may have been split to mirror
			{
				error = abs(geo_mode_c.real()-geo_mode/2.0)+abs(igwest_mode_c.real()-igwest_mode/2.0)+abs(igeast_mode_c.real()-igeast_mode/2.0);
				if (error < eps)
				{
					std::cout << " PASSED " << std::endl;
				}
				else
				{
					std::cout << " FAIL " << std::endl;
					std::cerr << "NORMAL_MODES: Failed to convert back and forward to/from normal modes!" << error << std::endl;
					exit(-1);
				}
			}
		}
	}	

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
