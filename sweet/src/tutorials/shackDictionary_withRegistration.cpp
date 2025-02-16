/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <iostream>

#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>


int main(int i_argc, char *i_argv[])
{
	/*
	 * We start by setting up the class to parse the program arguments
	 */
	std::cout << " + ProgramArguments()" << std::endl;
	sweet::Tools::ProgramArguments pa;
	if (!pa.setup(i_argc, i_argv))
	{
		std::cout << "Error: " << pa.error.get() << std::endl;
		return 1;
	}


	{
		/*
		 * Now we warmup the dictionary which can be stored to store arbitrary
		 * classes.
		 *
		 * This is specialized for processing our program arguments and using it
		 * later on in the SWEET programs.
		 */
		std::cout << " + ShackDictionary()" << std::endl;
		sweet::Shacks::Dictionary varShackDict;


		/*
		 * Register new classes
		 */
		std::cout << "   + registerFirstTime<...>()" << std::endl;
		varShackDict.registerFirstTime<sweet::IO::Shack>();
		varShackDict.registerFirstTime<sweet::TimeTree::Shack>();

		/*
		 * Now we close the registration
		 *
		 * This will avoid performance bugs!
		 */
		varShackDict.closeRegistration();

		{
			/*
			 * If we now try to register a new class, this should raise an error!
			 */
			bool retval = varShackDict.registerFirstTime<sweet::TimeTree::Shack>();

			if (!retval)
			{
				// Just get the error (deleting it) and continue
				varShackDict.error.get();
			}
			else
			{
				std::cerr << "varShackDict.registerFirstTime<> should have raised an error!" << std::endl;
				return EXIT_FAILURE;
			}
		}

		/*
		 * After registering all classes, we can check whether we should output the help information
		 */
		if (pa.argumentWithKeyExists("-h") || pa.argumentWithKeyExists("--help"))
		{
			varShackDict.printProgramArguments();
			return EXIT_FAILURE;
		}

		/*
		 * Now its time to process all program arguments with all registered classes
		 */
		varShackDict.processProgramArguments(pa);

		/*
		 * Get handler to new class ShackIOData
		 */
		sweet::IO::Shack *shackIOData = varShackDict.get<sweet::IO::Shack>();
		if (shackIOData == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varShackDict.error.get() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Get handler to new class ShackIOData
		 */
		sweet::TimeTree::Shack *shackTimeTree = varShackDict.get<sweet::TimeTree::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(varShackDict);

		/*
		 * Now we close getting the parameter class
		 *
		 * This will avoid performance bugs!
		 */
		varShackDict.closeGet();

		{
			/*
			 * If we now access a class instance, this should raise an error!
			 */

			sweet::TimeTree::Shack *ioDataParameters = varShackDict.get<sweet::TimeTree::Shack>();
			if (ioDataParameters == nullptr)
			{
				// Just get the error (deleting it) and continue
				varShackDict.error.get();
			}
			else
			{
				std::cerr << "This should have raised an error!" << std::endl;
				return EXIT_FAILURE;
			}
		}

		/*
		 * Now its time to output all program arguments of all registered classes
		 */
		std::cout << " + varShackDict.printShack()" << std::endl;
		varShackDict.printShackData("    ");

		/*
		 * And we can also print them individually
		 */
		std::cout << " + shackIOData->printShack()" << std::endl;
		shackIOData->printShack("    ");

		std::cout << " + shackParallelization->printShack()" << std::endl;
		shackTimeTree->printShack("    ");

		/*
		 * Of course, we can also access them directly
		 */
		std::cout << " + direct access:" << std::endl;
		std::cout << "   output_file_name: " << shackIOData->outputFileName << std::endl;
	}

	return 0;
}

