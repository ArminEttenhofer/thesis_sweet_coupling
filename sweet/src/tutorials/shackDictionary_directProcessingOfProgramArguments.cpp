/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/IO/Shack.hpp>
#include <iostream>

#include <sweet/Shacks/Dictionary.hpp>

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

	/*
	 * We instantiate a new class PDESWEParametersSphere2D
	 */
	sweet::IO::Shack shackIOData;

	/*
	 * After registering all classes, we can check whether we should output the help information
	 */
	if (pa.argumentWithKeyExists("-h") || pa.argumentWithKeyExists("--help"))
	{
		shackIOData.printProgramArguments();
		return EXIT_FAILURE;
	}

	shackIOData.processProgramArguments(pa);

	std::cout << " + sweParametersSphere2D->printShack()" << std::endl;
	shackIOData.printShack("    ");

	/*
	 * If you activate this, this should trigger an error
	 */
	{
		bool dummy;
		if (!pa.getArgumentValueByKey("-doesntexist", dummy, true))
		{
			std::cout << "Key not found, but this is on purpose" << std::endl;
			pa.error.print();
		}
		else
		{
			std::cerr << "This should have triggered an error, but it didn't. Stopping here." << std::endl;
			return EXIT_FAILURE;
		}
	}

	return 0;
}

