/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_SHACKPDEADVECTIONSPHERE2D_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_SHACKPDEADVECTIONSPHERE2D_HPP


#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>


/**
 * Coefficients for the PDE Advection equation on the sphere2D
 */
class ShackPDEAdvectionSphere2D	:
		public sweet::Shacks::Base
{
public:
	/*
	 * Compute errors compared to analytical solution
	 */
	bool compute_errors = false;


	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << "PDE advection sphere2D:" << std::endl;
		std::cout << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);

		return error.forwardWithPositiveReturn(i_pa.error);
	}


	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << i_prefix << "PDE advection sphere2D:" << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << std::endl;
	}
};


#endif
