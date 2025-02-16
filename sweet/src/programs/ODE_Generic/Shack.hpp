/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_ODE_GENERIC_SHACK_HPP
#define PROGRAMS_ODE_GENERIC_SHACK_HPP

#include <string>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Base.hpp>

namespace ODE_Generic {

class Shack	:
		public sweet::Shacks::Base
{
public:
	//! String of time stepping method
	std::string timestepping_method;

	//! ODE string
	std::string ode = "";

	/*
	 * Compute errors compared to analytical solution
	 */
	bool compute_errors = false;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "ShackODEGeneric parameters:" << std::endl;
		std::cout << i_prefix << "	--ode [str]	string of ODE" << std::endl;
		std::cout << i_prefix << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << i_prefix << "	--timestepping-method [string] String of time stepping method" << std::endl;
		std::cout << i_prefix << "	                               Use 'help' or 'helpall' to get more information" << std::endl;
	}


	void printShack(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "ShackODEGeneric parameters:" << std::endl;
		std::cout << i_prefix << " + ode: " << ode << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << i_prefix << " + timestepping_method: " << timestepping_method << std::endl;
	}


	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--ode", ode);
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);
		i_pa.getArgumentValueByKey("--timestepping-method", timestepping_method);
		i_pa.getArgumentValueByKey("--tm", timestepping_method);
		
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_pa);
		return true;
	}
};

}

#endif
