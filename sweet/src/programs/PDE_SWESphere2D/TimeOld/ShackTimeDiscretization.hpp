/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_SHACKTIMEDISCRETIZATION_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_SHACKTIMEDISCRETIZATION_HPP



#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>


class ShackTimeDiscretization	:
		public sweet::Shacks::Base
{
public:
	//! String of time stepping method
	std::string timestepping_method;

	//! Order of time stepping
	int timestepping_order = -1;

	//! Order of 2nd time stepping which might be used
	int timestepping_order2 = -1;

	//! Crank-Nicolson filter
	double timestepping_crank_nicolson_filter = 0.5;

	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << "Time Discretization:" << std::endl;
		std::cout << i_prefix << "	--timestepping-method [string] String of time stepping method" << std::endl;
		std::cout << i_prefix << "	                               Use 'help' or 'helpall' to get more information" << std::endl;
		std::cout << "	--timestepping-order [int]			Specify the order of the time stepping" << std::endl;
		std::cout << "	--timestepping-order2 [int]			Specify the order of the time stepping" << std::endl;
		std::cout << "	--crank-nicolson-filter [float]		Crank-Nicolson filter (default=0.5)" << std::endl;

	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--timestepping-method", timestepping_method);
		i_pa.getArgumentValueByKey("--tm", timestepping_method);
		i_pa.getArgumentValueByKey("-R", timestepping_order);
		i_pa.getArgumentValueByKey("--timestepping-order", timestepping_order);
		i_pa.getArgumentValueByKey("--timestepping-order2", timestepping_order2);
		i_pa.getArgumentValueByKey("--crank-nicolson-filter", timestepping_crank_nicolson_filter);

		if (i_pa.error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		return true;

	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << "TIME DISCRETIZATION:" << std::endl;
		std::cout << " + timestepping_method: " << timestepping_method << std::endl;
		std::cout << " + timestepping_order: " << timestepping_order << std::endl;
		std::cout << " + timestepping_order2: " << timestepping_order2 << std::endl;
		std::cout << " + timestepping_crank_nicolson_filter: " << timestepping_crank_nicolson_filter << std::endl;
		std::cout << std::endl;
	}
};




#endif
