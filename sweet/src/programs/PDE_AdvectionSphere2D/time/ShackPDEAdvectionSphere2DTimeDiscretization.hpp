/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_SHACKPDEADVECTIONSPHERE2DTIMEDISCRETIZATION_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_SHACKPDEADVECTIONSPHERE2DTIMEDISCRETIZATION_HPP


#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>


class ShackPDEAdvectionSphere2DTimeDiscretization	:
		public sweet::Shacks::Base
{
public:
	//! String of time stepping method
	std::string timestepping_method;

	//! Order of time stepping
	int timestepping_order = -1;

	//! Order of 2nd time stepping which might be used
	int timestepping_order2 = -1;

	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << "Time Discretization:" << std::endl;
		std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;
		std::cout << "	--timestepping-order [int]			Specify the order of the time stepping" << std::endl;
		std::cout << "	--timestepping-order2 [int]			Specify the order of the time stepping" << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--timestepping-method", timestepping_method);
		i_pa.getArgumentValueByKey("-R", timestepping_order);
		i_pa.getArgumentValueByKey("--timestepping-order", timestepping_order);
		i_pa.getArgumentValueByKey("--timestepping-order2", timestepping_order2);

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
		std::cout << std::endl;
	}
};




#endif
