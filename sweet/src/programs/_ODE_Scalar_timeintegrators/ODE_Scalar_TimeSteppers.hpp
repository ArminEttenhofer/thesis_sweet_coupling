/*
 * ODE_Scalar_TimeSteppers.hpp
 *
 *  Created on: 08 Jun 2022
 * Author: Joao Steinstraessrt <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS__ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TIMESTEPPERS_HPP
#define PROGRAMS__ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TIMESTEPPERS_HPP



#include <programs/_ODE_Scalar_timeintegrators/ODE_Scalar_TS_interface.hpp>
#include <sweet/Shacks/Dictionary.hpp>


class ODE_Scalar_TimeSteppers
{
public:
	ODE_Scalar_TS_interface *master = nullptr;

	ODE_Scalar_TimeSteppers()
	{
	}

	void reset()
	{
		if (master != nullptr)
		{
			delete master;
			master = nullptr;
		}
	}

	void setup(
			//const std::string &i_timestepping_method,
			///int &i_timestepping_order,

			sweet::Shacks::Dictionary &i_shackDict
	)
	{
		reset();
		master = new ODE_Scalar_TS_interface;
		master->setup(
				atof(i_shackDict.user_defined.var[1].c_str()),
				atof(i_shackDict.user_defined.var[2].c_str())
			);
	}

	~ODE_Scalar_TimeSteppers()
	{
		reset();
	}
};


#endif
