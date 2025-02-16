/*
 * ShackTimestepControl.hpp
 *
 *  Created on: Feb 21, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_TIMETREE_SHACK_HPP
#define INCLUDE_SWEET_TIMETREE_SHACK_HPP


#include <sweet/Shacks/Base.hpp>
#include <string>
#include <cmath>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>




namespace sweet {
namespace TimeTree {


/*!
 * \brief Shack for TimeTree
 */
class Shack	:
		public sweet::Shacks::Base
{
private:
	static constexpr double SHACK_TIMESTEP_CONTROL_EPSILON = 1e-13;

public:
	/*!
	 * Continue running simulation time stepping.
	 *
	 * This is beneficial to pause simulations if driven interactively.
	 */
	bool runSimulationTimesteps = true;

	//! Number of simulated time steps
	int currentTimestepNr = 0;

	//! Current time step size
	double currentTimestepSize = -1;

	//! Time in simulation
	double currentSimulationTime = 0;

	//! Multiplier for simulation time output
	double outputSimulationTimeMultiplier = 1;

	//! Maximum number of time steps to simulate
	int maxTimestepsNr = -1;

	//! Maximum simulation time to execute the simulation for
	double maxSimulationTime = -1;

	//! Maximum wallclock time to execute the simulation for
	double maxWallclockTime = -1;


	/*!
	 * Check a valid argument for a maximum time step number
	 */
	bool validateMaxTimestepNr()
	{
		if (!(maxTimestepsNr >= 0))
			return error.set("You need to set the maximal number of time steps");

		return true;
	}

	/*!
	 * Check for a valid maximum simulation time
	 */
	bool validateMaxSimulationTime()
	{
		if (!(maxSimulationTime >= 0))
			return error.set("You need to set the maximum simulation time using -t [float]");

		return true;
	}

	/*!
	 * Check for a valid maximum simulation time OR time step number
	 */
	bool validateMaxSimulationTimeOrTimestepNr()
	{
		if (!validateMaxSimulationTime() && !validateMaxTimestepNr())
			return error.set("You need to set the maximum simulation time using -t [float] or time step numbers using -T [int]");

		return true;
	}

	/*!
	 * Check for a valid time step size
	 */
	bool validateTimestepSize()
	{
		if (!(currentTimestepSize > 0))
			return error.set("Timestep size not set, use --dt=[float]");

		return true;
	}

	/*!
	 * Code which we require again and again to be executed before each time step.
	 *
	 * We have a special function for this since we need to keep round-off errors in mind.
	 *
	 * \return true if the current time step size was modified
	 */
	bool timestepHelperStart()
	{
		// If we didn't set max_simulation_time we just continue
		if (maxSimulationTime == -1)
			return false;

		// Check whether there might be some numerical issues of round-off errors in the last time step
		double diff = maxSimulationTime - (currentSimulationTime + currentTimestepSize);

		// Check if we're not close to the maximum simulation time and return
		if (diff > SHACK_TIMESTEP_CONTROL_EPSILON*maxSimulationTime)
			return false;

		/*
		 * This might also change the time step size if it's not necessary,
		 * but we avoid yet another if condition.
		 */
		currentTimestepSize = maxSimulationTime - currentSimulationTime;

		return true;
	}

	/*!
	 * Finish the current time step
	 */
	bool timestepHelperEnd()
	{
		// advance in time
		currentSimulationTime += currentTimestepSize;
		currentTimestepNr++;

		return false;
	}

	/*!
	 * Check whether final time step is reached (already finished).
	 *
	 * This includes various tests such as
	 * - timestep number itself
	 * - maximum simulation time
	 *
	 * There are no rounding errors included here.
	 * The time stepping itself has to care about this.
	 */
	bool isFinalTimestepReached()
	{
		if (maxTimestepsNr >= 0)
		{
			SWEET_ASSERT(maxTimestepsNr >= 0);
			SWEET_ASSERT(currentTimestepNr <= maxTimestepsNr);

			if (maxTimestepsNr == currentTimestepNr)
				return true;
		}

		if (maxSimulationTime >= 0)
		{
			double diff = maxSimulationTime - currentSimulationTime;

			if (diff < 0)
				SWEETErrorFatal("Internal error: This should never happen (diff < 0)");

			if (diff == 0)
			{
				SWEET_ASSERT(maxSimulationTime == currentSimulationTime);
				return true;
			}

#if SWEET_DEBUG
			if (diff < SHACK_TIMESTEP_CONTROL_EPSILON*maxSimulationTime)
				SWEETErrorFatal("Internal error: This should never happen (diff > epsilon)");
#endif
		}

		return false;
	}

	void printProgramArguments(const std::string& i_prefix = "")	override
	{
		std::cout << "" << std::endl;
		std::cout << "Timecontrol:" << std::endl;
		std::cout << "	--dt [float]	timestep size, default=?" << std::endl;
		std::cout << "	--max-wallclock-time [float]	wallclock time limitation, default=-1" << std::endl;
		std::cout << "	-t [float]	maximum simulation time, default=-1 (infinity)" << std::endl;
		std::cout << "	-T [int]	maximum number of time steps, default=-1 (infinity)" << std::endl;
		std::cout << "	-o [float]	time interval at which output should be written, (set to 0 for output at every time step), default=-1 (no output) " << std::endl;
	}


	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)	override
	{
		i_pa.getArgumentValueByKey("--dt", currentTimestepSize);
		i_pa.getArgumentValueByKey("--max-wallclock-time", maxWallclockTime);
		i_pa.getArgumentValueByKey("-t", maxSimulationTime);
		i_pa.getArgumentValueByKey("-T", maxTimestepsNr);

		currentTimestepNr = 0;
		currentSimulationTime = 0;

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(i_pa);
	}

	void printShack(
		const std::string& i_prefix = ""
	)	override
	{
		std::cout << std::endl;
		std::cout << "TIMECONTROL:" << std::endl;
		std::cout << " + runSimulationTimesteps: " << runSimulationTimesteps << std::endl;
		std::cout << " + currentTimestepNr: " << currentTimestepNr << std::endl;
		std::cout << " + currentTimestepSize: " << currentTimestepSize << std::endl;
		std::cout << " + currentSimulationTime: " << currentSimulationTime << std::endl;
		std::cout << " + maxTimestepsNr: " << maxTimestepsNr << std::endl;
		std::cout << " + maxSimulationTime: " << maxSimulationTime << std::endl;
		std::cout << " + maxWallclockTime: " << maxWallclockTime << std::endl;
		std::cout << std::endl;
	}
};

}}

#endif
