/*
 * ShackIOData.hpp
 *
 *  Created on: Feb 21, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_IO_SHACK_HPP
#define INCLUDE_SWEET_IO_SHACK_HPP


#include <getopt.h>
#include <iomanip>
#include <string>
#include <iostream>
#include <limits>

#include <sweet/Shacks/Base.hpp>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#include "../TimeTree/Shack.hpp"



namespace sweet {
namespace IO {


/*!
 * \brief Shack for Input and output data
 */
class Shack	:
		public sweet::Shacks::Base
{
public:

	//! prefix of filename for outputConfig of data
public:
	std::string outputFileName = "X";

	//! output mode of variables
public:
	std::string outputFileMode = "default";

	//! output each of the given simulation seconds
public:
	double outputEachSimTime = -1;

	//! Simulation seconds for next outputConfig
public:
	double _outputNextSimTime = 0;

	//! Last timestep number for output (to avoid duplicates)
private:
	int _outputLastSimTimestepNr = -1;

	//! time scaling for outputConfig
	//! e.g. use scaling by 1.0/(60*60) to output days instead of seconds
public:
	double outputFormatTimeScale = 1.0;

	//! inverse time scaling for outputConfig
	//! If this is set, the output_time_scale is overwritten with its inverse
public:
	double outputFormatTimeScaleInv = 1.0;

	//! precision for floating point outputConfig to std::cout and std::endl
public:
	int outputFormatFPPrecision = std::numeric_limits<double>::digits10 + 1;


	//! set verbosity of simulation
public:
	int verbosity = 0;

	//! activate GUI mode?
public:
	bool guiEnabled = (SWEET_GUI == 0 ? false : true);


	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << std::endl;
		std::cout << "IOData:" << std::endl;
		std::cout << "	--output-file-name [string]		String specifying the name of the output file" << std::endl;
		std::cout << "	--output-file-mode [string]		Format of output file, default: default" << std::endl;
		std::cout << "	--output-time-scale [float]		Output time scale, default: 1" << std::endl;
		std::cout << "	--output-time-scale-inv [float]		Inverse of output time scale, default: 1" << std::endl;
		std::cout << "	-v [int]			verbosity level" << std::endl;
		std::cout << "	-G [0/1]			graphical user interface" << std::endl;

		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--output-file-mode", outputFileMode);
		i_pa.getArgumentValueBy2Keys("--output-file-name", "-O", outputFileName);

		if (i_pa.getArgumentValueByKey("--output-time-scale", outputFormatTimeScale))
			outputFormatTimeScaleInv = 1.0/outputFormatTimeScale;

		if (i_pa.getArgumentValueByKey("--output-time-scale-inv", outputFormatTimeScaleInv))
			outputFormatTimeScale = 1.0/outputFormatTimeScaleInv;

		i_pa.getArgumentValueByKey("-d", outputFormatFPPrecision);
		i_pa.getArgumentValueByKey("-o", outputEachSimTime);

		if (i_pa.error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		if (outputFileMode == "default")
		{
			outputFileMode = "bin";

			if (outputFileName == "X")
				outputFileName = "output_%s_t%020.8f.sweet";
		}
		else
		{
			if (outputFileName == "X")
			{
				if (outputFileMode == "csv")
					outputFileName = "output_%s_t%020.8f.csv";
				else if (outputFileMode == "bin")
					outputFileName = "output_%s_t%020.8f.sweet";
				else if (outputFileMode == "csv_spec_evol")
					outputFileName = "output_%s_t%020.8f.txt";
				else
					return error.set("Unknown filemode '"+outputFileMode+"'");
			}
		}

		if (outputFileName == "-")
			outputFileName = "";

		if (outputFormatFPPrecision >= 0)
		{
			std::cout << std::setprecision(outputFormatFPPrecision);
			std::cerr << std::setprecision(outputFormatFPPrecision);
		}

		i_pa.getArgumentValueByKey("-G", guiEnabled);
		i_pa.getArgumentValueByKey("-v", verbosity);

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(i_pa);
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << "INPUT/OUTPUT:" << std::endl;
		std::cout << " + output_file_name: " << outputFileName << std::endl;
		std::cout << " + output_file_mode: " << outputFileMode << std::endl;
		std::cout << " + output_each_sim_seconds: " << outputEachSimTime << std::endl;
		std::cout << " + output_next_sim_seconds: " << _outputNextSimTime << std::endl;
		std::cout << " + output_time_scale: " << outputFormatTimeScale << std::endl;
		std::cout << " + output_time_scale_inv: " << outputFormatTimeScaleInv << std::endl;
		std::cout << " + output_floating_point_precision: " << outputFormatFPPrecision << std::endl;

		std::cout << " + verbosity: " << verbosity << std::endl;
		std::cout << " + gui_enabled: " << guiEnabled << std::endl;

		std::cout << std::endl;
	}

	/*!
	 * Check whether some output should be done.
	 *
	 * Output can be requested
	 *  - each simulation interval
	 *  - (...) more features coming soon
	 *
	 * \return true if an output should be done
	 */
	bool outputHelper(
			const sweet::TimeTree::Shack *i_shackTimeTree
	)
	{
		double eps = 1e-5;
		if (outputEachSimTime >= 0)
		{
			if (_outputLastSimTimestepNr != i_shackTimeTree->currentTimestepNr)
			{
				/*
				 * We could suffer of some roundoff errors.
				 * Hence, we check whehter we are within a certain interval given the time step size.
				 *
				 * curtime - 0.5*dt < t < curtime + 0.5*dt
				 */
				if (_outputNextSimTime < i_shackTimeTree->currentSimulationTime + i_shackTimeTree->currentTimestepSize*eps)
				{

					_outputNextSimTime += outputEachSimTime;

					if (_outputNextSimTime < i_shackTimeTree->currentSimulationTime)
					{
						int fac = (i_shackTimeTree->currentSimulationTime - _outputNextSimTime)/outputEachSimTime;
						_outputNextSimTime += fac*outputEachSimTime;

						if (_outputNextSimTime < i_shackTimeTree->currentSimulationTime)
							_outputNextSimTime += outputEachSimTime;

					}

					// Make sure that we're not beyond the maximum simulation time
					if (_outputNextSimTime > i_shackTimeTree->maxSimulationTime)
						_outputNextSimTime = i_shackTimeTree->maxSimulationTime;

					// Backup last timestep number to avoid duplicate outputs
					_outputLastSimTimestepNr = i_shackTimeTree->currentTimestepNr;

					if (i_shackTimeTree->maxSimulationTime == outputEachSimTime)
					{
						/*
						 * It seems that only the final time step should be printed
						 * => Check if we're at timestamp=0 and if it is so, then skip
						 */
						if (i_shackTimeTree->currentSimulationTime == 0)
							return false;
					}

					return true;
				}
			}
		}

		return false;
	}


	// DEPRECATED, switch to outputHelperShouldOutput
	__attribute__((deprecated))
	bool checkDoOutput(
			double i_current_simulation_time
	)
	{
		// output each time step
		if (outputEachSimTime < 0)
			return false;

		if (_outputNextSimTime-_outputNextSimTime*(1e-12) > i_current_simulation_time)
			return false;

		return true;
	}

	// DEPRECATED, switch to outputHelperShouldOutput
	__attribute__((deprecated))
	void advanceNextOutput(
			double i_current_simulation_time,
			double i_max_simulation_time
	)
	{
		if (_outputNextSimTime == i_max_simulation_time)
		{
			_outputNextSimTime = std::numeric_limits<double>::infinity();
		}
		else
		{
			while (_outputNextSimTime-_outputNextSimTime*(1e-12) <= i_current_simulation_time)
				_outputNextSimTime += outputEachSimTime;

			if (_outputNextSimTime > i_max_simulation_time)
				_outputNextSimTime = i_max_simulation_time;
		}

		return;
	}
};

}}

#endif
