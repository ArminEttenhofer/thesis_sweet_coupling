/*
 * ODEGeneric_DE_Dahlquist_Main.hpp
 *
 *  Created on: Jun 13, 2023
 *      Author: martin
 */

#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DE_DAHLQUIST_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DE_DAHLQUIST_HPP

#include <iostream>
#include <vector>
#include <cstdio>

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Data/GenericContainer/ConfigBase.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/IO/Shack.hpp>

#include "Benchmarks/BenchmarkRegistry.hpp"
#include "TimeTree/TimeTreeIR.hpp"
#include "DataContainer/Simulation.hpp"
#include "../DE_Base.hpp"
#include "DataContainer/Config.hpp"
#include "FileOutput.hpp"


namespace ODE_Generic {
namespace DE_Dahlquist {

class DE_Dahlquist :
		public DE_Base,
		public sweet::Data::GenericContainer::CastHelper<DataContainer::Simulation, DataContainer::Config>
{
	DataContainer::Config *config;

	// time integrators
	TimeTree *timeTree;

	sweet::IO::Shack *shackIO;
	sweet::TimeTree::Shack *shackTimeTree;

	// Handler to all benchmarks
	Benchmarks::BenchmarkRegistry *benchmarks;

	FileOutput fileOutput;

public:
	DE_Dahlquist()	:
		config(nullptr),
		timeTree(nullptr),
		benchmarks(nullptr)
	{
	}


	bool setup(
			sweet::Shacks::ProgramArgumentsDictionary *i_progArgShackDict,
			const std::string &i_timesteppingMethod
	)	override
	{
		SWEET_ASSERT(config == nullptr);
		config = new DataContainer::Config;

		/*
		 * Time Tree
		 */
		timeTree = new TimeTree;

#if !SWEET_XBRAID
		timeTree->setup_1_registerAllTimesteppers();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeTree);

		timeTree->setup_2_shackRegistration(i_progArgShackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeTree);

		i_progArgShackDict->processProgramArguments();

		if (i_timesteppingMethod == "")
			return error.set("Use --timestepping-method=... to choose time stepper");

		timeTree->setup_3_timestepper(
			i_timesteppingMethod,
			config
		);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeTree);
#endif

		/*
		 * Benchmarks
		 */
		benchmarks = new Benchmarks::BenchmarkRegistry;

		benchmarks->setup_1_registerAllBenchmark();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*benchmarks);

		benchmarks->setup_2_shackRegistration(i_progArgShackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*benchmarks);

		i_progArgShackDict->processProgramArguments();

		benchmarks->setup_3_benchmarkDetection();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*benchmarks);

		benchmarks->setup_4_benchmarkSetup();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*benchmarks);

		i_progArgShackDict->processProgramArguments();

		/*
		 * IO Shack
		 */
		i_progArgShackDict->getAutoRegistration(&shackIO);
		i_progArgShackDict->getAutoRegistration(&shackTimeTree);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i_progArgShackDict);


		fileOutput.setup(shackIO, shackTimeTree);

		return true;
	}

	bool clear()	override
	{
		if (benchmarks != nullptr)
		{
			delete benchmarks;
			benchmarks = nullptr;
		}

		if (timeTree != nullptr)
		{
			delete timeTree;
			timeTree = nullptr;
		}

		if (config != nullptr)
		{
			delete config;
			config = nullptr;
		}

		return true;
	}


	~DE_Dahlquist()
	{
		clear();
	}

	bool setTimeStepSize(
			double i_dt
	) override
	{
		SWEET_ASSERT(timeTree != nullptr);
		SWEET_ASSERT(timeTree->timeIntegrator != nullptr);
		timeTree->timeIntegrator->setTimeStepSize(i_dt);
		return true;
	}

	bool getInitialState(
		sweet::Data::GenericContainer::Base &o_U
	) override
	{
		SWEET_ASSERT(benchmarks != nullptr);
		benchmarks->benchmark->getInitialState(o_U);
		return true;
	}


	bool runTimestep(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
		) override
	{
		return timeTree->runIntegration(i_U, o_U, i_simulationTime);
	}


	bool print(
			const sweet::Data::GenericContainer::Base &i_data,
			std::string &i_type,
			std::iostream &io_iostream
	) override
	{
		io_iostream << cast(i_data).U;
		return true;
	}


	bool runFileOutput(
			const sweet::Data::GenericContainer::Base &i_data
	) override
	{
		///////std::ostringstream oss;

		///////std::vector<char> buffer;
		///////buffer.resize(shackIO->outputFileName.length()+1000);

		///////std::sprintf(
		///////		buffer.data(),
		///////		shackIO->outputFileName.c_str(),
		///////		"U",
		///////		shackTimeTree->currentSimulationTime*shackIO->outputFormatTimeScale
		///////	);

		///////std::string filename = buffer.data();

		const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &data = cast(i_data);

		////std::cout << "Writing to file '" << buffer.data() << "'" << std::endl;
		///////data.fileSave(filename);
		////fileOutput.fileSave(data.U, filename);
		fileOutput.fileSave(data.U, "prog_u");

		return true;
	}


	double computeError(
			const sweet::Data::GenericContainer::Base &i_data,
			double i_timeStamp,
			const std::string &i_errorNorm
	) override
	{
		//if (i_errorNorm == "lmax")
		{
			sweet::Data::GenericContainer::Base *b = getNewDataContainerInstance();
			benchmarks->benchmark->getReferenceState(*b, i_timeStamp);

			b->op_subVector(i_data);

			return std::abs(cast(*b).U);
		}

		return error.set("Norm '"+i_errorNorm+"' not supported");
		return true;
	}

	sweet::Data::GenericContainer::Base* getNewDataContainerInstance(
			int i_id = DataContainer::Simulation::DATA_SIMULATION
	) override
	{
		SWEET_ASSERT(config != nullptr);
		return config->getNewDataContainerInstance(i_id);
	}
};


}

}

#endif
