/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#ifndef INCLUDE_SWEET_TOOLS_STOPWATCHBOX_HPP
#define INCLUDE_SWEET_TOOLS_STOPWATCHBOX_HPP

#ifndef SWEET_BENCHMARK_TIMINGS
#define SWEET_BENCHMARK_TIMINGS 1
#endif

#include <iostream>
#include <sweet/Tools/Stopwatch.hpp>

namespace sweet {
namespace Tools {

class StopwatchBox
{
public:
	sweet::Tools::Stopwatch main;
	sweet::Tools::Stopwatch main_setup;
	sweet::Tools::Stopwatch main_timestepping;

#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::Stopwatch main_timestepping_nonlinearities;

	sweet::Tools::Stopwatch rexi;
	sweet::Tools::Stopwatch rexi_setup;
	sweet::Tools::Stopwatch rexi_shutdown;
	sweet::Tools::Stopwatch rexi_timestepping;
	sweet::Tools::Stopwatch rexi_timestepping_solver;
	sweet::Tools::Stopwatch rexi_timestepping_broadcast;
	sweet::Tools::Stopwatch rexi_timestepping_reduce;
	sweet::Tools::Stopwatch rexi_timestepping_miscprocessing;


	sweet::Tools::Stopwatch main_timestepping_semi_lagrangian;
#endif



	static StopwatchBox& getInstance()
	{
		static StopwatchBox instance;
		return instance;
	}

	void reset()
	{
		main.reset();
		main_setup.reset();
		main_timestepping.reset();
#if SWEET_BENCHMARK_TIMINGS
		main_timestepping_nonlinearities.reset();


		rexi.reset();
		rexi_setup.reset();
		rexi_shutdown.reset();
		rexi_timestepping.reset();
		rexi_timestepping_solver.reset();
		rexi_timestepping_broadcast.reset();
		rexi_timestepping_reduce.reset();
		rexi_timestepping_miscprocessing.reset();

		main_timestepping_semi_lagrangian.reset();
#endif
	}


	void output()
	{
		if (main() != 0 || main_setup() != 0 || main_timestepping() != 0)
		{
			std::cout << "[MULE] simulation_benchmark_timings.main: " << main() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.main_setup: " << main_setup() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.main_timestepping: " << main_timestepping() << std::endl;
#if SWEET_BENCHMARK_TIMINGS
			std::cout << "[MULE] simulation_benchmark_timings.main_timestepping_nonlinearities: " << main_timestepping_nonlinearities() << std::endl;
#endif
		}

#if SWEET_BENCHMARK_TIMINGS
		if (
				rexi() != 0 ||
				rexi_setup() != 0 ||
				rexi_shutdown() != 0 ||
				rexi_timestepping() != 0 ||
				rexi_timestepping_solver() != 0 ||
				rexi_timestepping_broadcast() != 0 ||
				rexi_timestepping_reduce() != 0 ||
				rexi_timestepping_miscprocessing() != 0
		)
		{
			std::cout << "[MULE] simulation_benchmark_timings.rexi: " << rexi() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_setup: " << rexi_setup() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_shutdown: " << rexi_shutdown() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_timestepping: " << rexi_timestepping() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_timestepping_solver: " << rexi_timestepping_solver() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_timestepping_broadcast: " << rexi_timestepping_broadcast() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_timestepping_reduce: " << rexi_timestepping_reduce() << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.rexi_timestepping_miscprocessing: " << rexi_timestepping_miscprocessing() << std::endl;

			std::cout << "[MULE] simulation_benchmark_timings.semi_lagrangian: " << main_timestepping_semi_lagrangian() << std::endl;
		}
#endif
	}


	StopwatchBox()
	{
		reset();
	}
};

}}

#endif
