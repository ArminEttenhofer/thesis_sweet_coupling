/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */
#ifndef PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_DEFAULT_HPP
#define PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_DEFAULT_HPP


#include <programs/ODE_Scalar/benchmarks/ODEScalarBench_BaseInterface.hpp>
#include <programs/ODE_Scalar/benchmarks/ShackODEScalarBenchmarks.hpp>
#include <stdlib.h>
#include <cmath>


/**
 * Setup Default benchmark
 */
class ODEScalarBenchDefault	:
		public ODEScalarBench_BaseInterface
{

public:
	bool setupBenchmark(
			double &o_u
	)
	{

		o_u = shackBenchmarks->u0;

		return true;
	}
};


#endif
