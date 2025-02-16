/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_SCALAR_TIME_ODESCALARTS_LN_ERK_HPP
#define PROGRAMS_ODE_SCALAR_TIME_ODESCALARTS_LN_ERK_HPP

#include <programs/ODE_Scalar/time/ODEScalarTS_BaseInterface.hpp>
#include <programs/ODE_Scalar/time/ShackODEScalarTimeDiscretization.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../benchmarks/ShackODEScalarBenchmarks.hpp"

class ODEScalarTS_ln_erk	:
		public ODEScalarTS_BaseInterface
{
public:
	double param_a;
	double param_b;

public:
	ODEScalarTS_ln_erk();

	virtual ~ODEScalarTS_ln_erk();

	bool setup();


public:
	void runTimestep(
			double &io_u,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);
};

#endif
