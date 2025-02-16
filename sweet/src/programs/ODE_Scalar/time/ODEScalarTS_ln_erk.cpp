/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#include <programs/ODE_Scalar/time/ODEScalarTS_ln_erk.hpp>


ODEScalarTS_ln_erk::ODEScalarTS_ln_erk()
{
}


ODEScalarTS_ln_erk::~ODEScalarTS_ln_erk()
{
}


/*
 * Setup
 */
bool ODEScalarTS_ln_erk::setup()
{
	param_a = shackODEScalarBenchmark->ode_parameters[0];
	param_b = shackODEScalarBenchmark->ode_parameters[1];
	return true;
}

void ODEScalarTS_ln_erk::runTimestep(
		double &io_u,		//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	SWEET_ASSERT(i_dt > 0);

	io_u += i_dt * (param_a * std::sin(io_u) + param_b * std::sin(i_simulation_timestamp));

}


