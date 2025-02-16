/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_SCALAR_TIME_ODESCALARTS_BASEINTERFACE_HPP
#define PROGRAMS_ODE_SCALAR_TIME_ODESCALARTS_BASEINTERFACE_HPP

#include <programs/ODE_Scalar/time/ShackODEScalarTimeDiscretization.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include <sweet/TimeTree/Shack.hpp>

#include "../benchmarks/ShackODEScalarBenchmarks.hpp"
#include "../ShackODEScalar.hpp"

#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
#endif

class ODEScalarTS_BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackODEScalar *shackODEScalar;
	ShackODEScalarTimeDiscretization *shackODEScalarTimeDisc;
	ShackODEScalarBenchmarks *shackODEScalarBenchmark;

	ODEScalarTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackODEScalar(nullptr),
		shackODEScalarTimeDisc(nullptr),
		shackODEScalarBenchmark(nullptr)
	{

	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackODEScalar = io_shackDict->getAutoRegistration<ShackODEScalar>();
		shackODEScalarTimeDisc = io_shackDict->getAutoRegistration<ShackODEScalarTimeDiscretization>();
		shackODEScalarBenchmark = io_shackDict->getAutoRegistration<ShackODEScalarBenchmarks>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
	)
	{
		return true;
	}


public:
	virtual void runTimestep(
			double &io_u,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	) = 0;

	~ODEScalarTS_BaseInterface() {}

/*
 * Martin@Joao Please let's discuss this.
 */

public:

#if (SWEET_PARAREAL_SCALAR) || (SWEET_XBRAID_SCALAR)
	void runTimestep(
			sweet::DEPRECATED_pint::Parareal_GenericData* io_data,

			double i_dt,		//!< time step size
			double i_sim_timestamp
	)
	{
		double u = io_data->get_pointer_to_data_Scalar()->simfields[0];

		runTimestep(	u,
				i_dt,
				i_sim_timestamp
			);

		io_data->get_pointer_to_data_Scalar()->simfields[0] = u;

	}



	// for parareal SL
	virtual void set_previous_solution(
				double &i_u
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			sweet::DEPRECATED_pint::Parareal_GenericData* i_data
	)
	{
		double u = i_data->get_pointer_to_data_Scalar()->simfields[0];

		set_previous_solution(u);
	};
#endif


};



#endif
