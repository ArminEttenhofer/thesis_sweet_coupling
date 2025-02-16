#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DATACONTAINER_CONFIG_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DATACONTAINER_CONFIG_HPP

#include <sweet/Data/GenericContainer/ConfigBase.hpp>

#include "Simulation.hpp"
#include "SemiLagPositions.hpp"

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DataContainer {

/*!
 * A special class which is forwarded to all
 *
 *  - time stepper instances and
 *  - DE term instances
 *
 * to *set up data buffers* and other things which are required.
 */
class Config:
		public sweet::Data::GenericContainer::ConfigBase
{
public:

	Config()
	{
	}


	/*!
	 * Return a new instance of a data container.
	 *
	 * This is what will be used by time steppers and/or DE term implementations
	 */
	sweet::Data::GenericContainer::Base* getNewDataContainerInstance(
			int i_id = DataContainer::Simulation::DATA_SIMULATION
	) const override
	{
		if (i_id == DataContainer::Simulation::DATA_SIMULATION)
		{
			/*
			 * Real valued
			 */
			ODE_Generic::DE_Dahlquist::DataContainer::Simulation *retval = new ODE_Generic::DE_Dahlquist::DataContainer::Simulation;
			retval->setup();
			return retval;
		}


		if (i_id == DataContainer::Simulation::DATA_SEMI_LAGRANGIAN_POSITIONS)
		{
			ODE_Generic::DE_Dahlquist::DataContainer::SemiLagPositions *retval = new ODE_Generic::DE_Dahlquist::DataContainer::SemiLagPositions;
			retval->setup();
			return retval;
		}

		SWEETErrorFatal("Invalid data type id");
		return 0;
	}
};

}}}

#endif
