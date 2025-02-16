#ifndef PROGRAMS_PDE_SWESPHERE2D_DATACONTAINER_CONFIG_HPP
#define PROGRAMS_PDE_SWESPHERE2D_DATACONTAINER_CONFIG_HPP

#include <sweet/Data/GenericContainer/ConfigBase.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2DComplex/Operators.hpp>

// Two different data containers
#include "SemiLagPositions.hpp"
#include "Simulation.hpp"
#include "Topography.hpp"


namespace PDE_SWESphere2D {
namespace DataContainer {

class Config:
		public sweet::Data::GenericContainer::ConfigBase
{
public:
	/*
	 * Just a pointer to an existing data container
	 */
	const Simulation *myDataContainer;

	sweet::Data::Sphere2D::Operators *ops;
	sweet::Data::Sphere2DComplex::Operators *opsComplex;


	Config()
	{
		myDataContainer = nullptr;
		//myDataContainerComplex = nullptr;
		ops = nullptr;
		opsComplex = nullptr;
	}


	/*!
	 * Return a new instance of a data container.
	 *
	 * This is what will be used by time steppers and/or DE term implementations
	 */
	sweet::Data::GenericContainer::Base* getNewDataContainerInstance(
			int i_id = -1
	) const override
	{
		if (i_id == -1 || i_id == sweet::Data::GenericContainer::Base::DATA_SIMULATION)
		{
			Simulation *retval = new Simulation;
			retval->setup(ops->sphere2DDataConfig);
			return retval;
		}


		if (i_id == -1 || i_id == sweet::Data::GenericContainer::Base::DATA_SEMI_LAGRANGIAN_POSITIONS)
		{
			/*
			 * Real valued
			 */
			SemiLagPositions *retval = new SemiLagPositions;
			retval->setup(ops->sphere2DDataConfig);
			return retval;
		}

		if (i_id == -1 || i_id == sweet::Data::GenericContainer::Base::DATA_TOPOGRAPHY)
		{
			/*
			 * Real valued
			 */
			Topography *retval = new Topography;
			retval->setup(ops->sphere2DDataConfig);
			return retval;
		}


		SWEETErrorFatal("Invalid data type id");
		return 0;
	}
};

}}

#endif
