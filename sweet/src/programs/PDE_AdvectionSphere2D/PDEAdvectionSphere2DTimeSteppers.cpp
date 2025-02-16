/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DTimeSteppers.hpp>
#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_na_erk.hpp>
#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_na_sl.hpp>
#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_na_trajectories.hpp>




PDEAdvectionSphere2DTimeSteppers::PDEAdvectionSphere2DTimeSteppers()
{
}

void PDEAdvectionSphere2DTimeSteppers::setup_1_registerAllTimesteppers()
{
	/*
	 * Register time integrators
	 */
	_registered_integrators.push_back(static_cast<PDEAdvectionSphere2DTS_BaseInterface*>(new PDEAdvectionSphere2DTS_na_erk));
	_registered_integrators.push_back(static_cast<PDEAdvectionSphere2DTS_BaseInterface*>(new PDEAdvectionSphere2DTS_na_sl));
	_registered_integrators.push_back(static_cast<PDEAdvectionSphere2DTS_BaseInterface*>(new PDEAdvectionSphere2DTS_na_trajectories));
}



bool PDEAdvectionSphere2DTimeSteppers::setup_2_shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->shackRegistration(io_shackDict);
	}
	return true;
}



void PDEAdvectionSphere2DTimeSteppers::printImplementedTimesteppingMethods(
	std::ostream &o_ostream,
	const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	std::string prefix = i_prefix+"  ";
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->printImplementedTimesteppingMethods(o_ostream, prefix);
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << std::endl;
}


bool PDEAdvectionSphere2DTimeSteppers::setup_3_timestepper(
		const std::string &i_timestepping_method,
		sweet::Shacks::Dictionary *io_shackDict,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	if (i_timestepping_method == "")
	{
		printImplementedTimesteppingMethods();
		return error.set("Please set time stepping method using --timestepping-method=...");
	}
	/*
	 * Find right one
	 */
	timestepper = nullptr;

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphere2DTS_BaseInterface *ts = _registered_integrators[i];

		if (ts->testImplementsTimesteppingMethod(i_timestepping_method))
		{
			if (timestepper != nullptr)
			{
				//std::cout << "Processing " << i+1 << "th element" << std::endl;
				return error.set("Duplicate implementation for method "+i_timestepping_method);
			}

			//std::cout << "Found matching time stepping method at " << i+1 << "th element" << std::endl;
			ts->setup(io_ops);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*ts);
			timestepper = ts;
		}
	}

	if (timestepper == nullptr)
		return error.set("No valid --timestepping-method '"+i_timestepping_method+"' provided");

	// Found integrator, freeing others
	_timesteppersFreeAll(timestepper);

	return true;
}


void PDEAdvectionSphere2DTimeSteppers::_timesteppersFreeAll(
		PDEAdvectionSphere2DTS_BaseInterface *i_skip_this_timestepper
)
{

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphere2DTS_BaseInterface *ts = _registered_integrators[i];

		if (ts == i_skip_this_timestepper)
			continue;

		delete ts;
	}

	_registered_integrators.clear();
}


void PDEAdvectionSphere2DTimeSteppers::clear()
{
	delete timestepper;
	timestepper = nullptr;

	_timesteppersFreeAll();
}


PDEAdvectionSphere2DTimeSteppers::~PDEAdvectionSphere2DTimeSteppers()
{
	clear();
}
