#include "null.hpp"

#include <vector>
#include <string>

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DE_Terms {

null::null()	:
	shackODEGeneric_DE_Dahlquist(nullptr)
{
	setEvalAvailable(EVAL_INTEGRATION);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_TENDENCIES);
}

null::~null()
{
}


null::null(
		const null &i_val
)	:
	TimeTree_Node_LeafHelper(i_val)
{
	shackODEGeneric_DE_Dahlquist = i_val.shackODEGeneric_DE_Dahlquist;
}


bool null::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	shackODEGeneric_DE_Dahlquist = io_shackDict->getAutoRegistration<Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> null::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("null");
	return retval;

}


bool null::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun *o_timeStepper
)
{
	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
	return true;
}


void null::clear()
{
	TimeTree_Node_LeafHelper::clear();
}



bool null::setTimeStepSize(double i_dt)
{
	_timestepSize = i_dt;
	return true;
}



bool null::setupTreeNodeByFunction(
		std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
		sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
)
{
	for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
	{
		sweet::TimeTree::TimeTreeIR::Argument *a = iter->get();

		switch(a->argType)
		{
		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
			error.set("Time steppers inside this time stepper are not allowed, yet"+a->getNewLineDebugMessage());
			return false;
			break;

		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
			error.set("Just values as arguments are not supported, yet"+a->getNewLineDebugMessage());
			break;

		default:
			SWEETErrorFatal("Internal error");
			return error.set("Internal error");
		}
	}

	// provide debug message in case that something goes wrong with the arguments
	setDebugMessage(i_function->getDebugMessage());
	return true;
}


bool null::_eval_integration(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

	o_U.U = i_U.U;
	return true;
}


bool null::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	//const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

	o_U.U = 0;
	return true;
}


bool null::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

	o_U.U = i_U.U;
	return true;
}

}}}

