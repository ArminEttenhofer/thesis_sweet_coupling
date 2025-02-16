#include "direct.hpp"

#include <vector>
#include <string>

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DE_Terms {

direct::direct()	:
	shackODEGeneric_DE_Dahlquist(nullptr)
{
	setEvalAvailable(EVAL_INTEGRATION);
	setEvalAvailable(EVAL_EXPONENTIAL);
}

direct::~direct()
{
}


direct::direct(
		const direct &i_val
)	:
	TimeTree_Node_LeafHelper(i_val)
{
	shackODEGeneric_DE_Dahlquist = i_val.shackODEGeneric_DE_Dahlquist;
}


bool direct::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	shackODEGeneric_DE_Dahlquist = io_shackDict->getAutoRegistration<Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> direct::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("direct");
	return retval;

}


bool direct::setupConfigAndForwardTimeStepperEval(
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

	_lambda = shackODEGeneric_DE_Dahlquist->lambda1
			+ shackODEGeneric_DE_Dahlquist->lambda2
			+ shackODEGeneric_DE_Dahlquist->lambda3;

	_mu = shackODEGeneric_DE_Dahlquist->mu;
	_phi = shackODEGeneric_DE_Dahlquist->phi;

	return true;
}


void direct::clear()
{
	TimeTree_Node_LeafHelper::clear();
}



bool direct::setTimeStepSize(double i_dt)
{
	_timestepSize = i_dt;
	return true;
}



bool direct::setupTreeNodeByFunction(
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


bool direct::_eval_integration(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

	double t0 = i_timeStamp;
	double t1 = i_timeStamp + _dt;

	const T &U0 = i_U.U;
	T &U1 = o_U.U;


	/*
	 * https://www.wolframalpha.com/input?i=solve+d%2Fdt+u%28t%29+%3D+%5Clambda+u%28t%29+%2B+%5Cmu+*+%5Cexp%28%5Cphi+*+t%29
	 */

	// Solving for c1 using (t0, U0) yields
	T c = (U0 - (_mu * std::exp(_phi*t0)) / (_phi - _lambda) ) * std::exp(-_lambda*t0);

	U1 = c*std::exp(_lambda*t1) + (_mu*std::exp(_phi*t1))/(_phi - _lambda);

	return true;
}


bool direct::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	return _eval_integration(i_U_, o_U_, i_timeStamp);
}

}}}

