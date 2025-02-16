#include <sstream>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTreeIR_2_TimeTreeNodes.hpp>

namespace sweet {
namespace TimeTree {

class TimeTreeIR_2_TimeTreeNodes;


const std::string TimeTree_Node_Base::evalTypeToString(TIME_STEPPER_TYPES i_evalType)
{
	switch(i_evalType)
	{
		case EVAL_NONE:		return "[NONE]";

		case EVAL_INTEGRATION:	return "integration";
		case EVAL_TENDENCIES:	return "tendencies";
		case EVAL_EULER_BACKWARD:	return "euler_backward";
		case EVAL_REXI_TERM:	return "rexi_term";
		//case EVAL_EULER_FORWARD:	return "euler_forward";
		case EVAL_EXPONENTIAL:	return "exponential";
		case EVAL_SEMI_LAGRANGIAN:	return "semi_lagrangian";
	}
	return ">ERROR<";
}

TimeTree_Node_Base::TimeTree_Node_Base()	:
	_evalTypeRequested(EVAL_NONE)
{
}


TimeTree_Node_Base::~TimeTree_Node_Base()
{
}

/*!
 * Copy constructor
 */

TimeTree_Node_Base::TimeTree_Node_Base(
		const TimeTree_Node_Base &i_src
)
{
	_registeredEvalTypes = i_src._registeredEvalTypes;
	_evalTypeRequested = i_src._evalTypeRequested;
}



bool TimeTree_Node_Base::setupTreeNodeByFunction(
		std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
		sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
)
{
	return error.set("setupTreeNodeByFunction not supported for this node");
}




bool TimeTree_Node_Base::setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		EvalFun *o_timeStepper = nullptr
	)
{
	return setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_INTEGRATION, o_timeStepper);
}



bool TimeTree_Node_Base::setupByKeyValue(
		const std::string &i_key,
		const std::string &i_value
){
	return error.set("setupByKeyValue(std::string,std::string) not available in tree node '"+getNodeNames()[0]+"'");
};



bool TimeTree_Node_Base::setupByKeyValue(
		const std::string &i_key,
		const double &i_value
){
	return error.set("setupByKeyValue(std::string,double) not available in tree node '"+getNodeNames()[0]+"'");
};


bool TimeTree_Node_Base::setupByKeyValue(
		const std::string &i_key,
		const std::complex<double> &i_value
){
	return error.set("setupByKeyValue(std::string,std::complex<double>) not available in tree node '"+getNodeNames()[0]+"'");
};


bool TimeTree_Node_Base::registerTimeStepperEval(
		TIME_STEPPER_TYPES i_evalType	//!< Type of time integrator
)
{
	if (_evalTypeRequested != EVAL_NONE)
		return error.set("Eval already requested for this node! Only one eval request is supported!");

	_evalTypeRequested = i_evalType;

	return true;
}


bool TimeTree_Node_Base::registerTimeStepperEval(
		TIME_STEPPER_TYPES i_evalType,	//!< Type of time integrator
		EvalFun *o_timeStepper	//!< Callback of time integrator
)
{
	if (_evalTypeRequested != EVAL_NONE)
		return error.set("Eval already requested for this node! Only one eval request is supported!");

	_evalTypeRequested = i_evalType;

	if (!isEvalAvailable(i_evalType))
	{
		std::ostringstream oss;
		oss << "Time stepper evaluation '" << evalTypeToString(i_evalType) << "' not available or not registered";
		return error.set(oss.str());
	}

	if (o_timeStepper == nullptr)
		return true;

	if (i_evalType == EVAL_INTEGRATION)
	{
		*o_timeStepper = &TimeTree_Node_Base::_eval_integration;
		return true;
	}

	if (i_evalType == EVAL_TENDENCIES)
	{
		*o_timeStepper = &TimeTree_Node_Base::_eval_tendencies;
		return true;
	}

	if (i_evalType == EVAL_EULER_BACKWARD)
	{
		*o_timeStepper = &TimeTree_Node_Base::_eval_eulerBackward;
		return true;
	}

	if (i_evalType == EVAL_REXI_TERM)
	{
		*o_timeStepper = &TimeTree_Node_Base::_eval_rexiTerm;
		return true;
	}

#if 0
	if (i_evalType == EVAL_EULER_FORWARD)
	{
		*o_timeStepper = &TimeTree_Node_Base::_eval_eulerForward;
		return true;
	}
#endif

	if (i_evalType == EVAL_EXPONENTIAL)
	{
		*o_timeStepper = &TimeTree_Node_Base::_eval_exponential;
		return true;
	}

	if (i_evalType == EVAL_SEMI_LAGRANGIAN)
	{
		if (o_timeStepper != nullptr)
			return error.set("evalType EVAL_SEMI_LAGRANGIAN doesn't support the return of an eval function!");
		return true;
	}

	std::ostringstream oss;
	oss << "Invalid Time stepper evaluation '" << i_evalType << "'";
	return error.set(oss.str());
}


void TimeTree_Node_Base::clear()
{
	_evalTypeRequested = EVAL_NONE;

#if SWEET_XBRAID
	if (U_prev_solution)
	{
		U_prev_solution->clear();
		U_prev_solution = nullptr;
	}
#endif

}



bool TimeTree_Node_Base::isEvalAvailable(TIME_STEPPER_TYPES i_evalType)
{
	for (auto e: _registeredEvalTypes)
	{
		if (e == i_evalType)
			return true;
	}

	return false;
}


std::string TimeTree_Node_Base::getSupportedEvalsAsString()
{
	std::ostringstream oss;

	for (auto &e: _registeredEvalTypes)
		oss << evalTypeToString(e) << ", ";

	std::string retval = oss.str();
	if (retval.size() > 2)
		retval = retval.substr(0, retval.size()-2);

	return retval;
}

void TimeTree_Node_Base::setEvalAvailable(TIME_STEPPER_TYPES i_evalType)
{
	_registeredEvalTypes.push_back(i_evalType);
}


bool TimeTree_Node_Base::_eval_integration(
		const sweet::Data::GenericContainer::Base &i_U,
		sweet::Data::GenericContainer::Base &o_U,
		double i_simulationTime
) {
	return error.set("_eval_integration() not available");
};


bool TimeTree_Node_Base::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U,
		sweet::Data::GenericContainer::Base &o_U,
		double i_timeStamp
){
	return error.set("_eval_tendencies() not available");
};


bool TimeTree_Node_Base::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U,
		sweet::Data::GenericContainer::Base &o_U,
		double i_timeStamp
){
	return error.set("_eval_eulerBackward() not available");
};


bool TimeTree_Node_Base::_eval_rexiTerm(
		const sweet::Data::GenericContainer::Base &i_U,
		sweet::Data::GenericContainer::Base &o_U,
		double i_timeStamp
){
	return error.set("_eval_rexiTerm() not available");
};


bool TimeTree_Node_Base::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U,
		sweet::Data::GenericContainer::Base &o_U,
		double i_timeStamp
){
	return error.set("_eval_exponential() not available");
};


bool TimeTree_Node_Base::evalNA_getNumStates(
		int *o_numStates
){
	return error.set("_evalNA_getNumStates() not available");
};


bool TimeTree_Node_Base::evalNA_departurePoints(
		const sweet::Data::GenericContainer::Base* i_states[],
		double i_timestepSize,
		sweet::Data::GenericContainer::Base &o_departurePositions
){
	return error.set("_evalNA_departurePoints() not available");
};


bool TimeTree_Node_Base::evalNA_interpolate(
		const sweet::Data::GenericContainer::Base &i_U_input,
		const sweet::Data::GenericContainer::Base &i_samplingPositions,
		sweet::Data::GenericContainer::Base &i_U_samples
){
	return error.set("_evalNA_interpolate() not available");
};


std::string TimeTree_Node_Base::setDebugMessage(
		const std::string &i_debugMessage	//!< Message for debugging
)
{
	_debugMessage = i_debugMessage;
	return _debugMessage;
}


std::string TimeTree_Node_Base::getNewLineDebugMessage()
{
	if (_debugMessage != "")
		return "\n"+_debugMessage;

	return _debugMessage;
}


}}
