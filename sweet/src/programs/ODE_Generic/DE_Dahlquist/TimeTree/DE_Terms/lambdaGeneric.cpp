#include "lambdaGeneric.hpp"

#include <vector>
#include <string>

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DE_Terms {

lambdaGeneric::lambdaGeneric(const std::string &i_termId)	:
	shackODEGeneric_DE_Dahlquist(nullptr),
	_termIdStr(i_termId)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EXPONENTIAL);
	setEvalAvailable(EVAL_REXI_TERM);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_SEMI_LAGRANGIAN);

	_semiLagrangian_order = -1;
}

lambdaGeneric::~lambdaGeneric()
{
}


lambdaGeneric::lambdaGeneric(
		const lambdaGeneric &i_val
)	:
	TimeTree_Node_LeafHelper(i_val)
{
	shackODEGeneric_DE_Dahlquist = i_val.shackODEGeneric_DE_Dahlquist;
	_termIdStr = i_val._termIdStr;
	_expFunction = i_val._expFunction;
	_semiLagrangian_order = i_val._semiLagrangian_order;
}


bool lambdaGeneric::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	shackODEGeneric_DE_Dahlquist = io_shackDict->getAutoRegistration<Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> lambdaGeneric::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back(_termIdStr);
	return retval;

}


bool lambdaGeneric::setupConfigAndForwardTimeStepperEval(
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

	if (_termIdStr == "l1")
			_lambda = shackODEGeneric_DE_Dahlquist->lambda1;
	else if (_termIdStr == "l2")
			_lambda = shackODEGeneric_DE_Dahlquist->lambda2;
	else if (_termIdStr == "l3")
			_lambda = shackODEGeneric_DE_Dahlquist->lambda3;
	else
		return error.set("Unknown term id '"+_termIdStr+"'");

	if (!_expFunction.isSetup())
		_expFunction.setup("phi0");

	return true;
}


void lambdaGeneric::clear()
{
	TimeTree_Node_LeafHelper::clear();
}





bool lambdaGeneric::setTimeStepSize(double i_dt)
{
	_timestepSize = i_dt;
	return true;
}




bool lambdaGeneric::setupByKeyValue(
		const std::string &i_key,
		const std::string &i_value
)
{
	if (i_key == "ExpIntegrationFunction")
	{
		if (_expFunction.functionName != "")
			return error.set("Function name for expFunction is already set ('"+_expFunction.functionName+"')");

		if (i_value == "")
			return error.set("Empty string for expFunction given");

		_expFunction.setup(i_value);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(_expFunction);

		return true;
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


bool lambdaGeneric::setupByKeyValue(
		const std::string &i_key,
		const std::complex<double> &i_value
)
{
	if (i_key == "rexiTermAlpha")
	{
		_rexiTermAlpha = i_value;
		return true;
	}

	if (i_key == "rexiTermBeta")
	{
		_rexiTermBeta = i_value;
		return true;
	}

	if (i_key == "rexiTermGamma")
	{
		_rexiTermGamma = i_value;
		_rexiTermGammaActive = true;
		return true;
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}



bool lambdaGeneric::setupTreeNodeByFunction(
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

		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
			if (a->key == "order" || a->key == "sl_order" || a->key == "o")
			{
				a->getValue(_semiLagrangian_order);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
				break;
			}

			return error.set("Key not supported"+a->getNewLineDebugMessage());
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




bool lambdaGeneric::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

#if !SWEET_XBRAID
	std::cout << i_timeStamp << ": " << i_U.U << std::endl;
#endif
	o_U.U = _lambda * i_U.U;

	return true;
}

bool lambdaGeneric::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

	_expFunction.eval(_dt*_lambda, o_U.U);
	o_U.U *= i_U.U;

	return true;
}


bool lambdaGeneric::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const ODE_Generic::DE_Dahlquist::DataContainer::Simulation &i_U = cast(i_U_);
	ODE_Generic::DE_Dahlquist::DataContainer::Simulation &o_U = cast(o_U_);

	o_U.U = i_U.U / (1.0 - _dt * _lambda);

	return true;
}

bool lambdaGeneric::_eval_rexiTerm(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	/*
	 * We can reuse the backward Euler time stepper which has the form
	 *
	 * U1 = (I - dt*L)^{-1} U0
	 *
	 * For REXI, we need to solve terms of the form
	 *
	 * U1 = \beta(dt*L - \alpha)^{-1} U0
	 *
	 * and rewrite it to
	 *
	 * U1 = -\beta / \alpha (I - (dt / \alpha)*L)^{-1} U0
	 */
	std::complex<double> dtComplex = _dt/_rexiTermAlpha;

	std::complex<double> foo = -_rexiTermBeta/_rexiTermAlpha;

	o_U.U = foo / (1.0 - dtComplex * _lambda) * i_U.U;

	if (_rexiTermGammaActive)
		o_U.U += _rexiTermGamma*i_U.U;

	return true;
}


bool lambdaGeneric::evalNA_getNumStates(
		int *o_numStates
)
{
	*o_numStates = 2;
	return true;
}


bool lambdaGeneric::evalNA_departurePoints(
		const sweet::Data::GenericContainer::Base* i_U[],	//!< Vector of states
		double i_timestepSize,
		sweet::Data::GenericContainer::Base &o_departurePositions		//!< Computed departure positions
)
{
	DataContainer::SemiLagPositions &departurePoints = static_cast<DataContainer::SemiLagPositions&>(o_departurePositions);

	departurePoints.data[0] = 1505;
	return true;
}

bool lambdaGeneric::evalNA_interpolate(
		const sweet::Data::GenericContainer::Base &i_U,		//!< Input simulation data
		const sweet::Data::GenericContainer::Base &i_samplingPositions,	//!< Sampling positions (computed by _evalNA_departurePoints)
		sweet::Data::GenericContainer::Base &o_U_sampled			//!< Output samples
)
{
	const DataContainer::SemiLagPositions &departurePoints = static_cast<const DataContainer::SemiLagPositions&>(i_samplingPositions);

	/*
	 * Ugly hack: We take into account this effect by simply using an exponential integration
	 */
	if (!_expFunction.isSetup())
		_expFunction.setup("phi0");

	_eval_exponential(i_U, o_U_sampled, 0);
	return true;
	if (departurePoints.data[0] != 1505)
		SWEETErrorFatal("Magic code not found");

	o_U_sampled.op_setVector(i_U);

	return true;
}


}}}

