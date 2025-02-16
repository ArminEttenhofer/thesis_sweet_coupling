#include <sweet/Data/Sphere2D/Operators.hpp>
#include "null.hpp"

#include <complex>
#include <cmath>
#include <vector>


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


null::null()	:
	_shackPDESWESphere2D(nullptr),
	_shackSphere2DDataOps(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_EXPONENTIAL);
}


null::~null()
{
}


null::null(
		const null &i_src
)	:
	TimeTree_Node_LeafHelper(i_src)
{
	_shackSphere2DDataOps = i_src._shackSphere2DDataOps;
	_shackPDESWESphere2D = i_src._shackPDESWESphere2D;
	_ops = i_src._ops;
}


bool null::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> null::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("null");
	return retval;

}



bool null::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Null term - doing nothing :-)':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}


bool null::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun *o_timeStepper
)
{
	const DataContainer::Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	SWEET_ASSERT(_ops != nullptr);

	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);

	return true;
}


void null::clear()
{
	TimeTree_Node_LeafHelper::clear();
}


bool null::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	o_U_.op_setZero();

	return true;
}



/*
 * Backward Euler evaluation of the term
 *
 * U1' = L U1 = (U1-U0)/dt
 * <=> L U1 dt = U1 - U0
 * <=> U1 - L U1 dt = U0
 * <=> (I - dt*L) U1 = U0
 *
 * U1 = (I - dt*L)^{-1} U0
 */
bool null::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	o_U_.op_setVector(i_U_);

	return true;
}



/*!
 * Evaluate exponential of linear term
 */
bool null::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	o_U_.op_setVector(i_U_);

	return true;
}

}}}
