#include <sweet/Data/Cart2D/Operators.hpp>
#include "n.hpp"


#include <vector>


namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {

n::n()	:
	_shackPDESWECart2D(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

n::~n()
{
}

n::n(
		const n &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWECart2D = i_value._shackPDESWECart2D;
	_ops = i_value._ops;
}


bool n::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWECart2D = io_shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> n::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("n");
	return retval;

}


bool n::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun* o_timeStepper
)
{
	const DataContainer::Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}

void n::clear()
{
	TimeTree_Node_LeafHelper::clear();
}

/*
 * Return the time tendencies of the PDE term
 */
bool n::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWECart2D != nullptr);

	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_U.u = -i_U.u*_ops->diff_c_x(i_U.u) - i_U.v*_ops->diff_c_y(i_U.u);
	o_U.v = -i_U.u*_ops->diff_c_x(i_U.v) - i_U.v*_ops->diff_c_y(i_U.v);

	if (_shackPDESWECart2D->use_only_linear_divergence) //only nonlinear advection left to solve
		o_U.h_pert = - (i_U.u*_ops->diff_c_x(i_U.h_pert) + i_U.v*_ops->diff_c_y(i_U.h_pert));
	else //full nonlinear equation on h
		o_U.h_pert = -_ops->diff_c_x(i_U.u*i_U.h_pert) - _ops->diff_c_y(i_U.v*i_U.h_pert);

	return true;
}

}}}
