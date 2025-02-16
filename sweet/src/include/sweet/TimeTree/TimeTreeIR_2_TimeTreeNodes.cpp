#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTreeIR_2_TimeTreeNodes.hpp>

namespace sweet {
namespace TimeTree {

TimeTreeIR_2_TimeTreeNodes::TimeTreeIR_2_TimeTreeNodes()	:
	deTermsRegistry(nullptr),
	interiorNodesRegistry(nullptr)
{
}


void TimeTreeIR_2_TimeTreeNodes::clear()
{
}


bool TimeTreeIR_2_TimeTreeNodes::setup(
	TimeTree_Node_Registry &i_deTerms,		//!< Leaf nodes (DE terms)
	TimeTree_Node_Registry &i_timeSteppers	//!< Interior nodes (time integrators)
)
{
	deTermsRegistry = &i_deTerms;
	interiorNodesRegistry = &i_timeSteppers;

	return true;
}



bool TimeTreeIR_2_TimeTreeNodes::assembleTimeTree(
		TimeTreeIR &i_tree,
		std::shared_ptr<TimeTree_Node_Base> &o_timestepper
)
{
	if (deTermsRegistry == nullptr || interiorNodesRegistry == nullptr)
		return error.set("You need to call setup(...) before assembleTimeTree(...)");

	return assembleTimeTreeNodeByFunction(
			i_tree.mainFunction,
			o_timestepper
		);
}


bool TimeTreeIR_2_TimeTreeNodes::assembleTimeTreeNodeByFunction(
	std::shared_ptr<TimeTreeIR::Function> &i_function,
	std::shared_ptr<TimeTree_Node_Base> &o_timestepper
)
{
	/*
	 * Step 1) Search for implementation of time stepper of this particular function
	 *
	 * We first search the interior nodes registry.
	 * If we don't find anything, we also search the leaf node registry.
	 */
	bool retval = interiorNodesRegistry->getTimeTreeNodeNewInstance(
			i_function->function_name,
			o_timestepper,
			false
		);

	if (!retval)
	{
		deTermsRegistry->getTimeTreeNodeNewInstance(
				i_function->function_name,
				o_timestepper
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*deTermsRegistry);
	}

	/*
	 * Step 2) Call the setup routine of the time stepper because
	 * Only the time steppers knows about how to process its arguments.
	 *
	 * We also hand over this class since there could be other
	 *
	 */
	o_timestepper->setupTreeNodeByFunction(i_function, *this);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*(o_timestepper.get()));

	return true;
}


bool TimeTreeIR_2_TimeTreeNodes::assembleTimeTreeNodeByName(
	const std::string i_nodeName,
	std::shared_ptr<TimeTree_Node_Base> &o_timestepper
)
{
	if (i_nodeName == "")
		SWEETErrorFatal("Internal error: Congratulations! You broke SWEET. This shouldn't happen.");

	/*
	 * Step 1) Search for implementation of time stepper of this particular function
	 */

	// we first search for the de terms
	if (!deTermsRegistry->getTimeTreeNodeNewInstance(
			i_nodeName,
			o_timestepper,
			false
		)
	)
	{
		// Then we search for interior node handlers
		interiorNodesRegistry->getTimeTreeNodeNewInstance(
				i_nodeName,
				o_timestepper
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*interiorNodesRegistry);
	}

	return true;
}

}}
