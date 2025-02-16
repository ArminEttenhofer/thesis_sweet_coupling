#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_REGISTRY_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_REGISTRY_HPP


#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/InteriorNodes/AddIntegration.hpp>
#include <sweet/TimeTree/InteriorNodes/AddTendencies.hpp>
#include <sweet/TimeTree/InteriorNodes/ETDRK.hpp>
#include <sweet/TimeTree/InteriorNodes/SLETDRK.hpp>
#include <sweet/TimeTree/InteriorNodes/ExplicitRungeKutta.hpp>
#include <sweet/TimeTree/InteriorNodes/Exponential.hpp>
#include <sweet/TimeTree/InteriorNodes/ImplicitRungeKutta.hpp>
#include <sweet/TimeTree/InteriorNodes/Lawson.hpp>
#include <sweet/TimeTree/InteriorNodes/NegIntegration.hpp>
#include <sweet/TimeTree/InteriorNodes/NegTendencies.hpp>
#include <sweet/TimeTree/InteriorNodes/REXI.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_FP.hpp>
// #include <sweet/TimeTree/InteriorNodes/SDC_Martin.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_Classic.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_ETD.hpp>
#include <sweet/TimeTree/InteriorNodes/SETTLS.hpp>
#include <sweet/TimeTree/InteriorNodes/StrangSplitting.hpp>
#include <sweet/TimeTree/InteriorNodes/SubCycling.hpp>
#include <sweet/TimeTree/TimeTree_Node_Registry.hpp>



namespace sweet {
namespace TimeTree {

/*!
 * \brief Interior nodes of time tree (time steppers and composers)
 */
namespace InteriorNodes {}


/*!
 * This is a class to register all time steppers
 */
class InteriorNodes_Registry
{
public:
	Error::Base error;

public:
	InteriorNodes_Registry()
	{
	}

	bool registerAll(
			sweet::TimeTree::TimeTree_Node_Registry &o_timeStepper_registry
	)
	{
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::AddTendencies>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::NegTendencies>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::AddIntegration>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::NegIntegration>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::ExplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::ImplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::StrangSplitting>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SubCycling>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::Exponential>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::REXI>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::ETDRK>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SLETDRK>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SETTLS>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SDC_FP>();
		// o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SDC_Martin>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SDC_Classic>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::SDC_ETD>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::TimeTree::InteriorNodes::LawsonRungeKutta>();

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(o_timeStepper_registry);
		return true;
	}
};

}}

#endif
