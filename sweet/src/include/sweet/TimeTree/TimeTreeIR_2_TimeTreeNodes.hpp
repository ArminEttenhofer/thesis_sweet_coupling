

/*
 * Always include these classes due to forward delaration stuff
 */

#ifndef INCLUDE_SWEET_TIMETREE_TIMETREE_ASSEMBLATION_HPP
#define INCLUDE_SWEET_TIMETREE_TIMETREE_ASSEMBLATION_HPP


#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_Registry.hpp>
#include <sweet/TimeTree/TimeTreeIR.hpp>
#include <string>
#include <vector>
#include <ostream>
#include <memory>


namespace sweet {
namespace TimeTree {

/*
 * Forward declaration of Base.
 * It's included at the end of this file.
 */

class TimeTree_Node_Base;

/*!
 * \brief Assembles the final time stepper for a given time stepping tree
 *
 * This searches and creates the corresponding instances of all tree nodes.
 */
class TimeTreeIR_2_TimeTreeNodes
{
public:
	Error::Base error;

private:
	/*!
	 * Registry for leaf nodes (evaluating terms of differential equations)
	 */
	TimeTree_Node_Registry *deTermsRegistry;

	/*!
	 * Registry for interior nodes: Time integrators and composers
	 */
	TimeTree_Node_Registry *interiorNodesRegistry;

public:
	TimeTreeIR_2_TimeTreeNodes();

public:
	void clear();

	/*!
	 * Provide registries of interior and leaf nodes
	 */
public:
	bool setup(
		TimeTree_Node_Registry &i_deTerms,		//!< Leaf nodes (DE terms)
		TimeTree_Node_Registry &i_timeSteppers	//!< Interior nodes (time integrators)
	);


	/*!
	 * Setup the full time stepping tree
	 */
public:
	bool assembleTimeTree(
			TimeTreeIR &i_tree,
			std::shared_ptr<TimeTree_Node_Base> &o_timestepper
	);


	/*!
	 * Setup the time stepper for a given function and return it
	 */
public:
	bool assembleTimeTreeNodeByFunction(
		std::shared_ptr<TimeTreeIR::Function> &i_function,
		std::shared_ptr<TimeTree_Node_Base> &o_timestepper
	);


	/*!
	 * Setup the time stepper for a given node name and return it
	 */
public:
	bool assembleTimeTreeNodeByName(
		const std::string i_nodeName,
		std::shared_ptr<TimeTree_Node_Base> &o_timestepper
	);
};

}}


#endif
