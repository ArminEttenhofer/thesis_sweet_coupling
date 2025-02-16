/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIMETREE_TIMETREE_HPP
#define PROGRAMS_PDE_SWECART2D_TIMETREE_TIMETREE_HPP

#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/InteriorNodes_Registry.hpp>
#include <sweet/TimeTree/TimeString_2_TimeTreeIR.hpp>
#include <sweet/TimeTree/TimeTreeIR_2_TimeTreeNodes.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_Registry.hpp>
#include <sweet/TimeTree/TimeTreeIR.hpp>
#include "../DataContainer/Config.hpp"
#include "../DataContainer/Simulation.hpp"

namespace PDE_SWECart2D {
namespace TimeTree {


/*!
 * SWE Cart2D time tree
 */
class TimeTree
{
public:
	sweet::Error::Base error;

	//! Base node of time integration
	std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base> timeIntegrator;

	//! Evaluation function to time integrate
	sweet::TimeTree::TimeTree_Node_Base::EvalFun evalFun;

	//! A registry of all PDE terms
	sweet::TimeTree::TimeTree_Node_Registry pdeTerm_registry;

	//! A registry of all time steppers
	sweet::TimeTree::TimeTree_Node_Registry timeStepper_registry;

	//! Config file for DE solver
	DataContainer::Config deSolver_Config;

public:
	bool setup_1_registerAllTimesteppers();

public:
	TimeTree();

	bool setup_2_shackRegistration(
			sweet::Shacks::Dictionary *i_shackDict
	);

	bool setup_3_timestepper(
			const std::string &i_timestepping_method,
			sweet::Shacks::ProgramArgumentsDictionary *i_progArgShackDict,
			sweet::Data::Cart2D::Operators *io_ops,
			sweet::Data::Cart2DComplex::Operators *io_opsComplex,
			const DataContainer::Simulation &i_U
	);

	bool runIntegration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	);

	void clear();

	~TimeTree();
};

}}

#endif
