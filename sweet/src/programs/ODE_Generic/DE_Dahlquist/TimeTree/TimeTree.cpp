#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/DE_Terms/lambdaGeneric.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/DE_Terms/direct.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/DE_Terms/forcing.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/DE_Terms/null.hpp>
#include "TimeTreeIR.hpp"


namespace ODE_Generic {
namespace DE_Dahlquist {


bool TimeTree::setup_1_registerAllTimesteppers()
{

	/*
	 * Registration of all possible PDE terms
	 */
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<ODE_Generic::DE_Dahlquist::DE_Terms::lambdaGeneric>("l1")
		);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<ODE_Generic::DE_Dahlquist::DE_Terms::lambdaGeneric>("l2")
		);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<ODE_Generic::DE_Dahlquist::DE_Terms::lambdaGeneric>("l3")
		);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<ODE_Generic::DE_Dahlquist::DE_Terms::forcing>()
		);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<ODE_Generic::DE_Dahlquist::DE_Terms::direct>()
		);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<ODE_Generic::DE_Dahlquist::DE_Terms::null>()
		);

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeTerm_registry);

	/*
	 * Register time steppers
	 */
	sweet::TimeTree::InteriorNodes_Registry timeSolvers_registry;
	timeSolvers_registry.registerAll(timeStepper_registry);

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSolvers_registry);

	return true;
}


TimeTree::TimeTree()
{
}

bool TimeTree::setup_2_shackRegistration(
		sweet::Shacks::Dictionary *i_shackDict
)
{
	pdeTerm_registry.shackRegistration(i_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeTerm_registry);

	timeStepper_registry.shackRegistration(i_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeStepper_registry);

	return true;
}

/////#if SWEET_XBRAID
/////bool TimeTree::setup_3_timestepper(
/////		const std::string &i_timestepping_method,	//!< String with time stepping method such as SS(ERK(lg,order=4),ERK(lc,order=2),order=2)
/////		sweet::Shacks::ProgramArgumentsDictionary* i_progArhShackDict,
/////		sweet::Data::Scalar::Operators* io_ops,
/////		sweet::Data::Scalar::Operators* io_ops_complex,
/////		const DataContainer::Simulation &i_U
/////)
/////{
/////	return setup_3_timestepper(i_timestepping_method, nullptr);
/////}
/////#endif


bool TimeTree::setup_3_timestepper(
		const std::string &i_timestepping_method,	//!< String with time stepping method such as SS(ERK(lg,order=4),ERK(lc,order=2),order=2)
		const sweet::Data::GenericContainer::ConfigBase *i_config
		/////const ODE_Generic::DE_Dahlquist::DataContainer::Config *i_config
)
{
	sweet::TimeTree::TimeTreeIR timeTree;

	/*
	 * Setup time stepping string parser and parse it
	 */
	sweet::TimeTree::TimeString_2_TimeTreeIR tsStringParser;
	tsStringParser.createTimeTree(
			i_timestepping_method,
			timeTree
		);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tsStringParser);
	timeTree.print();

	/*
	 * Ready to assemble time stepper
	 */
	sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes tssa;
	tssa.setup(pdeTerm_registry, timeStepper_registry);
	tssa.assembleTimeTree(
		timeTree,
		timeIntegrator
	);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tssa);

	timeIntegrator->setupConfigAndForwardTimeStepperEval(
			deSolver_Config,
			&evalFun
		);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeIntegrator);

	return true;
}

bool TimeTree::runIntegration(
		const sweet::Data::GenericContainer::Base &i_U,
		sweet::Data::GenericContainer::Base &o_U,
		double i_simulationTime
)
{
	return (timeIntegrator.get()->*evalFun)(i_U, o_U, i_simulationTime);
}

void TimeTree::clear()
{
	pdeTerm_registry.clear();
	timeStepper_registry.clear();

	timeIntegrator.reset();
}


TimeTree::~TimeTree()
{
}


}}
