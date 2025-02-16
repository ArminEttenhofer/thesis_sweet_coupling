#include "TimeTreeIR.hpp"


#include "TimeStepper/ln.hpp"

#include "TimeStepper/l.hpp"
#include "TimeStepper/lg.hpp"
#include "TimeStepper/lc.hpp"

#include "TimeStepper/n.hpp"
#include "TimeStepper/na.hpp"
#include "TimeStepper/nr.hpp"

#include "TimeStepper/null.hpp"



namespace PDE_SWECart2D {
namespace TimeTree {


bool TimeTree::setup_1_registerAllTimesteppers()
{

	/*
	 * Registration of all possible PDE terms
	 */
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::ln>();

	pdeTerm_registry.registerTimeTreeNode<TimeStepper::l>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::lg>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::lc>();

	pdeTerm_registry.registerTimeTreeNode<TimeStepper::n>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::na>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::nr>();

	pdeTerm_registry.registerTimeTreeNode<TimeStepper::null>();

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeTerm_registry);

	/*
	 * Register time steppers
	 */
	sweet::TimeTree::InteriorNodes_Registry registryAll;
	registryAll.registerAll(timeStepper_registry);

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(registryAll);

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

bool TimeTree::setup_3_timestepper(
		const std::string &i_timestepping_method,	//!< String with time stepping method such as SS(ERK(lg,order=4),ERK(lc,order=2),order=2)
		sweet::Shacks::ProgramArgumentsDictionary *i_progArgShackDict,
		sweet::Data::Cart2D::Operators *io_ops,
		sweet::Data::Cart2DComplex::Operators *io_opsComplex,
		const DataContainer::Simulation &i_U
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

#if 0
	timeIntegrator->shackRegistration(i_progArgShackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeIntegrator);

	i_progArgShackDict->processProgramArguments();
#endif

	deSolver_Config.myDataContainer = &i_U;
	deSolver_Config.ops = io_ops;
	deSolver_Config.opsComplex = io_opsComplex;
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
