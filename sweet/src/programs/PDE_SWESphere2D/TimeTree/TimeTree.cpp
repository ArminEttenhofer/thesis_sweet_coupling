#include <programs/PDE_SWESphere2D/TimeTree/TimeTree.hpp>
#include "TimeStepper/ln.hpp"

#include "TimeStepper/l.hpp"
#include "TimeStepper/lg.hpp"
#include "TimeStepper/lc.hpp"
#include "TimeStepper/lb.hpp"

#include "TimeStepper/n.hpp"
#include "TimeStepper/na_uv.hpp"
#include "TimeStepper/nr_uv.hpp"
#include "TimeStepper/na_vd.hpp"
#include "TimeStepper/nr_vd.hpp"

#include "TimeStepper/visc.hpp"

#include "TimeStepper/null.hpp"



namespace PDE_SWESphere2D {
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
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::lb>();

	pdeTerm_registry.registerTimeTreeNode<TimeStepper::n>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::na_uv>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::nr_uv>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::na_vd>();
	pdeTerm_registry.registerTimeTreeNode<TimeStepper::nr_vd>();

	pdeTerm_registry.registerTimeTreeNode<TimeStepper::visc>();
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<TimeStepper::visc>("visc_phi_pert")
	);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<TimeStepper::visc>("visc_vrt")
	);
	pdeTerm_registry.registerTimeTreeNode(
			std::make_shared<TimeStepper::visc>("visc_div")
	);

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

	i_shackDict->getAutoRegistration(&shackParallelization);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i_shackDict);

	////i_shackDict->getAutoRegistration(&shackBenchmarks);
	////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i_shackDict);

	return true;
}

bool TimeTree::setup_3_timestepper(
		const std::string &i_timesteppingMethod,	//!< String with time stepping method such as SS(ERK(lg,order=4),ERK(lc,order=2),order=2)
		sweet::Shacks::ProgramArgumentsDictionary *i_progArgShackDict,
		sweet::Data::Sphere2D::Operators *io_ops,
		sweet::Data::Sphere2DComplex::Operators *io_opsComplex,
		const DataContainer::Simulation &i_U
)
{
	/*
	 * Special handling of "help"
	 */
	if (i_timesteppingMethod == "help" || i_timesteppingMethod == "helpall")
		return false;

	/*
	 * Setup time stepping string parser and parse it
	 */
	sweet::TimeTree::TimeString_2_TimeTreeIR tsStringParser;
	tsStringParser.createTimeTree(
			i_timesteppingMethod,
			timeTreeIR
		);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tsStringParser);
	timeTreeIR.print();

	/*
	 * Ready to assemble time stepper
	 */
	sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes tssa;
	tssa.setup(pdeTerm_registry, timeStepper_registry);
	tssa.assembleTimeTree(
		timeTreeIR,
		timeIntegrator
	);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tssa);


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
	timeTreeIR.clear();

	timeIntegrator.reset();
}

void TimeTree::print()
{
	timeTreeIR.print();
}

TimeTree::~TimeTree()
{
}

bool TimeTree::outputHelp(
		std::ostream &o_ostringstream,
		const std::string &i_prefix,
		int i_verbosity
)
{
	std::string newPrefix = i_prefix + " | ";
	std::cout << "************************************************************" << std::endl;
	std::cout << "* TIME TREE: Interior node registry (Time steppers + tools)" << std::endl;
	std::cout << "************************************************************" << std::endl;
	timeStepper_registry.outputHelp(o_ostringstream, newPrefix, i_verbosity);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeStepper_registry);

	std::cout << "************************************************************" << std::endl;
	std::cout << "* TIME TREE: Leaf nodes registry (PDE terms)" << std::endl;
	std::cout << "************************************************************" << std::endl;
	pdeTerm_registry.outputHelp(o_ostringstream, newPrefix, i_verbosity);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeTerm_registry);

	return true;
}


}}
