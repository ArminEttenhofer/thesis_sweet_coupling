#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_N_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_N_HPP


#include <sweet/Error/Base.hpp>
#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"
#include "../../Shack.hpp"

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


class n	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<n>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	const sweet::Data::Sphere2D::Operators *_ops;

public:
	n();

	~n();

	n(
			const n &i_val
	);

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

	const std::vector<std::string> getNodeNames() override;

public:
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override;

	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,	//!< Solver-specific configuration
		TIME_STEPPER_TYPES i_evalType,						//!< Choose the kind of evaluation (forward Euler, exponential, etc.)
													//!< for which this node should be used for.
		TimeTree_Node_Base::EvalFun *o_timeStepper	//!< Return the evaluation function to call the time integrator.
													//!< In this case, only _eval_tendencies can be returned and
													//!< i_evalType has to be EVAL_TENDENCIES.
	) override;

	void clear() override;

	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;
};

}}}

#endif
