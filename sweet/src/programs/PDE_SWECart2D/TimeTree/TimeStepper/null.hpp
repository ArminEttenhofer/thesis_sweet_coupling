#ifndef PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_NULL_HPP
#define PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_NULL_HPP


/*
 * Generic includes
 */
#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../Shack.hpp"

/*
 * Time tree node related includes
 */
#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"

namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {

class null	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<null>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWECart2D;
	sweet::Data::Cart2D::Shack *_shackCart2DDataOps;
	const sweet::Data::Cart2D::Operators *_ops;


public:
	null();
	~null();

	null(const null &i_val);

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

	virtual
	const std::vector<std::string> getNodeNames() override;

	virtual
	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun* o_timeStepper
	) override;

	void clear() override;


	/*
	 * Return the time tendencies of the PDE term
	 */
public:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	) override;

	/*
	 * Return the backward Euler time step
	 */
public:
	bool _eval_eulerBackward(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	) override;

	/*
	 * Compute an exponential integration for a given exp term
	 */
private:
	bool _eval_exponential(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	) override;
};

}}}

#endif
