#ifndef PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_LN_HPP
#define PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_LN_HPP


/*
 * Generic includes
 */
#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
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


class ln	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<ln>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWECart2D;
	sweet::Data::Cart2D::Shack *_shackCart2DDataOps;
	const sweet::Data::Cart2D::Operators *_ops;

	////*
	/// * Coriolis effect
	/// */
	///sweet::Data::Cart2D::DataGrid _fg;

public:
	ln();
	~ln();

	ln(
			const ln &i_val
	);

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

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
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;
};

}}}

#endif
