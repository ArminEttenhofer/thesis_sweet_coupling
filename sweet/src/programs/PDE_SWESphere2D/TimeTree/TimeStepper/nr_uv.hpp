#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_NR_UV_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_NR_UV_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../DataContainer/Config.hpp"
#include "../../DataContainer/Simulation.hpp"
#include "../../Shack.hpp"

/*
 * Time tree node related includes
 */

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


class nr_uv	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<nr_uv>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	const sweet::Data::Sphere2D::Operators *_ops;

public:
	nr_uv();
	~nr_uv();
	nr_uv(
			const nr_uv &i_val
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
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
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
