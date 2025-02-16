#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_NR_VD_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_NR_VD_HPP


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


class nr_vd	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<nr_vd>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	const sweet::Data::Sphere2D::Operators *_ops;

public:
	nr_vd();
	~nr_vd();

	nr_vd(
			const nr_vd &i_val
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

private:
	void euler_timestep_update_na(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_phi_t,
			sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
			sweet::Data::Sphere2D::DataSpectral &o_div_t,

			double i_simulation_timestamp
	);

	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;
};

}}}

#endif
