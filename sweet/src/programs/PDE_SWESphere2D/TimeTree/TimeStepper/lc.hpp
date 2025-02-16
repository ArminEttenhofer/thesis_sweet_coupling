#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_LC_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_LC_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"
#include "../../Shack.hpp"

/*
 * Time tree node related includes
 */

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {

class lc :
	public sweet::TimeTree::TimeTree_Node_LeafHelper<lc>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	const sweet::Data::Sphere2D::Operators *_ops;

	/*
	 * Coriolis effect
	 */
	sweet::Data::Sphere2D::DataGrid fg;

	/*
	 * Temporary variables
	 */
	sweet::Data::Sphere2D::DataGrid ug;
	sweet::Data::Sphere2D::DataGrid vg;

public:
	lc();

	~lc() override;

	lc(
			const lc &i_val
	);

private:
	void _setupDataBuffers();

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)	override;

public:
	const std::vector<std::string> getNodeNames()	override;

public:
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override;

public:
	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override;


	void clear() override;


private:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_timeStamp
	)	override;
};

}}}

#endif
