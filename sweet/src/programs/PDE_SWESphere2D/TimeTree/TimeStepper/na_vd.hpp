#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_NA_VD_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_NA_VD_HPP


#include <vector>
#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include <sweet/SemiLagrangian/Sphere2D.hpp>

#include "../../DataContainer/Config.hpp"
#include "../../DataContainer/Simulation.hpp"
#include "../../Shack.hpp"


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


/*!
 * Provide the nonlinear advection part using vorticity and divergence as advected variables.
 */
class na_vd	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<na_vd>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	sweet::Data::Sphere2D::Shack *_shackSphere2DDataOps;
	
	const sweet::Data::Sphere2D::Operators *_ops;

	class _SemiLagrangian
	{
	public:
		_SemiLagrangian()	:
			shack(nullptr),
			helper(nullptr),
			order(-1)
		{}

		sweet::SemiLagrangian::Shack *shack;
		sweet::SemiLagrangian::Sphere2D *helper;
		int order;

	} _semiLagrangian;

public:
	na_vd();
	~na_vd();
	na_vd(
			const na_vd &i_val
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

	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	)	override;

	void clear() override;

	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;

	bool evalNA_getNumStates(
			int *o_numStates
	)	override;

	bool evalNA_departurePoints(
			const sweet::Data::GenericContainer::Base* i_states[],	//!< Vector of states
			double i_timestepSize,
			sweet::Data::GenericContainer::Base &o_departurePositions		//!< Computed departure positions
	) override;

	bool evalNA_interpolate(
			const sweet::Data::GenericContainer::Base &i_U_input,		//!< Input simulation data
			const sweet::Data::GenericContainer::Base &i_samplingPositions,	//!< Sampling positions (computed by _evalNA_departurePoints)
			sweet::Data::GenericContainer::Base &o_U_samples			//!< Output samples
	) override;

};

}}}

#endif
