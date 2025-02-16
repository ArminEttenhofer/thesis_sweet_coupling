#ifndef PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_NA_HPP
#define PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_NA_HPP


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
#include "../../DataContainer/SemiLagPositions.hpp"
#include "../../DataContainer/Config.hpp"

#include <sweet/SemiLagrangian/Cart2D.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>

namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {

class na	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<na>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWECart2D;
	sweet::Data::Cart2D::Shack *_shackCart2DDataOps;
	const sweet::Data::Cart2D::Operators *_ops;

	////sweet::SemiLagrangian::Cart2D _semiLagrangian;
	////sweet::Data::Cart2D::DataSampler _sampler2D;
	////sweet::Data::Vector::Vector<double> posx_a, posy_a;

	class _SemiLagrangian
	{
	public:
		_SemiLagrangian()	:
			shack(nullptr),
			helper(nullptr),
			_sampler2D(nullptr),
			order(-1)
		{}

		sweet::SemiLagrangian::Shack *shack;
		sweet::SemiLagrangian::Cart2D *helper;
		sweet::Data::Cart2D::DataSampler *_sampler2D;
		sweet::Data::Vector::Vector<double> posx_a, posy_a;
		int order;

	} _semiLagrangian;



public:
	na();
	~na();
	na(
			const na &i_val
	);

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

	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	)	override;


	void clear() override;

	/*
	 * Return the time tendencies of the PDE term
	 */
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

	void setup_SL();
};

}}}

#endif
