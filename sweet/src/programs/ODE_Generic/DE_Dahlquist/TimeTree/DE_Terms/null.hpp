#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_TIMETREE_DE_TERMS_NULL_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_TIMETREE_DE_TERMS_NULL_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../DataContainer/Config.hpp"
#include "../../DataContainer/Simulation.hpp"
#include "../../Shack.hpp"


namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DE_Terms {

/**
 * This is one of the \lambda_i terms in solving the ODE
 *
 * d/dt u(t) = \lambda_1 * u(t) + \lambda_2 * u(t) + \lambda_3 * u(t) + \mu * sin(t)
 */
class null	:
		public sweet::TimeTree::TimeTree_Node_LeafHelper<null>,
		public sweet::Data::GenericContainer::CastHelper<
				DataContainer::Simulation,
				DataContainer::Config
			>
{
private:
	Shack *shackODEGeneric_DE_Dahlquist;

public:
	null();

	~null();

	null(
			const null &i_val
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
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override;

	void clear() override;

	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	) override;

	bool setTimeStepSize(double i_dt)	override;

private:
	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;

private:
	bool _eval_eulerBackward(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;

private:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;
};

}}}

#endif
