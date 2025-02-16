#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_TIMETREE_DE_TERMS_LGENERIC_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_TIMETREE_DE_TERMS_LGENERIC_HPP


/*
 * Generic includes
 */
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
 * d/dt u(t) = \lambda_1 * u(t) + \lambda_2 * u(t) + \lambda_3 * u(t) + mu * sin(t)
 */
class lambdaGeneric	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<lambdaGeneric>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *shackODEGeneric_DE_Dahlquist;

private:
	std::string _termIdStr;

private:
	std::complex<double> _lambda;

	sweet::ExpIntegration::ExpFunction<double> _expFunction;

	std::complex<double> _rexiTermAlpha;
	std::complex<double> _rexiTermBeta;
	std::complex<double> _rexiTermGamma;
	bool _rexiTermGammaActive;

	int _semiLagrangian_order;

public:
	lambdaGeneric(const std::string &i_termId);

	~lambdaGeneric();

	lambdaGeneric(
			const lambdaGeneric &i_val
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

	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override;

	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	) override;

	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	) override;

	bool setTimeStepSize(double i_dt)	override;

private:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;

private:
	bool _eval_exponential(
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
	bool _eval_rexiTerm(
			const sweet::Data::GenericContainer::Base &i_U_,
			sweet::Data::GenericContainer::Base &o_U_,
			double i_timeStamp
	)	override;

public:
	bool evalNA_getNumStates(
			int *o_numStates
	)	override;

public:
	bool evalNA_departurePoints(
			const sweet::Data::GenericContainer::Base* i_states[],	//!< Vector of states
			double i_timestepSize,
			sweet::Data::GenericContainer::Base &o_departurePositions		//!< Computed departure positions
	) override;

public:
	bool evalNA_interpolate(
			const sweet::Data::GenericContainer::Base &i_U_input,		//!< Input simulation data
			const sweet::Data::GenericContainer::Base &i_samplingPositions,	//!< Sampling positions (computed by _evalNA_departurePoints)
			sweet::Data::GenericContainer::Base &o_U_samples			//!< Output samples
	) override;
};

}}}

#endif
