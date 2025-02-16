#ifndef PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_LG_HPP
#define PROGRAMS_PDE_SWECART2D_TIMETREE_TIMESTEPPER_LG_HPP


/*
 * Generic includes
 */
#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../Shack.hpp"
#include <sweet/ExpIntegration/Shack.hpp>

#include <programs/PDE_SWECart2D/NormalModes.hpp>

/*
 * Time tree node related includes
 */
#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"

namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {

class lg	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<lg>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWECart2D;
	sweet::Data::Cart2D::Shack *_shackCart2DDataOps;
	sweet::ExpIntegration::Shack *_shackExpIntegration;
	const sweet::Data::Cart2D::Operators *_ops;
	const sweet::Data::Cart2DComplex::Operators *_opsComplex;

	sweet::ExpIntegration::ExpFunction<double> _expFunction;
	PDE_SWECart2D::NormalModes::NormalModes pdeSWECart2DNormalModes;

	/*!
	 * Complex-valued time step size for complex-valued backward Euler
	 *
	 * This is used for REXI solvers of the form
	 * 	U_1 = (I-dt*L)^{-1} U_0
	 */
	std::complex<double> _rexiTermAlpha;
	std::complex<double> _rexiTermBeta;
	std::complex<double> _rexiTermGamma;
	bool _rexiTermGammaActive;

	std::vector<std::vector<std::array<std::array<std::complex<double>, 3>, 3>>> Z;  // Z[k1][k2][0,1,2][0,1,2];
	double _dt_precompute_phin;

public:
	lg();
	~lg();

	lg(const lg &i_val);

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

	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override;

	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
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

private:
	/*
	 * Return evaluation of backward Euler time step with complex data
	 */
	bool _eval_rexiTerm(
			const sweet::Data::GenericContainer::Base &i_U_,
			sweet::Data::GenericContainer::Base &o_U_,
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
