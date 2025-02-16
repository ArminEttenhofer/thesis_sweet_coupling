#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_LG_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_LG_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>
#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"
#include "../../Shack.hpp"
#include <sweet/ExpIntegration/Shack.hpp>

/*
 * Time tree node related includes
 */

namespace PDE_SWESphere2D {
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
	Shack *_shackPDESWESphere2D;
	sweet::Data::Sphere2D::Shack *_shackSphere2DDataOps;
	sweet::ExpIntegration::Shack *_shackExpIntegration;
	const sweet::Data::Sphere2D::Operators *_ops;
	const sweet::Data::Sphere2DComplex::Operators *_opsComplex;

	sweet::ExpIntegration::ExpFunction<double> _expFunction;

	//! Value if evaluating phi functions at 0
	std::complex<double> expSteadyStateValue;

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

	/*!
		Store evaluation of the direct exponential:
		Z = Q*e^{\Lambda*\Dt}*Q^{-1} = f(m,n,Dt)
		
		For each wavenumber pair (m,n), we store a 2x2 matrix
		These matrices need to be computed again only if the time step sizes is modified,
		which is checked by comparing the current Dt with _dt_precompute_phin
	*/
	std::vector<std::array<std::array<std::complex<double>, 2>, 2>> Z;  // Z[idx][0,1][0,1];
	std::vector<sweet::Data::Sphere2D::DataSpectral> exp_coef;
	double _dt_precompute_phin;

public:
	lg();
	~lg();

	lg(const lg &i_val);

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

public:
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override;

	const std::vector<std::string> getNodeNames() override;

	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override;

	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override;

	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	) override;

	bool setTimeStepSize(double i_dt)	override;

	void clear() override;


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

private:
	bool _computeExpDirectCoefficients();

};

}}}

#endif
