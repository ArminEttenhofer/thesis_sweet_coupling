#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_VISC_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_VISC_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
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

class visc	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<visc>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	sweet::Data::Sphere2D::Shack *_shackSphere2DDataOps;
	const sweet::Data::Sphere2D::Operators *_ops;
	const sweet::Data::Sphere2DComplex::Operators *_opsComplex;


	//! Order of (hyper)viscosity
	int _viscOrder;

	//! Viscosity
	double _viscosity;

	//! Viscosity based on normalization
	double _viscosityHighestModeNormalized;

	//! If true, we ensure that (hyper)viscosity always leads to a decreasing effect
	bool _alwaysNegative;

	//! Scale the viscosity in addition with " ts / '_timeStepDependingViscosityFactor'"
	//! This allows to have the viscosity vanish for smaller time step sizes (=> convergence)
	double _timeStepDependingViscosityFactor;

	//! Exponential function to be used
	sweet::ExpIntegration::ExpFunction<double> _expFunction;

	//! Value if evaluating phi functions at 0
	std::complex<double> expSteadyStateValue;

	/*!
	 * Type of viscosity
	 */
	enum VISC_TYPE
	{
		VISC_ALL = 1,
		VISC_PHI_PERT,
		VISC_VRT,
		VISC_DIV
	};

	VISC_TYPE _viscType;

	/*!
	 * Constructor which allows to specify also the type of viscosity
	 *
	 * Different types are supported:
	 *
	 * 'visc': Viscosity is applied to all fields
	 * 'visc_phi_pert': Viscosity is applied only to the phi_pert field
	 * 'visc_vrt': Viscosity is only applied to the vorticity field
	 * 'visc_div': Viscosity is only applied to the divergence field
	 */
public:
	visc(
			const std::string &i_viscType = "visc"	///< type of viscosity
	);
	~visc();

	visc(const visc &i_val);

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

private:
	bool _setupInternals();

public:
	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	)	override;

	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override;

	void clear() override;


	/*!
	 * Simply set the time step size
	 */
	bool setTimeStepSize(double i_dt)	override;

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
