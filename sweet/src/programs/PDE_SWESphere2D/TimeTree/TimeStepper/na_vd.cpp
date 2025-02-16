#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/SemiLagPositions.hpp"
#include "na_vd.hpp"


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


na_vd::na_vd()	:
	_shackPDESWESphere2D(nullptr),
	_shackSphere2DDataOps(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_SEMI_LAGRANGIAN);
}

na_vd::~na_vd()
{
}

na_vd::na_vd(
		const na_vd &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWESphere2D = i_value._shackPDESWESphere2D;
	_shackSphere2DDataOps = i_value._shackSphere2DDataOps;
	_semiLagrangian = i_value._semiLagrangian;
	
	_ops = i_value._ops;
}

bool na_vd::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_semiLagrangian.shack = io_shackDict->getAutoRegistration<sweet::SemiLagrangian::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> na_vd::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("na_vd");
	return retval;

}



bool na_vd::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Non-linear advection term (using vorticity-divergence formulation)':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool na_vd::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun *o_timeStepper
)
{
	const DataContainer::Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	SWEET_ASSERT(_evalTypeRequested != EVAL_NONE);

	if (_evalTypeRequested == EVAL_TENDENCIES)
	{
	}
	else if (_evalTypeRequested == EVAL_SEMI_LAGRANGIAN)
	{
		if (_semiLagrangian.order != 1 && _semiLagrangian.order != 2)
		{
			std::ostringstream oss;
			oss << "na_vd DE term: Invalid time stepping order " << _semiLagrangian.order;
			return error.set(oss.str());
		}

		SWEET_ASSERT(_semiLagrangian.helper == nullptr);

		_semiLagrangian.helper = new sweet::SemiLagrangian::Sphere2D;
		_semiLagrangian.helper->setup(
				_ops->sphere2DDataConfig,
				_semiLagrangian.shack,
				_semiLagrangian.order
			);
	}
	else
	{
		SWEETErrorFatal("Shouldn't happen :-(");
	}
	return true;
}

void na_vd::clear()
{
	TimeTree_Node_LeafHelper::clear();

	if (_semiLagrangian.helper != nullptr)
	{
		delete _semiLagrangian.helper;
		_semiLagrangian.helper = nullptr;
	}
}


bool na_vd::setupTreeNodeByFunction(
		std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
		sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
)
{
	for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
	{
		sweet::TimeTree::TimeTreeIR::Argument *a = iter->get();

		switch(a->argType)
		{
		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
			error.set("Time steppers inside this time stepper are not allowed, yet"+a->getNewLineDebugMessage());
			return false;
			break;

		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
			if (a->key == "order" || a->key == "sl_order" || a->key == "o")
			{
				a->getValue(_semiLagrangian.order);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
				break;
			}

			return error.set("Key not supported"+a->getNewLineDebugMessage());
			break;

		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
			error.set("Just values as arguments are not supported, yet"+a->getNewLineDebugMessage());
			break;

		default:
			SWEETErrorFatal("Internal error");
			return error.set("Internal error");
		}
	}

	// provide debug message in case that something goes wrong with the arguments
	setDebugMessage(i_function->getDebugMessage());
	//return _setupArgumentInternals();
	return true;
}



bool na_vd::_eval_tendencies(
	const sweet::Data::GenericContainer::Base &i_U_,
	sweet::Data::GenericContainer::Base &o_U_,
	double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);

	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, U_u_phys, U_v_phys);

	sweet::Data::Sphere2D::DataGrid U_div_phys = i_U.div.toGrid();
	o_U.phi_pert = -_ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.phi_pert.toGrid());
	o_U.vrt = -_ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.vrt.toGrid());
	o_U.div = -_ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.div.toGrid());

	return true;
}

bool na_vd::evalNA_getNumStates(
		int *o_numStates
)
{
	*o_numStates = 2;
	return true;
}


bool na_vd::evalNA_departurePoints(
		const sweet::Data::GenericContainer::Base* i_U[],	//!< Vector of states
		double i_timestepSize,
		sweet::Data::GenericContainer::Base &o_departurePositions		//!< Computed departure positions
)
{
	/*
	 * Step 1) SL
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 *
	 * See Hortal's paper for equation.
	 */
	const DataContainer::Simulation &U0 = cast(*i_U[1]);
	const DataContainer::Simulation &U1 = cast(*i_U[0]);

	const sweet::Data::Sphere2D::DataSpectral &U0_vrt = U0.vrt;
	const sweet::Data::Sphere2D::DataSpectral &U0_div = U0.div;

	const sweet::Data::Sphere2D::DataSpectral &U1_vrt = U1.vrt;
	const sweet::Data::Sphere2D::DataSpectral &U1_div = U1.div;

	DataContainer::SemiLagPositions &departurePoints = static_cast<DataContainer::SemiLagPositions&>(o_departurePositions);

	sweet::Data::Vector::Vector<double> &departurePointsLon = departurePoints.lon;
	sweet::Data::Vector::Vector<double> &departurePointsLat = departurePoints.lat;



	sweet::Data::Sphere2D::DataGrid U_u_lon_prev, U_v_lat_prev;
	_ops->vrtdiv_2_uv(U0_vrt, U0_div, U_u_lon_prev, U_v_lat_prev);

	sweet::Data::Sphere2D::DataGrid U_u_lon, U_v_lat;
	_ops->vrtdiv_2_uv(U1_vrt, U1_div, U_u_lon, U_v_lat);

	double dt_div_radius = i_timestepSize / _shackSphere2DDataOps->sphere2d_radius;

	// Calculate departure points
	sweet::Data::Vector::Vector<double> _pos_lon_d, _pos_lat_d;
	_semiLagrangian.helper->semi_lag_departure_points_settls_specialized(
			// previous velocities
			dt_div_radius*U_u_lon_prev, dt_div_radius*U_v_lat_prev,
			// current velocities
			dt_div_radius*U_u_lon, dt_div_radius*U_v_lat,
			departurePointsLon, departurePointsLat		// OUTPUT
	);

	return true;
}

bool na_vd::evalNA_interpolate(
		const sweet::Data::GenericContainer::Base &i_U,		//!< Input simulation data
		const sweet::Data::GenericContainer::Base &i_samplingPositions,	//!< Sampling positions (computed by _evalNA_departurePoints)
		sweet::Data::GenericContainer::Base &o_U_sampled			//!< Output samples
)
{
	const DataContainer::Simulation &U = cast(i_U);
	const sweet::Data::Sphere2D::DataSpectral &U_phi_pert = U.phi_pert;
	const sweet::Data::Sphere2D::DataSpectral &U_vrt = U.vrt;
	const sweet::Data::Sphere2D::DataSpectral &U_div = U.div;

	DataContainer::Simulation &U_sampled = cast(o_U_sampled);
	sweet::Data::Sphere2D::DataSpectral &U_sampled_phi_pert = U_sampled.phi_pert;
	sweet::Data::Sphere2D::DataSpectral &U_sampled_vrt = U_sampled.vrt;
	sweet::Data::Sphere2D::DataSpectral &U_sampled_div = U_sampled.div;

	const DataContainer::SemiLagPositions &samplingPoints = static_cast<const DataContainer::SemiLagPositions&>(i_samplingPositions);

	sweet::Data::Vector::Vector<double> &samplingPointsLon = samplingPoints.lon;
	sweet::Data::Vector::Vector<double> &samplingPointsLat = samplingPoints.lat;

	_semiLagrangian.helper->apply_sl_timeintegration_vd(
				_ops,
				U_phi_pert, U_vrt, U_div,
				samplingPointsLon, samplingPointsLat,
				U_sampled_phi_pert, U_sampled_vrt, U_sampled_div
			);

	return true;
}

}}}
