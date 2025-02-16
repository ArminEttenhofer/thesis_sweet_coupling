#include <sweet/Data/Cart2D/Operators.hpp>
#include "na.hpp"
#include <vector>



namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {


na::na()	:
	_shackPDESWECart2D(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_SEMI_LAGRANGIAN);
}

na::~na()
{
}


na::na(
		const na &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWECart2D = i_value._shackPDESWECart2D;
	_shackCart2DDataOps = i_value._shackCart2DDataOps;
	_semiLagrangian = i_value._semiLagrangian;
	_ops = i_value._ops;
}

bool na::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWECart2D = io_shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_semiLagrangian.shack = io_shackDict->getAutoRegistration<sweet::SemiLagrangian::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> na::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("na");
	return retval;

}


bool na::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun* o_timeStepper
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

	if (_evalTypeRequested == EVAL_TENDENCIES)
	{
	}
	else if (_evalTypeRequested == EVAL_SEMI_LAGRANGIAN)
	{
		if (_semiLagrangian.order != 1 && _semiLagrangian.order != 2)
		{
			std::ostringstream oss;
			oss << "na DE term: Invalid time stepping order " << _semiLagrangian.order;
			return error.set(oss.str());
		}

		SWEET_ASSERT(_semiLagrangian.helper == nullptr);

		_semiLagrangian.helper = new sweet::SemiLagrangian::Cart2D;
		_semiLagrangian._sampler2D = new sweet::Data::Cart2D::DataSampler;
		setup_SL();
	}
	else
	{
		SWEETErrorFatal("Shouldn't happen :-(");
	}


	return true;
}

bool na::setupTreeNodeByFunction(
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




void na::clear()
{
	TimeTree_Node_LeafHelper::clear();

	if (_semiLagrangian.helper != nullptr)
	{
		delete _semiLagrangian.helper;
		_semiLagrangian.helper = nullptr;
	}

	if (_semiLagrangian._sampler2D != nullptr)
	{
		delete _semiLagrangian._sampler2D;
		_semiLagrangian._sampler2D = nullptr;
	}

}

/*
 * Return the time tendencies of the PDE term
 */
bool na::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	// TODO
	SWEETErrorFatal("Not implemented (yet!)");

	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_U.u = -i_U.u*_ops->diff_c_x(i_U.u) - i_U.v*_ops->diff_c_y(i_U.u);
	o_U.v = -i_U.u*_ops->diff_c_x(i_U.v) - i_U.v*_ops->diff_c_y(i_U.v);

	o_U.h_pert = - (i_U.u*_ops->diff_c_x(i_U.h_pert) + i_U.v*_ops->diff_c_y(i_U.h_pert));

	return true;
}

bool na::evalNA_getNumStates(
		int *o_numStates
)
{
	*o_numStates = 2;
	return true;
}


bool na::evalNA_departurePoints(
		const sweet::Data::GenericContainer::Base* i_U[],	//!< Vector of states
		double i_timestepSize,
		sweet::Data::GenericContainer::Base &o_departurePositions		//!< Computed departure positions
)
{

	///SWEETErrorFatal("Not implemented (yet!)");

	const DataContainer::Simulation &U0 = cast(*i_U[1]);
	const DataContainer::Simulation &U1 = cast(*i_U[0]);

	const sweet::Data::Cart2D::DataSpectral &U0_u = U0.u;
	const sweet::Data::Cart2D::DataSpectral &U0_v = U0.v;

	const sweet::Data::Cart2D::DataSpectral &U1_u = U1.u;
	const sweet::Data::Cart2D::DataSpectral &U1_v = U1.v;

	DataContainer::SemiLagPositions &departurePoints = static_cast<DataContainer::SemiLagPositions&>(o_departurePositions);

	sweet::Data::Vector::Vector<double> &departurePointsX = departurePoints.pos_x;
	sweet::Data::Vector::Vector<double> &departurePointsY = departurePoints.pos_y;
	departurePointsX = _semiLagrangian.posx_a;
	departurePointsY = _semiLagrangian.posy_a;

	sweet::Data::Cart2D::Staggering staggering;
	SWEET_ASSERT(staggering.staggering_type == 'a');

	// Calculate departure points
	_semiLagrangian.helper->semi_lag_departure_points_settls(
			U0_u.toGrid(),		U0_v.toGrid(),
			U1_u.toGrid(),		U1_v.toGrid(),
			_semiLagrangian.posx_a,			_semiLagrangian.posy_a,
			i_timestepSize,
			departurePointsX,	departurePointsY,
			_shackCart2DDataOps->cart2d_domain_size,
			&staggering,
			_semiLagrangian.order,

			_semiLagrangian.shack->semi_lagrangian_max_iterations,
			_semiLagrangian.shack->semi_lagrangian_convergence_threshold

	);

	return true;
}

bool na::evalNA_interpolate(
		const sweet::Data::GenericContainer::Base &i_U,		//!< Input simulation data
		const sweet::Data::GenericContainer::Base &i_samplingPositions,	//!< Sampling positions (computed by _evalNA_departurePoints)
		sweet::Data::GenericContainer::Base &o_U_sampled			//!< Output samples
)
{

	///SWEETErrorFatal("Not implemented (yet!)");

	const DataContainer::Simulation &U = cast(i_U);
	const sweet::Data::Cart2D::DataSpectral &U_h_pert = U.h_pert;
	const sweet::Data::Cart2D::DataSpectral &U_u = U.u;
	const sweet::Data::Cart2D::DataSpectral &U_v = U.v;

	DataContainer::Simulation &U_sampled = cast(o_U_sampled);
	sweet::Data::Cart2D::DataSpectral &U_sampled_h_pert = U_sampled.h_pert;
	sweet::Data::Cart2D::DataSpectral &U_sampled_u = U_sampled.u;
	sweet::Data::Cart2D::DataSpectral &U_sampled_v = U_sampled.v;

	const DataContainer::SemiLagPositions &samplingPoints = static_cast<const DataContainer::SemiLagPositions&>(i_samplingPositions);

	const sweet::Data::Vector::Vector<double> &samplingPointsX = samplingPoints.pos_x;
	const sweet::Data::Vector::Vector<double> &samplingPointsY = samplingPoints.pos_y;

	const sweet::Data::Cart2D::DataGrid U_h_pert_phys = U_h_pert.toGrid();
	const sweet::Data::Cart2D::DataGrid U_u_phys = U_u.toGrid();
	const sweet::Data::Cart2D::DataGrid U_v_phys = U_v.toGrid();

	_semiLagrangian._sampler2D->bicubic_scalar(U_h_pert_phys, samplingPointsX, samplingPointsY, U_sampled_h_pert, -0.5, -0.5);
	_semiLagrangian._sampler2D->bicubic_scalar(U_u_phys, samplingPointsX, samplingPointsY, U_sampled_u, -0.5, -0.5);
	_semiLagrangian._sampler2D->bicubic_scalar(U_v_phys, samplingPointsX, samplingPointsY, U_sampled_v, -0.5, -0.5);

	return true;
}


void na::setup_SL()
{
	_semiLagrangian.posx_a.setup(_ops->cart2DDataConfig->grid_number_elements);
	_semiLagrangian.posy_a.setup(_ops->cart2DDataConfig->grid_number_elements);

	// Setup sampler for future interpolations
	_semiLagrangian._sampler2D->setup(_shackCart2DDataOps->cart2d_domain_size, _ops->cart2DDataConfig);

	// Setup semi-lag
	_semiLagrangian.helper->setup(_shackCart2DDataOps->cart2d_domain_size, _ops->cart2DDataConfig);

	sweet::Data::Cart2D::DataGrid tmp_x(_ops->cart2DDataConfig);
	tmp_x.grid_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)i)*_shackCart2DDataOps->cart2d_domain_size[0]/(double)_shackCart2DDataOps->space_res_physical[0];
			},
			false
	);
	sweet::Data::Cart2D::DataGrid tmp_y(_ops->cart2DDataConfig);
	tmp_y.grid_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)j)*_shackCart2DDataOps->cart2d_domain_size[1]/(double)_shackCart2DDataOps->space_res_physical[1];
			},
			false
	);

	//Initialize arrival points with h position
	sweet::Data::Vector::Vector<double> pos_x = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(tmp_x);
	sweet::Data::Vector::Vector<double> pos_y = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(tmp_y);

	double cell_size_x = _shackCart2DDataOps->cart2d_domain_size[0]/(double)_shackCart2DDataOps->space_res_physical[0];
	double cell_size_y = _shackCart2DDataOps->cart2d_domain_size[1]/(double)_shackCart2DDataOps->space_res_physical[1];

	// Initialize arrival points with h position
	_semiLagrangian.posx_a = pos_x+0.5*cell_size_x;
	_semiLagrangian.posy_a = pos_y+0.5*cell_size_y;

}

}}}
