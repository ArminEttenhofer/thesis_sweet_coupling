#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include <sweet/Error/Fatal.hpp>
#include "visc.hpp"

#include <complex>
#include <cmath>
#include <vector>

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


visc::visc(const std::string &i_viscType)	:
	_shackPDESWESphere2D(nullptr),
	_shackSphere2DDataOps(nullptr),
	_ops(nullptr),
	_opsComplex(nullptr),

	_viscOrder(2),
	_viscosity(0),
	_viscosityHighestModeNormalized(0),
	_alwaysNegative(true),

	_timeStepDependingViscosityFactor(0)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_EXPONENTIAL);

	if (i_viscType == "visc")
	{
		_viscType = VISC_ALL;
	}
	else if (i_viscType == "visc_phi_pert")
	{
		_viscType = VISC_PHI_PERT;
	}
	else if (i_viscType == "visc_vrt")
	{
		_viscType = VISC_VRT;
	}
	else if (i_viscType == "visc_div")
	{
		_viscType = VISC_DIV;
	}
	else
	{
		SWEETErrorFatal("Wrong visc type '"+i_viscType+"'");
	}
}


visc::~visc()
{
}


visc::visc(
		const visc &i_src
)	:
	TimeTree_Node_LeafHelper(i_src)
{
	_shackSphere2DDataOps = i_src._shackSphere2DDataOps;
	_shackPDESWESphere2D = i_src._shackPDESWESphere2D;
	_ops = i_src._ops;
	_opsComplex = i_src._opsComplex;

	_viscType = i_src._viscType;

	_viscOrder = i_src._viscOrder;
	_viscosity = i_src._viscosity;
	_viscosityHighestModeNormalized = i_src._viscosityHighestModeNormalized;

	_timeStepDependingViscosityFactor = i_src._timeStepDependingViscosityFactor;
	_alwaysNegative = i_src._alwaysNegative;
}


bool visc::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> visc::getNodeNames()
{
	std::vector<std::string> retval;

	if (_viscType == VISC_ALL)
	{
		retval.push_back("visc");
	}
	else if (_viscType == VISC_PHI_PERT)
	{
		retval.push_back("viscPhiPert");
	}
	else if (_viscType == VISC_VRT)
	{
		retval.push_back("viscVrt");
	}
	else if (_viscType == VISC_DIV)
	{
		retval.push_back("viscDiv");
	}
	else
	{
		SWEETErrorFatal("Internal error");
	}

	return retval;

}


bool visc::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Viscosity (";

	switch(_viscType)
	{
		case VISC_ALL:		o_ostream << "all";		break;
		case VISC_PHI_PERT:	o_ostream << "phi_pert";		break;
		case VISC_VRT:		o_ostream << "vrt";		break;
		case VISC_DIV:		o_ostream << "div";		break;
	}

	o_ostream << ")':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  WARNING: This can be also used for hyperviscosity." << std::endl;
	o_ostream << i_prefix << "  WARNING: In case of hyperviscosity and depending on the order, the minus signs" << std::endl;
	o_ostream << i_prefix << "  WARNING: are automatically swapped so that hyperviscosity acts always similar" << std::endl;
	o_ostream << i_prefix << "  WARNING: to viscosity, hence never amplifying." << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Parameters:" << std::endl;
	o_ostream << i_prefix << "    - order=[int]" << std::endl;
	o_ostream << i_prefix << "        Order of hyperviscosity (2, 4, 6, 8, ...)." << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "    - visc=[float]" << std::endl;
	o_ostream << i_prefix << "        Value of hyperviscosity." << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "    - nvisc=[float]" << std::endl;
	o_ostream << i_prefix << "        Normalized viscosity" << std::endl;
	o_ostream << i_prefix << "        (Normalized = highest mode damps with factor -1)" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "    - an=[bool]" << std::endl;
	o_ostream << i_prefix << "    - alwaysNegative=[bool]" << std::endl;
	o_ostream << i_prefix << "        Ensure that (hyper)viscosity in spectral space is always negative (default)" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "    - tsvisc=[float]" << std::endl;
	o_ostream << i_prefix << "        Time step dependent viscosity. Scales the viscosity by 'ts / thisValue'" << std::endl;

	return true;
}


bool visc::_setupInternals()
{
	if (_viscOrder != 2 && _viscOrder != 4 && _viscOrder != 6 && _viscOrder != 8)
		return error.set("Only orders 2, 4, 6, 8 are supported"+getNewLineDebugMessage());

	if (_viscosityHighestModeNormalized != 0)
	{
		double inv_r2 = 1.0/(_shackSphere2DDataOps->sphere2d_radius*_shackSphere2DDataOps->sphere2d_radius);

		int n = _ops->sphere2DDataConfig->spectral_modes_n_max;

		double D = 0;

		if (_viscOrder == 2)
		{
			D = (double)n*((double)n+1)*inv_r2;
		}
		else if (_viscOrder == 4)
		{
			D = (double)n*((double)n+1)*inv_r2;
			D = D*D;
		}
		else if (_viscOrder == 6)
		{
			D = (double)n*((double)n+1)*inv_r2;
			D = D*D*D;
		}
		else if (_viscOrder == 8)
		{
			D = (double)n*((double)n+1)*inv_r2;
			D = D*D*D*D;
		}

		SWEET_ASSERT(D >= 0);

		_viscosity = _viscosityHighestModeNormalized/D;
	}

	if (_timeStepDependingViscosityFactor != 0)
	{
		_viscosity *= _dt/_timeStepDependingViscosityFactor;

		double dt2 = _dt*_dt;

		if (_viscOrder == 2)
			_viscosity /= dt2;
		else if (_viscOrder == 4)
			_viscosity /= dt2*dt2;
		else if (_viscOrder == 6)
			_viscosity /= dt2*dt2*dt2;
		else if (_viscOrder == 8)
			_viscosity /= dt2*dt2*dt2*dt2;
	}

	return true;
}

bool visc::setupTreeNodeByFunction(
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
			return error.set("Functions are not supported here"+a->getNewLineDebugMessage());

		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
			return error.set("Values are not supported here"+a->getNewLineDebugMessage());
			break;

		case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
			if (a->key == "order" || a->key == "o")
			{
				a->getValue(_viscOrder);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
				break;
			}

			if (a->key == "visc" || a->key == "viscosity" || a->key == "v")
			{
				a->getValue(_viscosity);
				break;
			}

			if (a->key == "nvisc" || a->key == "normalizedViscosity")
			{
				a->getValue(_viscosityHighestModeNormalized);
				break;
			}

			if (a->key == "tsvisc" || a->key == "timeSteppingViscosity")
			{
				a->getValue(_timeStepDependingViscosityFactor);
				break;
			}

			if (a->key == "alwaysNegative" || a->key == "an")
			{
				a->getValue(_alwaysNegative);
				break;
			}

			return error.set("Key not supported"+a->getNewLineDebugMessage());
			break;

		default:
			SWEETErrorFatal("Internal error");
			return error.set("Internal error");
		}
	}

	// provide debug message in case that something goes wrong with the arguments
	setDebugMessage(i_function->getDebugMessage());

	return true;
}



bool visc::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun *o_timeStepper
)
{
	const DataContainer::Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;
	_opsComplex = myConfig.opsComplex;


	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_opsComplex != nullptr);

	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);

	if (i_evalType == EVAL_TENDENCIES)
	{
	}
	else if (i_evalType == EVAL_EULER_BACKWARD)
	{
	}
	else if (i_evalType == EVAL_EXPONENTIAL)
	{
		if (_expFunction.functionName == "")
		{
			// set default to phi0
			_expFunction.setup("phi0");
		}

		_expFunction.eval(std::complex<double>(0,0), expSteadyStateValue);
	}
	else
	{
		SWEETErrorFatal("Evaluation not supported - this error should have been caught earlier");
	}
	
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void visc::clear()
{
	TimeTree_Node_LeafHelper::clear();
}

/*!
 * Helper class.
 *
 * We use templated functions here since these implementations
 * can be SIMDizised (without if branchings)
 */
class ViscHelperClass
{
	/*
	 * Laplacian in spectral space:
	 *
	 * \nabla^2 = -n*(n+1)/(r^2)
	 */

	/*!
	 * Helper function to evaluate tendencies
	 *
	 * We use a templated function here since this implementation allows
	 * to be SIMDizised (without if branchings)
	 */
public:
	template <int order>
	static
	void eval_tendencies(
			sweet::Data::Sphere2D::DataSpectral &io_u,
			double i_inv_r2,
			double i_viscosity,
			bool i_alwaysNegative
	) {
		if (i_alwaysNegative)
		{
			/*
			 * Warning: This doesn't follow the standard hyperviscosity,
			 * but allows to simply increase the order and to get the similar effect
			 * without worring about the minus sign.
			 */
			io_u.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &io_data)
				{
					double D;
					if (order == 2)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
					}
					else if (order == 4)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D;
					}
					else if (order == 6)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D*D;
					}
					else if (order == 8)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D*D*D;
					}

					SWEET_ASSERT(D <= 0);

					io_data *= D*i_viscosity;
				}
			);
		}
		else
		{
			io_u.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &io_data)
				{
					double D;
					if (order == 2)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
					}
					else if (order == 4)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D;
					}
					else if (order == 6)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D*D;
					}
					else if (order == 8)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D*D*D;
					}

					io_data *= D*i_viscosity;
				}
			);
		}
	}

	/*!
	 * Helper function to evaluate backward Euler
	 *
	 * U1 = (I - dt*L)^{-1} U0
	 */
public:
	template <int order>
	static
	void eval_eulerBackward(
			sweet::Data::Sphere2D::DataSpectral &io_u,
			double i_inv_r2,
			double i_viscosity_dt,
			bool i_alwaysNegative
	) {
		if (i_alwaysNegative)
		{
			/*
			 * Warning: This doesn't follow the standard hyperviscosity,
			 * but allows to simply increase the order and to get the similar effect
			 * without worring about the minus sign.
			 */
			io_u.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &io_data)
				{
					double D;
					if (order == 2)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
					}
					else if (order == 4)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D;
					}
					else if (order == 6)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D*D;
					}
					else if (order == 8)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D*D*D;
					}

					/*
					 * (I - dt*L)^{-1} U0
					 */
					SWEET_ASSERT(1.0/(1.0 - i_viscosity_dt*D) <= 1.0);

					io_data = io_data/(1.0 - i_viscosity_dt*D);
				}
			);
		}
		else
		{
			io_u.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &io_data)
				{
					double D;
					if (order == 2)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
					}
					else if (order == 4)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D;
					}
					else if (order == 6)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D*D;
					}
					else if (order == 8)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D*D*D;
					}

					io_data = io_data/(1.0 - i_viscosity_dt*D);
				}
			);
		}
	}

	/*!
	 * Helper function to evaluate exponential
	 *
	 * U1 = (I - dt*L)^{-1} U0
	 */
public:
	template <int order>
	static
	void eval_exponential(
			sweet::Data::Sphere2D::DataSpectral &io_u,
			double i_inv_r2,
			double i_viscosity_dt,
			bool i_alwaysNegative
	) {
		if (i_alwaysNegative)
		{
			/*
			 * Warning: This doesn't follow the standard hyperviscosity,
			 * but allows to simply increase the order and to get the similar effect
			 * without worring about the minus sign.
			 */
			io_u.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &io_data)
				{
					double D;
					if (order == 2)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
					}
					else if (order == 4)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D;
					}
					else if (order == 6)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D*D;
					}
					else if (order == 8)
					{
						D = (double)n*((double)n+1)*i_inv_r2;
						D = -D*D*D*D;
					}

					io_data *= std::exp(i_viscosity_dt*D);
				}
			);
		}
		else
		{
			io_u.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &io_data)
				{
					double D;
					if (order == 2)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
					}
					else if (order == 4)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D;
					}
					else if (order == 6)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D*D;
					}
					else if (order == 8)
					{
						D = -(double)n*((double)n+1)*i_inv_r2;
						D = D*D*D*D;
					}

					io_data *= std::exp(i_viscosity_dt*D);
				}
			);
		}
	}
};


/*!
 * Simply set the time step size
 */
bool visc::setTimeStepSize(double i_dt)
{
	TimeTree_Node_LeafHelper::setTimeStepSize(i_dt);

	_setupInternals();

	return true;
}

bool visc::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);

	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	double inv_r2 = 1.0/(_shackSphere2DDataOps->sphere2d_radius*_shackSphere2DDataOps->sphere2d_radius);

	if (_viscType == VISC_PHI_PERT || _viscType == VISC_ALL)
	{
		o_U.phi_pert = i_U.phi_pert;

		if (_viscOrder == 2)
			ViscHelperClass::eval_tendencies<2>(o_U.phi_pert, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_tendencies<4>(o_U.phi_pert, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_tendencies<6>(o_U.phi_pert, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_tendencies<8>(o_U.phi_pert, inv_r2, _viscosity, _alwaysNegative);
	}
	else
	{
		o_U.phi_pert.spectral_setZero();
	}

	if (_viscType == VISC_VRT || _viscType == VISC_ALL)
	{
		o_U.vrt = i_U.vrt;

		if (_viscOrder == 2)
			ViscHelperClass::eval_tendencies<2>(o_U.vrt, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_tendencies<4>(o_U.vrt, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_tendencies<6>(o_U.vrt, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_tendencies<8>(o_U.vrt, inv_r2, _viscosity, _alwaysNegative);
	}
	else
	{
		o_U.vrt.spectral_setZero();
	}

	if (_viscType == VISC_DIV || _viscType == VISC_ALL)
	{
		o_U.div = i_U.div;

		if (_viscOrder == 2)
			ViscHelperClass::eval_tendencies<2>(o_U.div, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_tendencies<4>(o_U.div, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_tendencies<6>(o_U.div, inv_r2, _viscosity, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_tendencies<8>(o_U.div, inv_r2, _viscosity, _alwaysNegative);
	}
	else
	{
		o_U.div.spectral_setZero();
	}

	return true;
}



bool visc::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(i_U.phi_pert.sphere2DDataConfig->spectral_modes_n_max == _shackSphere2DDataOps->space_res_spectral[1]-1);


	/*
	 * Laplacian in spectral space:
	 *
	 * \nabla^2 = -n*(n+1)/(r^2)
	 */
	double inv_r2 = 1.0/(_shackSphere2DDataOps->sphere2d_radius*_shackSphere2DDataOps->sphere2d_radius);

	o_U.phi_pert = i_U.phi_pert;
	if (_viscType == VISC_PHI_PERT || _viscType == VISC_ALL)
	{
		if (_viscOrder == 2)
			ViscHelperClass::eval_eulerBackward<2>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_eulerBackward<4>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_eulerBackward<6>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_eulerBackward<8>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
	}

	o_U.vrt = i_U.vrt;
	if (_viscType == VISC_VRT || _viscType == VISC_ALL)
	{
		if (_viscOrder == 2)
			ViscHelperClass::eval_eulerBackward<2>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_eulerBackward<4>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_eulerBackward<6>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_eulerBackward<8>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
	}

	o_U.div = i_U.div;
	if (_viscType == VISC_DIV || _viscType == VISC_ALL)
	{
		if (_viscOrder == 2)
			ViscHelperClass::eval_eulerBackward<2>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_eulerBackward<4>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_eulerBackward<6>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_eulerBackward<8>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
	}

	return true;
}


bool visc::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(i_U.phi_pert.sphere2DDataConfig->spectral_modes_n_max == _shackSphere2DDataOps->space_res_spectral[1]-1);


	/*
	 * Laplacian in spectral space:
	 *
	 * \nabla^2 = -n*(n+1)/(r^2)
	 */
	double inv_r2 = 1.0/(_shackSphere2DDataOps->sphere2d_radius*_shackSphere2DDataOps->sphere2d_radius);

	o_U.phi_pert = i_U.phi_pert;
	if (_viscType == VISC_PHI_PERT || _viscType == VISC_ALL)
	{
		if (_viscOrder == 2)
			ViscHelperClass::eval_exponential<2>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_exponential<4>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_exponential<6>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_exponential<8>(o_U.phi_pert, inv_r2, _viscosity*_dt, _alwaysNegative);
	}

	o_U.vrt = i_U.vrt;
	if (_viscType == VISC_VRT || _viscType == VISC_ALL)
	{
		if (_viscOrder == 2)
			ViscHelperClass::eval_exponential<2>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_exponential<4>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_exponential<6>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_exponential<8>(o_U.vrt, inv_r2, _viscosity*_dt, _alwaysNegative);
	}

	o_U.div = i_U.div;
	if (_viscType == VISC_DIV || _viscType == VISC_ALL)
	{
		if (_viscOrder == 2)
			ViscHelperClass::eval_exponential<2>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 4)
			ViscHelperClass::eval_exponential<4>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 6)
			ViscHelperClass::eval_exponential<6>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
		else if (_viscOrder == 8)
			ViscHelperClass::eval_exponential<8>(o_U.div, inv_r2, _viscosity*_dt, _alwaysNegative);
	}

	return true;
}



}}}
