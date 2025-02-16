#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include "lg.hpp"

#include <complex>
#include <cmath>
#include <vector>

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


lg::lg()	:
	_shackPDESWESphere2D(nullptr),
	_shackSphere2DDataOps(nullptr),
	_ops(nullptr),
	_opsComplex(nullptr),
	_rexiTermAlpha(666,-666),
	_rexiTermBeta(666,-666),
	_rexiTermGamma(666,-666),
	_rexiTermGammaActive(false)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_EXPONENTIAL);
	setEvalAvailable(EVAL_REXI_TERM);
}


lg::~lg()
{
}


lg::lg(
		const lg &i_src
)	:
	TimeTree_Node_LeafHelper(i_src)
{
	_shackSphere2DDataOps = i_src._shackSphere2DDataOps;
	_shackPDESWESphere2D = i_src._shackPDESWESphere2D;
	_shackExpIntegration = i_src._shackExpIntegration;
	_ops = i_src._ops;
	_opsComplex = i_src._opsComplex;

	_rexiTermAlpha = i_src._rexiTermAlpha;
	_rexiTermBeta = i_src._rexiTermBeta;
	_rexiTermGamma = i_src._rexiTermGamma;
	_rexiTermGammaActive = i_src._rexiTermGammaActive;

}


bool lg::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ExpIntegration::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> lg::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lg");
	return retval;

}



bool lg::setupByKeyValue(
		const std::string &i_key,
		const std::string &i_value
)
{
	if (i_key == "ExpIntegrationFunction")
	{
		if (_expFunction.functionName != "")
			return error.set("Function name for expFunction is already set ('"+_expFunction.functionName+"')");

		if (i_value == "")
			return error.set("Empty string for expFunction given");

		_expFunction.setup(i_value);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(_expFunction);

		return true;
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


/*!
 * Setup a key which value is complex valued
 */
bool lg::setupByKeyValue(
		const std::string &i_key,
		const std::complex<double> &i_value
)
{
	if (i_key == "rexiTermAlpha")
	{
		_rexiTermAlpha = i_value;
		return true;
	}

	if (i_key == "rexiTermBeta")
	{
		_rexiTermBeta = i_value;
		return true;
	}

	if (i_key == "rexiTermGamma")
	{
		_rexiTermGamma = i_value;
		_rexiTermGammaActive = true;
		return true;
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


bool lg::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Fast gravity mode of SWE terms':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool lg::setupConfigAndForwardTimeStepperEval(
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

		if (_shackExpIntegration->direct_precompute_phin)
		{
			for (int i = 0; i < 5; i++)
				exp_coef.push_back(sweet::Data::Sphere2D::DataSpectral(_ops->sphere2DDataConfig));
		}

	}
	else if (i_evalType == EVAL_REXI_TERM)
	{
	}
	else
	{
		SWEETErrorFatal("Evaluation not supported - this error should have been caught earlier");
	}
	
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}

bool lg::setTimeStepSize(double i_dt)
{
	TimeTree_Node_LeafHelper::setTimeStepSize(i_dt);

	if (_shackExpIntegration->direct_precompute_phin)
		_computeExpDirectCoefficients();

	return true;
}


void lg::clear()
{
	TimeTree_Node_LeafHelper::clear();

	if (_shackExpIntegration->direct_precompute_phin)
		for (int i = 0; i < 5; i++)
			exp_coef[i].clear();

}


bool lg::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);

	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	double gh = _shackPDESWESphere2D->gravitation * _shackPDESWESphere2D->h0;

	// TODO: Write in a way which directly writes output to output array
	o_U.phi_pert = -gh*i_U.div;
	o_U.div = -_ops->laplace(i_U.phi_pert);
	o_U.vrt.spectral_setZero();

	return true;
}



bool lg::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	/*
	 * avg. geopotential
	 */
	double GH = _shackPDESWESphere2D->h0*_shackPDESWESphere2D->gravitation;

	sweet::Data::Sphere2D::DataSpectral rhs = i_U.div + _ops->implicit_L(i_U.phi_pert, _dt);
	o_U.div = _ops->implicit_helmholtz(rhs, -GH*_dt*_dt, _shackSphere2DDataOps->sphere2d_radius);
	o_U.phi_pert = i_U.phi_pert - _dt*GH*o_U.div;
	o_U.vrt = i_U.vrt;

	return true;
}


bool lg::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_expFunction.isSetup());


	/***********************
	 * !!! IMPORTANT !!!
	 *
	 * For evaluating phi2, phi3, etc. functions,
	 * the exponential of the steady state mode is *NOT* equal to 1
	 *
	 * Therefore, we need to have a special treatment for:
	 * 1. The D=0 case where we need to scale by the phi_n(0)
	 * 2. Copying over the vorticity field where we also need to scale by phi_n(0)
	 *
	 * !!! IMPORTANT !!!
	 */

	/*
	 * Using exponential integrators, we must compute an
	 */
	double ir = 1.0/_shackSphere2DDataOps->sphere2d_radius;

	/*
	 * See doc/time_integration/exponential_integration/swe_sphere2d_direct_exp_int_for_L_nonrotating_sphere2D/
	 */

	/*
	 * G value in Q matrix
	 */
	double G = -_shackPDESWESphere2D->gravitation*_shackPDESWESphere2D->h0;

	_expFunction.eval(std::complex<double>(0,0), expSteadyStateValue);


	bool compute_phin = false;
	// The matrices are computed if
	// - it is required by the user, or
	// - if it is the first time step, or
	// - if the time step size has changed w.r.t. the currently stored matrices
	if ( (! _shackExpIntegration->direct_precompute_phin) ||
		  ////i_timeStamp == 0                ||
		  _dt_precompute_phin != _dt       )
		compute_phin = true;

	//////if (_shackExpIntegration->direct_precompute_phin && i_timeStamp == 0)
	//////{
	//////	_dt_precompute_phin = _dt;
	//////	for (int m = 0; m <= _ops->sphere2DDataConfig->spectral_modes_m_max; m++)
	//////	{
	//////		std::size_t idx = _ops->sphere2DDataConfig->getArrayIndexByModes(m, m);
	//////		for (int n = m; n <= _ops->sphere2DDataConfig->spectral_modes_n_max; n++)
	//////		{
	//////			std::array<std::array<std::complex<double>, 2>, 2> aux = {{{0, 0}, {0, 0}}};
	//////			Z.push_back(aux);
	//////			idx++;
	//////		}
	//////	}
	//////}


	///if (compute_phin)
	///	_computeExpDirectCoefficients();


	if ( ! _shackExpIntegration->direct_precompute_phin )
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= _ops->sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = _ops->sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= _ops->sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				double D = (double)n*((double)n+1.0)*ir*ir;

				const std::complex<double> &i_phi_pert = i_U.phi_pert.spectral_space_data[idx];
				const std::complex<double> &i_div = i_U.div.spectral_space_data[idx];
				const std::complex<double> &i_vrt = i_U.vrt.spectral_space_data[idx];

				std::complex<double> &o_phi_pert = o_U.phi_pert.spectral_space_data[idx];
				std::complex<double> &o_div = o_U.div.spectral_space_data[idx];
				std::complex<double> &o_vrt = o_U.vrt.spectral_space_data[idx];

				if (D == 0)
				{
					/*
					 * the phi_n for n>=2 funcions are not evaluated to 1 at x=0
					 */
					o_phi_pert = i_phi_pert*expSteadyStateValue.real();
					o_div = i_div*expSteadyStateValue.real();
					o_vrt = i_vrt*expSteadyStateValue.real();

					idx++;
					continue;
				}

				// result will be imaginary only!
				std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

				std::complex<double> l0 = -sqrt_DG/(2*G) * i_phi_pert + 0.5*i_div;
				std::complex<double> l1 = +sqrt_DG/(2*G) * i_phi_pert + 0.5*i_div;

				std::complex<double> tmp;
				_expFunction.eval(-_dt*sqrt_DG, tmp);
				l0 *= tmp;

				_expFunction.eval(_dt*sqrt_DG, tmp);
				l1 *= tmp;

				o_phi_pert = G/sqrt_DG * (l1 - l0);
				o_div = l0 + l1;
				o_vrt = i_vrt*expSteadyStateValue.real();

				idx++;
			}
		}

	}
	else
	{

////		if (i_timeStamp == 0)
////			_computeExpDirectCoefficients();


		o_U.phi_pert = i_U.phi_pert.spectralElementwiseMultiplication(exp_coef[0]) + i_U.div.spectralElementwiseMultiplication(exp_coef[1]);
		o_U.div      = i_U.phi_pert.spectralElementwiseMultiplication(exp_coef[2]) + i_U.div.spectralElementwiseMultiplication(exp_coef[3]);
		o_U.vrt      = i_U.vrt.spectralElementwiseMultiplication(exp_coef[4]);

////		SWEET_THREADING_SPACE_PARALLEL_FOR
////		for (int m = 0; m <= _ops->sphere2DDataConfig->spectral_modes_m_max; m++)
////		{
////			std::size_t idx = _ops->sphere2DDataConfig->getArrayIndexByModes(m, m);
////			for (int n = m; n <= _ops->sphere2DDataConfig->spectral_modes_n_max; n++)
////			{
////
////				const std::complex<double> &i_phi_pert = i_U.phi_pert.spectral_space_data[idx];
////				const std::complex<double> &i_div = i_U.div.spectral_space_data[idx];
////				const std::complex<double> &i_vrt = i_U.vrt.spectral_space_data[idx];
////
////				std::complex<double> &o_phi_pert = o_U.phi_pert.spectral_space_data[idx];
////				std::complex<double> &o_div = o_U.div.spectral_space_data[idx];
////				std::complex<double> &o_vrt = o_U.vrt.spectral_space_data[idx];
////
////				double D = (double)n*((double)n+1.0)*ir*ir;
////
////				if (D == 0)
////				{
////
////					o_phi_pert = i_phi_pert*expSteadyStateValue.real();
////					o_div = i_div*expSteadyStateValue.real();
////					o_vrt = i_vrt*expSteadyStateValue.real();
////
////					idx++;
////					continue;
////				}
////
////				if (compute_phin)
////				{
////
////					///_computeExpDirectCoefficients();
////
////					std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));
////
////					std::complex<double> exp_p;
////					std::complex<double> exp_m;
////					_expFunction.eval(_dt*sqrt_DG, exp_p);
////					_expFunction.eval(-_dt*sqrt_DG, exp_m);
////					std::complex<double> exp_exp = .5 * (exp_p + exp_m);
////					std::complex<double> exp_exp_m = .5 * (exp_p - exp_m);
////
////					Z[idx][0][0] = exp_exp;
////					Z[idx][0][1] = exp_exp_m * sqrt_DG / D;
////					Z[idx][1][0] = exp_exp_m * sqrt_DG / G;
////					Z[idx][1][1] = exp_exp;
////     			}
////
////
////				o_phi_pert = Z[idx][0][0] * i_phi_pert + Z[idx][0][1] * i_div;
////				o_div      = Z[idx][1][0] * i_phi_pert + Z[idx][1][1] * i_div;
////				o_vrt = i_vrt*expSteadyStateValue.real();
////
////				idx++;
////			}
////		}
////	
	}

	return true;
}


bool lg::_computeExpDirectCoefficients()
{

	/*
	 * Using exponential integrators, we must compute an
	 */
	double ir = 1.0/_shackSphere2DDataOps->sphere2d_radius;

	/*
	 * See doc/time_integration/exponential_integration/swe_sphere2d_direct_exp_int_for_L_nonrotating_sphere2D/
	 */

	/*
	 * G value in Q matrix
	 */
	double G = -_shackPDESWESphere2D->gravitation*_shackPDESWESphere2D->h0;

	_expFunction.eval(std::complex<double>(0,0), expSteadyStateValue);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= _ops->sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		for (int n = m; n <= _ops->sphere2DDataConfig->spectral_modes_n_max; n++)
		{

			double D = (double)n*((double)n+1.0)*ir*ir;

			if (D == 0)
			{
				exp_coef[0].spectral_set(n, m, expSteadyStateValue.real());
				exp_coef[1].spectral_set(n, m, 0.);
				exp_coef[2].spectral_set(n, m, 0.);
				exp_coef[3].spectral_set(n, m, expSteadyStateValue.real());
				exp_coef[4].spectral_set(n, m, expSteadyStateValue.real());

				continue;
			}

			std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

			std::complex<double> exp_p;
			std::complex<double> exp_m;
			_expFunction.eval(_dt*sqrt_DG, exp_p);
			_expFunction.eval(-_dt*sqrt_DG, exp_m);
			std::complex<double> exp_exp = .5 * (exp_p + exp_m);
			std::complex<double> exp_exp_m = .5 * (exp_p - exp_m);

			exp_coef[0].spectral_set(n, m, exp_exp);
			exp_coef[1].spectral_set(n, m, exp_exp_m * sqrt_DG / D);
			exp_coef[2].spectral_set(n, m, exp_exp_m * sqrt_DG / G);
			exp_coef[3].spectral_set(n, m, exp_exp);
			exp_coef[4].spectral_set(n, m, expSteadyStateValue.real());

		}
	}

	return true;


}


bool lg::_eval_rexiTerm(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_Ur = cast(i_U_);
	DataContainer::Simulation &o_Ur = cast(o_U_);

	/*
	 * We can reuse the backward Euler time stepper which has the form
	 *
	 * U1 = (I - dt*L)^{-1} U0
	 *
	 * For REXI, we need to solve terms of the form
	 *
	 * U1 = \beta(dt*L - \alpha)^{-1} U0
	 *
	 * and rewrite it to
	 *
	 * U1 = -\beta / \alpha (I - (dt / \alpha)*L)^{-1} U0
	 */

	sweet::Data::Sphere2DComplex::DataSpectral o_U_phi_pert(i_Ur.phi_pert.sphere2DDataConfig);
	sweet::Data::Sphere2DComplex::DataSpectral o_U_vrt(i_Ur.phi_pert.sphere2DDataConfig);
	sweet::Data::Sphere2DComplex::DataSpectral o_U_div(i_Ur.phi_pert.sphere2DDataConfig);

	sweet::Data::Sphere2DComplex::DataSpectral i_U_phi_pert = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_Ur.phi_pert);
	sweet::Data::Sphere2DComplex::DataSpectral i_U_vrt = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_Ur.vrt);
	sweet::Data::Sphere2DComplex::DataSpectral i_U_div = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_Ur.div);

	double gh0 = _shackPDESWESphere2D->h0*_shackPDESWESphere2D->gravitation;

	std::complex<double> dtComplex = _dt/_rexiTermAlpha;

	sweet::Data::Sphere2DComplex::DataSpectral rhs = i_U_div + _opsComplex->implicit_L(i_U_phi_pert, dtComplex);
	o_U_div = _opsComplex->implicit_helmholtz(rhs, -gh0*dtComplex*dtComplex, _shackSphere2DDataOps->sphere2d_radius);
	o_U_phi_pert = i_U_phi_pert - (dtComplex*gh0)*o_U_div;
	o_U_vrt = i_U_vrt;

	std::complex<double> foo = -_rexiTermBeta/_rexiTermAlpha;
	o_U_phi_pert *= foo;
	o_U_vrt *= foo;
	o_U_div *= foo;

	if (_rexiTermGammaActive)
	{
		o_U_phi_pert += _rexiTermGamma*i_U_phi_pert;
		o_U_vrt += _rexiTermGamma*i_U_vrt;
		o_U_div += _rexiTermGamma*i_U_div;
	}

	o_Ur.phi_pert = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(o_U_phi_pert);
	o_Ur.vrt = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(o_U_vrt);
	o_Ur.div = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(o_U_div);

	return true;
}

}}}
