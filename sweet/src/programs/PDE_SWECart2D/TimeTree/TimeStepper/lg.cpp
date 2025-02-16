#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Cart2DComplex_DataSpectral.hpp>
//#include <sweet/Data/Cart2D/Convert/Complex_DataSpectral_2_DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataSpectral_2_Cart2D_DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include "lg.hpp"

#include <complex>
#include <cmath>
#include <vector>

namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {


lg::lg()	:
	_shackPDESWECart2D(nullptr),
	_shackCart2DDataOps(nullptr),
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
	_shackPDESWECart2D = i_src._shackPDESWECart2D;
	_shackCart2DDataOps = i_src._shackCart2DDataOps;
	_shackExpIntegration = i_src._shackExpIntegration;
	_ops = i_src._ops;
	_opsComplex = i_src._opsComplex;

	pdeSWECart2DNormalModes = i_src.pdeSWECart2DNormalModes;

	_rexiTermAlpha = i_src._rexiTermAlpha;
	_rexiTermBeta = i_src._rexiTermBeta;
	_rexiTermGamma = i_src._rexiTermGamma;
	_rexiTermGammaActive = i_src._rexiTermGammaActive;
}


bool lg::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWECart2D = io_shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ExpIntegration::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	pdeSWECart2DNormalModes.shackRegistration(*io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWECart2DNormalModes);

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

	if (i_key == "rexiTermPreallocation")
	{
		return true;
#if 0
		if (i_value == "true")
		{
			_rexiTermPreallocation = true;
			return true;
		}

		if (i_value == "false")
		{
			_rexiTermPreallocation = false;
			return true;
		}
#endif
		return error.set("Invalid key-value combination of '"+i_key+"' => '"+i_value+"'");
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


bool lg::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun* o_timeStepper
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

	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	if (i_evalType == EVAL_EXPONENTIAL)
	{
		if (_expFunction.functionName == "")
		{
			// set default to phi0
			_expFunction.setup("phi0");
		}
	}

	return true;
}


void lg::clear()
{
	TimeTree_Node_LeafHelper::clear();
}


/*
 * Return the time tendencies of the PDE term
 */
bool lg::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWECart2D != nullptr);

	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	// A- grid method
	if (!_shackCart2DDataOps->space_grid_use_c_staggering)
	{
		/*
		 * linearized non-conservative (advective) formulation:
		 *
		 * h_t = -h0*u_x - h0*v_ym
		 * u_t = -g * h_x + f*v
		 * v_t = -g * h_y - f*u
		 */

#if 1
		o_U.u = -_shackPDESWECart2D->gravitation*_ops->diff_c_x(i_U.h_pert);
		o_U.v = -_shackPDESWECart2D->gravitation*_ops->diff_c_y(i_U.h_pert);

		// standard update
		o_U.h_pert = -(_ops->diff_c_x(i_U.u) + _ops->diff_c_y(i_U.v))*_shackPDESWECart2D->h0;
#else

	#if 0
		// U-only
		o_u_t = -shackPDESWECart2D->gravitation*ops->diff_c_x(i_h) + shackPDESWECart2D->cart2d_rotating_f0*i_v;
		//o_v_t.grid_set_zero();
		o_v_t = - shackPDESWECart2D->cart2d_rotating_f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_x(i_u))*shackPDESWECart2D->h0;

	#else
		// V-only
		//o_u_t.spectral_set_zero();
		o_u_t = +shackPDESWECart2D->cart2d_rotating_f0*i_v;
		o_v_t = -shackPDESWECart2D->gravitation*ops->diff_c_y(i_h) - shackPDESWECart2D->cart2d_rotating_f0*i_u;// - shackDict.sim.f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_y(i_v))*shackPDESWECart2D->h0;
	#endif

#endif
	}
	else // shackDict.disc.use_staggering = true
	{
		// STAGGERED GRID

		/*
		 * Sadourny energy conserving scheme
		 *
		 * Note, that this grid does not follow the formulation
		 * in the paper of Robert Sadourny, but looks as follows:
		 *
		 *              ^
		 *              |
		 *       ______v0,1_____
		 *       |             |
		 *       |			   |
		 *       |             |
		 *  u0,0 |->  H/P0,0   |u1,0 ->
		 *(0,0.5)|			   |
		 *       |      ^      |
		 *   q0,0|______|______|
		 * (0,0)      v0,0
		 *           (0.5,0)
		 *
		 * V_t + q N x (P V) + grad( g P + 1/2 V*V) = 0
		 * P_t + div(P V) = 0
		 */

		sweet::Data::Cart2D::DataSpectral H = _shackPDESWECart2D->gravitation*i_U.h_pert;// + 0.5*(ops->avg_f_x(i_u*i_u) + ops->avg_f_y(i_v*i_v));

		sweet::Data::Cart2D::DataGrid o_u_t_phys(o_U.u.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid o_v_t_phys(o_U.v.cart2DDataConfig);
		o_u_t_phys =  - _ops->diff_b_x(H).toGrid();
		o_v_t_phys =  - _ops->diff_b_y(H).toGrid();
		o_U.u.loadCart2DDataGrid(o_u_t_phys);
		o_U.v.loadCart2DDataGrid(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		o_U.h_pert = -_ops->diff_f_x(_shackPDESWECart2D->h0*i_U.u) - _ops->diff_f_y(_shackPDESWECart2D->h0*i_U.v);
	}

	return true;
}


/*
 * Backward Euler evaluation of the term
 *
 * U1' = L U1 = (U1-U0)/dt
 * <=> L U1 dt = U1 - U0
 * <=> U1 - L U1 dt = U0
 * <=> (I - dt*L) U1 = U0
 *
 * U1 = (I - dt*L)^{-1} U0
 */
bool lg::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	// Helmholtz eq. with f = 0

#if SWEET_USE_CART2D_SPECTRAL_SPACE

	sweet::Data::Cart2D::DataSpectral &eta0 = i_U.h_pert;
	sweet::Data::Cart2D::DataSpectral &u0 = i_U.u;
	sweet::Data::Cart2D::DataSpectral &v0 = i_U.v;

	double alpha = 1.0/_dt;
	double f = 0.;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + f*f;

	double eta_bar = _shackPDESWECart2D->h0;
	double g = _shackPDESWECart2D->gravitation;

	sweet::Data::Cart2D::DataSpectral rhs =
			(kappa/alpha) * eta0
			- eta_bar*(_ops->diff_c_x(u0) + _ops->diff_c_y(v0))
			- (f*eta_bar/alpha) * (_ops->diff_c_x(v0) - _ops->diff_c_y(u0))
		;

	sweet::Data::Cart2D::DataSpectral lhs = (-g*eta_bar*(_ops->diff2_c_x + _ops->diff2_c_y)).spectral_addScalarAll(kappa);
	o_U.h_pert = rhs.spectral_div_element_wise(lhs);

	sweet::Data::Cart2D::DataSpectral uh = u0 - g*_ops->diff_c_x(o_U.h_pert);
	sweet::Data::Cart2D::DataSpectral vh = v0 - g*_ops->diff_c_y(o_U.h_pert);

	o_U.u = alpha/kappa * uh     + f/kappa * vh;
	o_U.v =    -f/kappa * uh + alpha/kappa * vh;

#else

	sweet::Data::Cart2D::DataSpectralComplex eta0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(i_U.h_pert);
	sweet::Data::Cart2D::DataSpectralComplex u0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(i_U.u);
	sweet::Data::Cart2D::DataSpectralComplex v0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(i_U.v);

	double alpha = 1.0/i_dt;
	double f = 0.;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + f*f;

	double eta_bar = _shackPDESWECart2D->h0;
	double g = _shackPDESWECart2D->gravitation;

	sweet::Data::Cart2D::DataSpectralComplex rhs =
			(kappa/alpha) * eta0
			- eta_bar*(opComplex.diff_c_x(u0) + opComplex.diff_c_y(v0))
			- (f/alpha) * (opComplex.diff_c_x(v0) - opComplex.diff_c_y(u0))
		;

	sweet::Data::Cart2D::DataSpectralComplex lhs = (-g*eta_bar*(opComplex.diff2_c_x + opComplex.diff2_c_y)).spectral_addScalarAll(kappa);
	sweet::Data::Cart2D::DataSpectralComplex eta = rhs.spectral_div_element_wise(lhs);

	sweet::Data::Cart2D::DataSpectralComplex uh = u0 - g*opComplex.diff_c_x(eta);
	sweet::Data::Cart2D::DataSpectralComplex vh = v0 - g*opComplex.diff_c_y(eta);

	sweet::Data::Cart2D::DataSpectralComplex u1 = alpha/kappa * uh     + f/kappa * vh;
	sweet::Data::Cart2D::DataSpectralComplex v1 =    -f/kappa * uh + alpha/kappa * vh;

	o_U.h_pert = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(eta);
	o_U.u = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(u1);
	o_U.v = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(v1);
#endif

	return true;
}



/*!
 * Evaluate exponential of linear term
 */
bool lg::_eval_exponential(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_expFunction.isSetup());

	#define T double

	typedef std::complex<T> complex;
	complex I(0.0, 1.0);

	T dt = _dt;

	bool compute_phin = false;
	if ( (! _shackExpIntegration->direct_precompute_phin) ||
		  i_timeStamp == 0                ||
		  _dt_precompute_phin != dt       )
		compute_phin = true;

	if (compute_phin)
	{
		if (_shackExpIntegration->direct_precompute_phin)
		{
			for (std::size_t ik1 = 0; ik1 < _ops->cart2DDataConfig->spectral_data_size[1]; ik1++)
			{
				//if (i_simulation_timestamp == 0)
				//{
					_dt_precompute_phin = _dt;
					std::vector<std::array<std::array<complex, 3>, 3>> aux = {};
					Z.push_back(aux);  // Z[k1][k2][0,1,2][0,1,2];
				//}
				for (std::size_t ik0 = 0; ik0 < _ops->cart2DDataConfig->spectral_data_size[0]; ik0++)
				{
					//if (i_simulation_timestamp == 0)
					//{
						std::array<std::array<complex, 3>, 3> aux2 = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
						Z[ik1].push_back(aux2);
								//}
				}
			}
		}
	}

#if SWEET_THREADING_SPACE
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
#endif
	for (std::size_t ik1 = 0; ik1 < _ops->cart2DDataConfig->spectral_data_size[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < _ops->cart2DDataConfig->spectral_data_size[0]; ik0++)
		{
			// Light alternative
			std::array<std::array<std::complex<T>, 3>, 3> Z_single_wavenumber_pair;  // to avoid too much memory allocation in very large simulations;

			T k1;
			if (ik1 < _ops->cart2DDataConfig->spectral_data_size[1]/2)
				k1 = (T)ik1;
			else
				k1 = (T)((int)ik1-(int)_ops->cart2DDataConfig->spectral_data_size[1]);

			T k0 = (T)ik0;

			complex U[3];
			U[0] = i_U.h_pert.spectral_get(ik1, ik0);
			U[1] = i_U.u.spectral_get(ik1, ik0);
			U[2] = i_U.v.spectral_get(ik1, ik0);

			if (compute_phin)
			{
				/*
				 * Matrix with Eigenvectors (column-wise) and its invers
				 */
				complex v[3][3];
				complex v_inv[3][3];

				/*
				 * Eigenvalues
				 */
				complex lambda[3];

				// Eigendecomposition with f = 0
				pdeSWECart2DNormalModes.sw_eigen_decomp(
									k0, k1, false, v, lambda,
									0,
									_shackPDESWECart2D->h0,
									_shackPDESWECart2D->gravitation
								);
				pdeSWECart2DNormalModes.sw_eigen_decomp(
									k0, k1, true, v_inv, lambda,
									0,
									_shackPDESWECart2D->h0,
									_shackPDESWECart2D->gravitation
								);

				complex v_lambda[3][3];

				for (int i = 0; i < 3; i++)
				{
					std::complex<T> &lam = lambda[i];

					std::complex<T> K;
					_expFunction.eval(lam*dt, K);
					for (int j = 0; j < 3; j++)
						v_lambda[j][i] = v[j][i] * K;
				}

				for (int j = 0; j < 3; j++)
					for (int i = 0; i < 3; i++)
					{
						std::complex<double> d = 0.;
						for (int k = 0; k < 3; k++)
							d += v_lambda[j][k] * v_inv[k][i];

						if (_shackExpIntegration->direct_precompute_phin)
							Z[ik1][ik0][j][i] = d;
						else // light version
							Z_single_wavenumber_pair[j][i] = d;
					}
			}

			complex U_copy[3];
			for (int k = 0; k < 3; k++)
			{
				U_copy[k] = U[k];
				U[k] = 0.0;
			}

			if (_shackExpIntegration->direct_precompute_phin)
			{
				for (int k = 0; k < 3; k++)
					for (int j = 0; j < 3; j++)
						U[k] += Z[ik1][ik0][k][j] * U_copy[j];
			}
			else // light version
			{
				for (int k = 0; k < 3; k++)
					for (int j = 0; j < 3; j++)
						U[k] += Z_single_wavenumber_pair[k][j] * U_copy[j];
			}


#if SWEET_QUADMATH
			std::complex<double> tmp0(U[0].real(), U[0].imag());
			o_U.h_pert.spectral_set(ik1, ik0, tmp0);

			std::complex<double> tmp1(U[1].real(), U[1].imag());
			o_U.u.spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U[2].real(), U[2].imag());
			o_U.v.spectral_set(ik1, ik0, tmp2);
#else
			o_U.h_pert.spectral_set(ik1, ik0, U[0]);
			o_U.u.spectral_set(ik1, ik0, U[1]);
			o_U.v.spectral_set(ik1, ik0, U[2]);
#endif
		}
	}
	
	o_U.h_pert.spectral_zeroAliasingModes();
	o_U.u.spectral_zeroAliasingModes();
	o_U.v.spectral_zeroAliasingModes();

	return true;

}



/*!
 * Evaluate REXI term
 *
 * o_U = \beta (dt*L - \alpha)^-1 i_U
 */
bool lg::_eval_rexiTerm(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_Ur = cast(i_U_);
	DataContainer::Simulation &o_Ur = cast(o_U_);

	// TODO
	SWEETErrorFatal("Not implemented (yet!)");

	return true;

//////	/*
//////	 * We can reuse the backward Euler time stepper which has the form
//////	 *
//////	 * U1 = (I - dt*L)^{-1} U0
//////	 *
//////	 * For REXI, we need to solve terms of the form
//////	 *
//////	 * U1 = \beta(dt*L - \alpha)^{-1} U0
//////	 *
//////	 * and rewrite it to
//////	 *
//////	 * U1 = -\beta / \alpha (I - (dt / \alpha)*L)^{-1} U0
//////	 */
//////
//////	sweet::Data::Sphere2DComplex::DataSpectral o_U_phi_pert(i_Ur.phi_pert.sphere2DDataConfig);
//////	sweet::Data::Sphere2DComplex::DataSpectral o_U_vrt(i_Ur.phi_pert.sphere2DDataConfig);
//////	sweet::Data::Sphere2DComplex::DataSpectral o_U_div(i_Ur.phi_pert.sphere2DDataConfig);
//////
//////	sweet::Data::Sphere2DComplex::DataSpectral i_U_phi_pert = sweet::Data::Sphere2D::Convert::DataSpectral_2_Complex_DataSpectral::grid_convert(i_Ur.phi_pert);
//////	sweet::Data::Sphere2DComplex::DataSpectral i_U_vrt = sweet::Data::Sphere2D::Convert::DataSpectral_2_Complex_DataSpectral::grid_convert(i_Ur.vrt);
//////	sweet::Data::Sphere2DComplex::DataSpectral i_U_div = sweet::Data::Sphere2D::Convert::DataSpectral_2_Complex_DataSpectral::grid_convert(i_Ur.div);
//////
//////	double GH = _shackPDESWESphere2D->h0*_shackPDESWESphere2D->gravitation;
//////
//////	std::complex<double> dtComplex = _dt/_rexiTermAlpha;
//////
//////	sweet::Data::Sphere2DComplex::DataSpectral rhs = i_U_div + _opsComplex->implicit_L(i_U_phi_pert, dtComplex);
//////	o_U_div = _opsComplex->implicit_helmholtz(rhs, -GH*dtComplex*dtComplex, _shackSphere2DDataOps->sphere2d_radius);
//////	o_U_phi_pert = i_U_phi_pert - (dtComplex*GH)*o_U_div;
//////	o_U_vrt = i_U_vrt;
//////
//////	o_U_phi_pert *= -_rexiTermBeta/_rexiTermAlpha;
//////	o_U_vrt *= -_rexiTermBeta/_rexiTermAlpha;
//////	o_U_div *= -_rexiTermBeta/_rexiTermAlpha;
//////
//////	if (_rexiTermGammaActive)
//////	{
//////		o_U_phi_pert += _rexiTermGamma*i_U_phi_pert;
//////		o_U_vrt += _rexiTermGamma*i_U_vrt;
//////		o_U_div += _rexiTermGamma*i_U_div;
//////	}
//////
//////	o_Ur.phi_pert = sweet::Data::Sphere2D::Convert::Complex_DataSpectral_2_DataSpectral::grid_convert_real(o_U_phi_pert);
//////	o_Ur.vrt = sweet::Data::Sphere2D::Convert::Complex_DataSpectral_2_DataSpectral::grid_convert_real(o_U_vrt);
//////	o_Ur.div = sweet::Data::Sphere2D::Convert::Complex_DataSpectral_2_DataSpectral::grid_convert_real(o_U_div);
//////
//////	return true;
}

}}}
