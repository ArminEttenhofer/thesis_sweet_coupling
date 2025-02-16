#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include "l.hpp"

#include <vector>

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {

l::l()	:
	_shackPDESWESphere2D(nullptr),
	_shackSphere2DDataOps(nullptr),
	_shackPDESWESphere2DBenchmarks(nullptr),
	_ops(nullptr),
	_opsComplex(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_REXI_TERM);
}

l::~l()
{
}


l::l(
		const l &i_src
)	:
	TimeTree_Node_LeafHelper(i_src),
	_evalDataREXITerm(i_src._evalDataREXITerm)
{
	_shackPDESWESphere2D = i_src._shackPDESWESphere2D;
	_shackSphere2DDataOps = i_src._shackSphere2DDataOps;
	_shackPDESWESphere2DBenchmarks = i_src._shackPDESWESphere2DBenchmarks;
	_ops = i_src._ops;
	_opsComplex = i_src._opsComplex;
}


bool l::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackPDESWESphere2DBenchmarks = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> l::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("l");
	return retval;

}


bool l::setupByKeyValue(
		const std::string &i_key,
		const std::string &i_value
)
{
	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


/*!
 * Setup a key which value is complex valued
 */
bool l::setupByKeyValue(
		const std::string &i_key,
		const std::complex<double> &i_value
)
{
	if (i_key == "rexiTermAlpha")
	{
		_evalDataREXITerm._rexiTermAlpha = i_value;
		return true;
	}

	if (i_key == "rexiTermBeta")
	{
		_evalDataREXITerm._rexiTermBeta = i_value;
		return true;
	}

	if (i_key == "rexiTermGamma")
	{
		_evalDataREXITerm._rexiTermGamma = i_value;
		_evalDataREXITerm._rexiTermGammaActive = true;
		return true;
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


bool l::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'All linear SWE terms (lg + lc)':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool l::setupConfigAndForwardTimeStepperEval(
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
		_evalDataTendencies._fg.setup(myConfig.ops->sphere2DDataConfig);
		_evalDataTendencies._ug.setup(myConfig.ops->sphere2DDataConfig);
		_evalDataTendencies._vg.setup(myConfig.ops->sphere2DDataConfig);
		_evalDataTendencies._tmpg1.setup(myConfig.ops->sphere2DDataConfig);
		_evalDataTendencies._tmpg2.setup(myConfig.ops->sphere2DDataConfig);

		if (_shackPDESWESphere2D->sphere2d_use_fsphere2D)
			_ops->getFG_fSphere2D(_shackPDESWESphere2D->sphere2d_fsphere2d_f0, _evalDataTendencies._fg);
		else
			_ops->getFG_rotatingSphere2D(_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega, _evalDataTendencies._fg);
	}
	else if (i_evalType == EVAL_EULER_BACKWARD)
	{
		if (_shackPDESWESphere2D->sphere2d_use_fsphere2D)
			SWEETErrorFatal("Not supported");

		_evalDataBackwardEuler._rhs.setup(myConfig.ops->sphere2DDataConfig);

		_evalDataBackwardEuler.sphBandedMatrix_GridReal.setup(_ops->sphere2DDataConfig, 4);
	}
	else if (i_evalType == EVAL_REXI_TERM)
	{
		if (_shackPDESWESphere2D->sphere2d_use_fsphere2D)
			SWEETErrorFatal("Not supported");

		_evalDataREXITerm._rhs.setup(myConfig.opsComplex->sphere2DDataConfig);

		_evalDataREXITerm.sphBandedMatrix_GridComplex.setup(_opsComplex->sphere2DDataConfig, 4);
	}
	else
	{
		error.set("The DE term evaluation '"+evalTypeToString(i_evalType)+"' is not supported for the DE term 'l'"+getNewLineDebugMessage());
	}

	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void l::clear()
{
	TimeTree_Node_LeafHelper::clear();
}



bool l::setTimeStepSize(double i_dt)
{
	TimeTree_Node_LeafHelper::setTimeStepSize(i_dt);

	SWEET_ASSERT(_evalTypeRequested != EVAL_NONE);

	if (_evalTypeRequested == EVAL_EULER_BACKWARD)
	{
		double two_coriolis = 2.0*_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		double dt_two_omega = _timestepSize*two_coriolis;
		double gh0 = _shackPDESWESphere2D->gravitation * _shackPDESWESphere2D->h0;

		SWEET_ASSERT(_ops->sphere2DDataConfig != nullptr);

		_evalDataBackwardEuler.sphBandedMatrix_GridReal.solver_setZero();
		_evalDataBackwardEuler.sphBandedMatrix_GridReal.solver_addComponent_implicit_J(dt_two_omega);
		_evalDataBackwardEuler.sphBandedMatrix_GridReal.solver_addComponent_implicit_FJinvF(dt_two_omega);
		_evalDataBackwardEuler.sphBandedMatrix_GridReal.solver_addComponent_implicit_L(gh0*_timestepSize, _timestepSize, _shackSphere2DDataOps->sphere2d_radius);
	}
	else if (_evalTypeRequested == EVAL_REXI_TERM)
	{
		std::complex<double> alpha = _evalDataREXITerm._rexiTermAlpha;

		double two_coriolis = 2.0*_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		std::complex<double> timestepSizeComplex = _dt/alpha;
		std::complex<double> dt_two_omega = timestepSizeComplex*two_coriolis;
		double gh0 = _shackPDESWESphere2D->gravitation * _shackPDESWESphere2D->h0;

		SWEET_ASSERT(_ops->sphere2DDataConfig != nullptr);
		
		_evalDataREXITerm.sphBandedMatrix_GridComplex.solver_setZero();
		_evalDataREXITerm.sphBandedMatrix_GridComplex.solver_addComponent_implicit_J(dt_two_omega);
		_evalDataREXITerm.sphBandedMatrix_GridComplex.solver_addComponent_implicit_FJinvF(dt_two_omega);
		_evalDataREXITerm.sphBandedMatrix_GridComplex.solver_addComponent_implicit_L(gh0*timestepSizeComplex, timestepSizeComplex, _shackSphere2DDataOps->sphere2d_radius);
	}

	return true;
}



bool l::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);


	if (!_shackPDESWESphere2D->sphere2d_use_fsphere2D)
	{
		double gh0 = _shackPDESWESphere2D->gravitation * _shackPDESWESphere2D->h0;

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
		 */
		sweet::Data::Sphere2D::DataGrid& _fg = _evalDataTendencies._fg;
		sweet::Data::Sphere2D::DataGrid& _ug = _evalDataTendencies._ug;
		sweet::Data::Sphere2D::DataGrid& _vg = _evalDataTendencies._vg;

		sweet::Data::Sphere2D::DataGrid& _tmpg1 = _evalDataTendencies._tmpg1;
		sweet::Data::Sphere2D::DataGrid& _tmpg2 = _evalDataTendencies._tmpg2;

		/*
		 * Step 1a
		 */
		_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, _ug, _vg);

		/*
		 * Step 1b
		 */
		_tmpg1 = _ug*_fg;
		_tmpg2 = _vg*_fg;

		/*
		 * Step 1c
		 */
		_ops->uv_2_vrtdiv(_tmpg1, _tmpg2, o_U.div, o_U.vrt);

		/*
		 * Step 1d
		 */
		o_U.vrt *= -1.0;

		/*
		 * Step 1e
		 */
		o_U.div += -_ops->laplace(i_U.phi_pert + _shackPDESWESphere2DBenchmarks->topography.topography*_shackPDESWESphere2D->gravitation);

		/*
		 * DIV on velocity field
		 */
		o_U.phi_pert = (-gh0)*i_U.div;
	}
	else
	{
		double gh = _shackPDESWESphere2D->gravitation * _shackPDESWESphere2D->h0;

		o_U.div = -_ops->laplace(i_U.phi_pert);

		o_U.vrt = -_shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_U.div;
		o_U.div += _shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_U.vrt;

		o_U.phi_pert = -gh*i_U.div;
	}

	return true;
}


bool l::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{

	if (_shackPDESWESphere2DBenchmarks->benchmark_name == "geostrophic_balance_topography" ||
            _shackPDESWESphere2DBenchmarks->benchmark_name == "williamson5" ||
            _shackPDESWESphere2DBenchmarks->benchmark_name == "flow_over_mountain"
            )
		SWEETErrorFatal("Not implemented.");

	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	sweet::Data::Sphere2D::DataSpectral& rhs = _evalDataBackwardEuler._rhs;

	// avg. geopotential
	double gh0 = _shackPDESWESphere2D->h0*_shackPDESWESphere2D->gravitation;
	double dt_two_omega = _dt*2.0*_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

	rhs = i_U.div + _ops->implicit_FJinv(i_U.vrt, dt_two_omega) + _ops->implicit_L(i_U.phi_pert, _dt);

	o_U.div = _evalDataBackwardEuler.sphBandedMatrix_GridReal.solve(rhs);
	o_U.phi_pert = i_U.phi_pert - _dt*gh0*o_U.div;
	o_U.vrt = _ops->implicit_Jinv(i_U.vrt - _ops->implicit_F(o_U.div, dt_two_omega), dt_two_omega);

	return true;
}


/*!
 *
 * We can reuse the backward Euler time stepper which has the form
 * \f[
 * 	U_1 = (I - \Delta t \cdot L)^{-1} U_0
 * \f]
 *
 * For REXI, we need to solve terms of the form
 *
 * \f[
 * 	U_1 = \beta(\Delta t \cdot L - \alpha)^{-1} U_0
 * \f]
 *
 * and rewrite it to
 *
 * \f[
 * 	U_1 = -\frac{\beta}{\alpha} (I - \frac{\Delta t}{\alpha}*L)^{-1} U_0
 * \f]
 *
 */
bool l::_eval_rexiTerm(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	if (_shackPDESWESphere2DBenchmarks->benchmark_name == "geostrophic_balance_topography" ||
            _shackPDESWESphere2DBenchmarks->benchmark_name == "williamson5" ||
            _shackPDESWESphere2DBenchmarks->benchmark_name == "flow_over_mountain"
            )
		SWEETErrorFatal("Not implemented.");

	const DataContainer::Simulation &i_Ur = cast(i_U_);
	DataContainer::Simulation &o_Ur = cast(o_U_);


	std::complex<double> alpha = _evalDataREXITerm._rexiTermAlpha;
	std::complex<double> beta = _evalDataREXITerm._rexiTermBeta;
	std::complex<double> gamma = _evalDataREXITerm._rexiTermGamma;

	sweet::Data::Sphere2DComplex::DataSpectral o_U_phi_pert(i_Ur.phi_pert.sphere2DDataConfig);
	sweet::Data::Sphere2DComplex::DataSpectral o_U_vrt(i_Ur.phi_pert.sphere2DDataConfig);
	sweet::Data::Sphere2DComplex::DataSpectral o_U_div(i_Ur.phi_pert.sphere2DDataConfig);

	const sweet::Data::Sphere2DComplex::DataSpectral i_U_phi_pert = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_Ur.phi_pert);
	const sweet::Data::Sphere2DComplex::DataSpectral i_U_vrt = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_Ur.vrt);
	const sweet::Data::Sphere2DComplex::DataSpectral i_U_div = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_Ur.div);


	sweet::Data::Sphere2DComplex::DataSpectral& rhs = _evalDataREXITerm._rhs;

	std::complex<double> dtComplex = _dt/alpha;

	//avg. geopotential
	double gh0 = _shackPDESWESphere2D->h0*_shackPDESWESphere2D->gravitation;
	std::complex<double> dt_two_omega = dtComplex*2.0*_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

	rhs = i_U_div + _opsComplex->implicit_FJinv(i_U_vrt, dt_two_omega) + _opsComplex->implicit_L(i_U_phi_pert, dtComplex);

	o_U_div = _evalDataREXITerm.sphBandedMatrix_GridComplex.solve(rhs);
	o_U_phi_pert = i_U_phi_pert - dtComplex*gh0*o_U_div;
	o_U_vrt = _opsComplex->implicit_Jinv(i_U_vrt - _opsComplex->implicit_F(o_U_div, dt_two_omega), dt_two_omega);

	std::complex<double> foo = -beta/alpha;
	o_U_phi_pert *= foo;
	o_U_vrt *= foo;
	o_U_div *= foo;

	if (_evalDataREXITerm._rexiTermGammaActive)
	{
		o_U_phi_pert += gamma*i_U_phi_pert;
		o_U_vrt += gamma*i_U_vrt;
		o_U_div += gamma*i_U_div;
	}

	o_Ur.phi_pert = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(o_U_phi_pert);
	o_Ur.vrt = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(o_U_vrt);
	o_Ur.div = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(o_U_div);

	return true;
}

}}}
