#include <sweet/Data/Sphere2D/Operators.hpp>
#include "ln.hpp"


#include <vector>


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


ln::ln()	:
	_shackPDESWESphere2D(nullptr),
	_shackPDESWESphere2DBenchmarks(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

ln::~ln()
{
}


ln::ln(
		const ln &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWESphere2D = i_value._shackPDESWESphere2D;
	_shackPDESWESphere2DBenchmarks = i_value._shackPDESWESphere2DBenchmarks;
	_ops = i_value._ops;
}


bool ln::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackPDESWESphere2DBenchmarks = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> ln::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("ln");
	return retval;

}


bool ln::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'All linear and nonlinear terms of SWE':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool ln::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun *o_timeStepper
)
{
	const DataContainer::Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	if (_shackPDESWESphere2D->sphere2d_use_fsphere2D)
		_fg = _ops->getFG_fSphere2D(_shackPDESWESphere2D->sphere2d_fsphere2d_f0);
	else
		_fg = _ops->getFG_rotatingSphere2D(_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega);

	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void ln::clear()
{
	TimeTree_Node_LeafHelper::clear();
}

bool ln::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);


	double gh0 = _shackPDESWESphere2D->gravitation * _shackPDESWESphere2D->h0;

	const sweet::Data::Sphere2D::DataGrid topography = _shackPDESWESphere2DBenchmarks->topography.topography_grid;

	/*
	 * NON-LINEAR SWE
	 *
	 * See
	 * 	Williamson, David L., Drake, John B., Hack, James J., Jakob, Rudiger, & Swarztrauber, Paul N. (1992).
	 * 	A standard test set for numerical approximations to the shallow water equations in spherical geometry.
	 * 	Journal of Computational Physics, 102(1), 211–224. https://doi.org/10.1016/S0021-9991(05)80016-6
	 *
	 * "2.3 Vorticity/Divergence Form"
	 */

	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
	 */
	sweet::Data::Sphere2D::DataGrid phi_pert_phys = i_U.phi_pert.toGrid();

	/*
	 * Step 1a
	 */
	sweet::Data::Sphere2D::DataGrid ug, vg;
	_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, ug, vg);

	/*
	 * Step 1b
	 */
	sweet::Data::Sphere2D::DataGrid vrtg = i_U.vrt.toGrid();

	/*
	 * Step 1c
	 */

	using namespace sweet;

	// left part of eq. (19)
	sweet::Data::Sphere2D::DataGrid u_nl = ug*(vrtg+_fg);

	// left part of eq. (20)
	sweet::Data::Sphere2D::DataGrid v_nl = vg*(vrtg+_fg);

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	_ops->uv_2_vrtdiv(u_nl, v_nl, o_U.div, o_U.vrt);


	/*
	 * Step 1e
	 */
	o_U.vrt *= -1.0;

	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	sweet::Data::Sphere2D::DataGrid tmpg = 0.5*(ug*ug+vg*vg);

	////sweet::Data::Sphere2D::DataSpectral e = phi_pert_phys+tmpg;
	sweet::Data::Sphere2D::DataSpectral e = phi_pert_phys+tmpg + topography*_shackPDESWESphere2D->gravitation ;

	/*
	 * Step 1g
	 */
	o_U.div -= _ops->laplace(e);

	/*
	 * Compute Phi geopotential tendencies
	 */

	/*
	 * Step 2a
	 */
	u_nl = ug*(phi_pert_phys + gh0);
	v_nl = vg*(phi_pert_phys + gh0);

	_ops->uv_2_vrtdiv(u_nl,v_nl, e, o_U.phi_pert);

	o_U.phi_pert *= -1.0;

	return true;
}

}}}
