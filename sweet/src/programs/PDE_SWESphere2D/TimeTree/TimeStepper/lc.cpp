#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Error/Assert.hpp>
#include "lc.hpp"


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


lc::lc()	:
	_shackPDESWESphere2D(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

lc::~lc()
{
}


lc::lc(
		const lc &i_value
)	:
	sweet::TimeTree::TimeTree_Node_LeafHelper<lc>(i_value)
{
	_shackPDESWESphere2D = i_value._shackPDESWESphere2D;
	_ops = i_value._ops;
}

bool lc::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);
	return true;
}

const std::vector<std::string> lc::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lc");
	return retval;

}

bool lc::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Coriolis term of SWE terms':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}


void lc::_setupDataBuffers()
{
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);
	SWEET_ASSERT(_ops != nullptr);

	if (_shackPDESWESphere2D->sphere2d_use_fsphere2D)
		fg = _ops->getFG_fSphere2D(_shackPDESWESphere2D->sphere2d_fsphere2d_f0);
	else
		fg = _ops->getFG_rotatingSphere2D(_shackPDESWESphere2D->sphere2d_rotating_coriolis_omega);

	ug.setup(_ops->sphere2DDataConfig);
	vg.setup(_ops->sphere2DDataConfig);
}


bool lc::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun *o_timeStepper
)
{
	const DataContainer::Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	_setupDataBuffers();

	// default setup
	TimeTree_Node_Base::registerTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void lc::clear()
{
	TimeTree_Node_LeafHelper::clear();
}


bool lc::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	/*
	 * step 1a
	 */
	_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, ug, vg);

	/*
	 * step 1c
	 */

	/*
	 * step 1b
	 */
	sweet::Data::Sphere2D::DataGrid tmp_u = ug*fg;
	sweet::Data::Sphere2D::DataGrid tmp_v = vg*fg;

	_ops->uv_2_vrtdiv(tmp_u, tmp_v, o_U.div, o_U.vrt);

	/*
	 * step 1d
	 */
	o_U.vrt *= -1.0;


	/*
	 * step 1e
	 * Nothing to do
	 */

	/*
	 * step 2a
	 * Zero tendencies
	 */
	o_U.phi_pert.spectral_setZero();
	return true;
}

}}}
