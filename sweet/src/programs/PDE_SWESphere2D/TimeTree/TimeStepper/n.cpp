#include "n.hpp"

#include <sweet/Data/Sphere2D/Operators.hpp>

#include <vector>


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {

n::n()	:
	_shackPDESWESphere2D(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

n::~n()
{
}

n::n(
		const n &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWESphere2D = i_value._shackPDESWESphere2D;
	_ops = i_value._ops;
}


bool n::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> n::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("n");
	return retval;

}


bool n::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Nonlinear terms of SWE terms (na + nr)':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool n::setupConfigAndForwardTimeStepperEval(
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

	return true;
}

void n::clear()
{
	TimeTree_Node_LeafHelper::clear();
}

bool n::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);


	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	sweet::Data::Sphere2D::DataGrid ug(i_U.phi_pert.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataGrid vg(i_U.phi_pert.sphere2DDataConfig);

	sweet::Data::Sphere2D::DataGrid vrtg = i_U.vrt.toGrid();
	sweet::Data::Sphere2D::DataGrid divg = i_U.div.toGrid();
	_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, ug, vg);

	sweet::Data::Sphere2D::DataGrid phig = i_U.phi_pert.toGrid();

	sweet::Data::Sphere2D::DataGrid tmpg1 = ug*(vrtg/*+fg*/);
	sweet::Data::Sphere2D::DataGrid tmpg2 = vg*(vrtg/*+fg*/);

	_ops->uv_2_vrtdiv(tmpg1, tmpg2, o_U.div, o_U.vrt);

	o_U.vrt *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	sweet::Data::Sphere2D::DataSpectral tmpspec(i_U.phi_pert.sphere2DDataConfig);
	_ops->uv_2_vrtdiv(tmpg1,tmpg2, tmpspec, o_U.phi_pert);

	o_U.phi_pert *= -1.0;

	sweet::Data::Sphere2D::DataGrid tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = /*phig+*/tmpg;

	o_U.div += -_ops->laplace(tmpspec);

	return true;
}

}}}
