#include <sweet/Data/Sphere2D/Operators.hpp>
#include "nr_uv.hpp"
#include <vector>


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


nr_uv::nr_uv()	:
	_shackPDESWESphere2D(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}


nr_uv::~nr_uv()
{
}

nr_uv::nr_uv(
		const nr_uv &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWESphere2D = i_value._shackPDESWESphere2D;
	_ops = i_value._ops;
}

bool nr_uv::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}


const std::vector<std::string> nr_uv::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("nr");
	retval.push_back("nr_uv");
	return retval;

}



bool nr_uv::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Non-linear remainder term (using velocity-based formulation)':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool nr_uv::setupConfigAndForwardTimeStepperEval(
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

void nr_uv::clear()
{
	TimeTree_Node_LeafHelper::clear();
}

bool nr_uv::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);

	o_U.phi_pert = -sweet::Data::Sphere2D::DataSpectral(i_U.phi_pert.toGrid()*i_U.div.toGrid());
	o_U.vrt.spectral_setZero();
	o_U.div.spectral_setZero();

	return true;
}

}}}
