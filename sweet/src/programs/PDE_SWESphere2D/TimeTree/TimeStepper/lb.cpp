#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include "lb.hpp"

#include <complex>
#include <cmath>
#include <vector>

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


lb::lb()	:
	_shackPDESWESphere2D(nullptr),
	_shackSphere2DDataOps(nullptr),
	_shackPDESWESphere2DBenchmarks(nullptr),
	_ops(nullptr),
	_opsComplex(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EULER_BACKWARD);
}


lb::~lb()
{
}


lb::lb(
		const lb &i_src
)	:
	TimeTree_Node_LeafHelper(i_src)
{
	_shackSphere2DDataOps = i_src._shackSphere2DDataOps;
	_shackPDESWESphere2D = i_src._shackPDESWESphere2D;
	_shackPDESWESphere2DBenchmarks = i_src._shackPDESWESphere2DBenchmarks;
	_ops = i_src._ops;
	_opsComplex = i_src._opsComplex;
}


bool lb::shackRegistration(
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

const std::vector<std::string> lb::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lb");
	return retval;

}



bool lb::setupByKeyValue(
		const std::string &i_key,
		const std::string &i_value
)
{
	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


/*!
 * Setup a key which value is complex valued
 */
bool lb::setupByKeyValue(
		const std::string &i_key,
		const std::complex<double> &i_value
)
{
	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


bool lb::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Fast gravity mode of SWE terms with bathymetry':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}

bool lb::setupConfigAndForwardTimeStepperEval(
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
	else
	{
		SWEETErrorFatal("Evaluation not supported - this error should have been caught earlier");
	}
	
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void lb::clear()
{
	TimeTree_Node_LeafHelper::clear();
}


bool lb::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);

	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	// TODO: Write in a way which directly writes output to output array
	o_U.phi_pert.spectral_setZero();
	o_U.div = -_ops->laplace(_shackPDESWESphere2DBenchmarks->topography.topography*_shackPDESWESphere2D->gravitation);
	o_U.vrt.spectral_setZero();

	return true;
}



bool lb::_eval_eulerBackward(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	_eval_tendencies(i_U_, o_U_, i_timeStamp);


	//SWEETErrorFatal("Not implemented!");
	//////*
	///// * avg. geopotential
	///// */
	/////double GH = _shackPDESWESphere2D->h0*_shackPDESWESphere2D->gravitation;

	/////sweet::Data::Sphere2D::DataSpectral rhs = i_U.div + _ops->implicit_L(i_U.phi_pert, _dt);
	/////o_U.div = _ops->implicit_helmholtz(rhs, -GH*_dt*_dt, _shackSphere2DDataOps->sphere2d_radius);
	/////o_U.phi_pert = i_U.phi_pert - _dt*GH*o_U.div;
	/////o_U.vrt = i_U.vrt;

	return true;
}


}}}
