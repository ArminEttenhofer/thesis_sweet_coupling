#include <sweet/Data/Sphere2D/Operators.hpp>
#include "nr_vd.hpp"


#include <vector>


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {

nr_vd::nr_vd()	:
	_shackPDESWESphere2D(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

nr_vd::~nr_vd()
{
}

nr_vd::nr_vd(
		const nr_vd &i_value
)	:
	TimeTree_Node_LeafHelper(i_value)
{
	_shackPDESWESphere2D = i_value._shackPDESWESphere2D;
	_ops = i_value._ops;
}


bool nr_vd::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> nr_vd::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("nr_vd");
	return retval;

}



bool nr_vd::outputHelp(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << "InteriorNode: 'Non-linear advection term (using vorticity-divergence-based formulation)':" << std::endl;
	o_ostream << i_prefix << std::endl;
	o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;

	return true;
}


bool nr_vd::setupConfigAndForwardTimeStepperEval(
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

void nr_vd::clear()
{
	TimeTree_Node_LeafHelper::clear();
}


void nr_vd::euler_timestep_update_na(
		const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_U_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	_ops->vrtdiv_2_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	sweet::Data::Sphere2D::DataGrid U_div_phys = i_U_div.toGrid();
	o_phi_t -= _ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_phi.toGrid());
	o_vrt_t -= _ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_vrt.toGrid());
	o_div_t -= _ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_div.toGrid());
}


bool nr_vd::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWESphere2D != nullptr);


	o_U.phi_pert.spectral_setZero();
	o_U.div.spectral_setZero();
	o_U.vrt.spectral_setZero();


	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, U_u_phys, U_v_phys);

	// dt calculation starts here

	sweet::Data::Sphere2D::DataGrid U_div_phys = i_U.div.toGrid();

	o_U.phi_pert -= sweet::Data::Sphere2D::DataSpectral(i_U.phi_pert.toGrid()*i_U.div.toGrid());

	if (0)
	{
		o_U.vrt -= o_U.vrt.toGrid()*U_div_phys;
	}
	else
	{
		/*
		 * N from UV formulation
		 */
//		double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;


//		const sweet::Data::Sphere2D::DataSpectral &U_phi = i_U.phi;
		//const sweet::Data::Sphere2D::DataSpectral &U_vrt = i_U.vrt;
		const sweet::Data::Sphere2D::DataSpectral &U_div = i_U.div;


		sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
		_ops->vrtdiv_2_uv(i_U.vrt, i_U.div, U_u_phys, U_v_phys);

		sweet::Data::Sphere2D::DataGrid U_div_phys = U_div.toGrid();

		/*
		 * Velocity
		 */
		sweet::Data::Sphere2D::DataGrid vrtg = i_U.vrt.toGrid();

		sweet::Data::Sphere2D::DataGrid u_nl = U_u_phys*vrtg;
		sweet::Data::Sphere2D::DataGrid v_nl = U_v_phys*vrtg;

		sweet::Data::Sphere2D::DataSpectral vrt, div;
		_ops->uv_2_vrtdiv(u_nl, v_nl, vrt, div);
		//o_U.vrt -= div;


		/*
		 * NA part to be subtracted
		 */
		sweet::Data::Sphere2D::DataSpectral phi_tmp(i_U.phi_pert.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vrt_tmp(i_U.vrt.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div_tmp(i_U.div.sphere2DDataConfig);

		phi_tmp.spectral_setZero();
		vrt_tmp.spectral_setZero();
		div_tmp.spectral_setZero();

		euler_timestep_update_na(
				i_U.phi_pert, i_U.vrt, i_U.div,
				phi_tmp, vrt_tmp, div_tmp,
				i_timeStamp
			);

		vrt_tmp = vrt_tmp.toGrid();

		o_U.vrt += -div - vrt_tmp;
	}

	const sweet::Data::Sphere2D::DataGrid U_vrt_phys = i_U.vrt.toGrid();
	o_U.div += _ops->uv_2_vort(U_vrt_phys*U_u_phys, U_vrt_phys*U_v_phys);
	o_U.div += _ops->uv_2_div(U_div_phys*U_u_phys, U_div_phys*U_v_phys);
	o_U.div -= 0.5*_ops->laplace(U_u_phys*U_u_phys + U_v_phys*U_v_phys);
	o_U.div -= U_div_phys*U_div_phys;

	return true;
}

}}}
