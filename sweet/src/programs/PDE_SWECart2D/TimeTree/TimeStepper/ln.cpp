#include <sweet/Data/Cart2D/Operators.hpp>
#include "ln.hpp"


#include <vector>


namespace PDE_SWECart2D {
namespace TimeTree {
namespace TimeStepper {


ln::ln()	:
	_shackPDESWECart2D(nullptr),
	_shackCart2DDataOps(nullptr),
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
	_shackPDESWECart2D = i_value._shackPDESWECart2D;
	_shackCart2DDataOps = i_value._shackCart2DDataOps;
	_ops = i_value._ops;
}


bool ln::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	_shackPDESWECart2D = io_shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> ln::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("ln");
	return retval;

}


bool ln::setupConfigAndForwardTimeStepperEval(
	const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
	TIME_STEPPER_TYPES i_evalType,
	TimeTree_Node_Base::EvalFun* o_timeStepper
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


void ln::clear()
{
	TimeTree_Node_LeafHelper::clear();
}

/*
 * Return the time tendencies of the PDE term
 */
bool ln::_eval_tendencies(
		const sweet::Data::GenericContainer::Base &i_U_,
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
)
{
	const DataContainer::Simulation &i_U = cast(i_U_);
	DataContainer::Simulation &o_U = cast(o_U_);

	SWEET_ASSERT(_ops != nullptr);
	SWEET_ASSERT(_shackPDESWECart2D != nullptr);

	// A-grid method
	if (!_shackCart2DDataOps->space_grid_use_c_staggering)
	{
		/*
		 * non-conservative (advective) formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */

		sweet::Data::Cart2D::DataSpectral total_h = i_U.h_pert + _shackPDESWECart2D->h0;

		o_U.u = -_shackPDESWECart2D->gravitation*_ops->diff_c_x(total_h) - i_U.u*_ops->diff_c_x(i_U.u) - i_U.v*_ops->diff_c_y(i_U.u);
		o_U.v = -_shackPDESWECart2D->gravitation*_ops->diff_c_y(total_h) - i_U.u*_ops->diff_c_x(i_U.v) - i_U.v*_ops->diff_c_y(i_U.v);

		o_U.u += _shackPDESWECart2D->cart2d_rotating_f0*i_U.v;
		o_U.v -= _shackPDESWECart2D->cart2d_rotating_f0*i_U.u;

		// standard update
		/*
		 * P UPDATE
		 */
		if (!_shackPDESWECart2D->use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			//o_h_t = -ops->diff_f_x(U) - ops->diff_f_y(V);
			o_U.h_pert = -_ops->diff_c_x(i_U.u*total_h) - _ops->diff_c_y(i_U.v*total_h);
		}
		else // use linear divergence
		{
			//o_h_t = -ops->diff_f_x(shackPDESWECart2D->h0*i_u) - ops->diff_f_y(shackPDESWECart2D->h0*i_v);
			o_U.h_pert = -i_U.u*_ops->diff_c_x(total_h) - i_U.v*_ops->diff_c_y(total_h) + //nonlinear adv
					-_ops->diff_c_x(i_U.u*_shackPDESWECart2D->h0) - _ops->diff_c_y(i_U.v*_shackPDESWECart2D->h0); //linear div
		}

	}
	else // shackDict.disc.use_staggering = true
	{
		// STAGGERED GRID

		sweet::Data::Cart2D::DataSpectral U(i_U.h_pert.cart2DDataConfig); // U flux
		sweet::Data::Cart2D::DataSpectral V(i_U.h_pert.cart2DDataConfig); // V flux
		sweet::Data::Cart2D::DataSpectral H(i_U.h_pert.cart2DDataConfig); //Bernoulli potential

		sweet::Data::Cart2D::DataGrid U_phys(i_U.h_pert.cart2DDataConfig); // U flux
		sweet::Data::Cart2D::DataGrid V_phys(i_U.h_pert.cart2DDataConfig); // V flux
		sweet::Data::Cart2D::DataGrid H_phys(i_U.h_pert.cart2DDataConfig); //Bernoulli potential

		sweet::Data::Cart2D::DataGrid i_u_phys = i_U.u.toGrid();
		sweet::Data::Cart2D::DataGrid i_v_phys = i_U.v.toGrid();

		sweet::Data::Cart2D::DataGrid total_h_phys = i_U.h_pert.toGrid() + _shackPDESWECart2D->h0;
		sweet::Data::Cart2D::DataSpectral total_h(i_U.h_pert.cart2DDataConfig);
		total_h.loadCart2DDataGrid(total_h_phys);


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
		/*
		 * U and V updates
		 */

		U_phys = _ops->avg_b_x(total_h_phys)*i_u_phys;
		V_phys = _ops->avg_b_y(total_h_phys)*i_v_phys;
		H_phys = _shackPDESWECart2D->gravitation*total_h_phys + 0.5*(_ops->avg_f_x(i_u_phys*i_u_phys) + _ops->avg_f_y(i_v_phys*i_v_phys));

		U.loadCart2DDataGrid(U_phys);
		V.loadCart2DDataGrid(V_phys);
		H.loadCart2DDataGrid(H_phys);


		// Potential vorticity
		sweet::Data::Cart2D::DataGrid total_h_pv_phys = total_h_phys;
		sweet::Data::Cart2D::DataSpectral total_h_pv = total_h_phys(i_U.h_pert.cart2DDataConfig);
		total_h_pv_phys = _ops->avg_b_x(_ops->avg_b_y(total_h_phys));
		total_h_pv.loadCart2DDataGrid(total_h_pv_phys);

#if 0
		if (total_h_pv.reduce_min() < 0.00000001)
		{
			std::cerr << "Test case not adequate for vector invariant formulation. Null or negative water height" << std::endl;
			std::cerr << "Min h_pv   : " << total_h_pv.reduce_min() << std::endl;
			std::cerr << "Min h_total: " << total_h.reduce_min() << std::endl;
			std::cerr << "Min h_pert : " << i_h.reduce_min() << std::endl;
			SWEETErrorFatal("SWE_Cart2D_TS_ln_erk: Methods unstable or inadequate for vector invariant swe");;
		}
#endif

		sweet::Data::Cart2D::DataSpectral q = (_ops->diff_b_x(i_U.v) - _ops->diff_b_y(i_U.u) + _shackPDESWECart2D->cart2d_rotating_f0) / total_h_pv;
		sweet::Data::Cart2D::DataGrid q_phys = q.toGrid();

		// u, v tendencies
		// Energy conserving scheme
		sweet::Data::Cart2D::DataGrid o_u_t_phys = _ops->avg_f_y(q_phys*_ops->avg_b_x(V_phys)) - _ops->diff_b_x(H).toGrid();
		sweet::Data::Cart2D::DataGrid o_v_t_phys = -_ops->avg_f_x(q_phys*_ops->avg_b_y(U_phys)) - _ops->diff_b_y(H).toGrid();
		o_U.u.loadCart2DDataGrid(o_u_t_phys);
		o_U.v.loadCart2DDataGrid(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		if (!_shackPDESWECart2D->use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			o_U.h_pert = -_ops->diff_f_x(U) - _ops->diff_f_y(V);
		}
		else // use linear divergence
		{
			o_U.h_pert = -i_U.u*_ops->diff_f_x(total_h) - i_U.v*_ops->diff_f_y(total_h) + //nonlinear adv
					-_ops->diff_f_x(i_U.u*_shackPDESWECart2D->h0) - _ops->diff_f_y(i_U.v*_shackPDESWECart2D->h0); //linear div
			//o_h_t = -ops->diff_f_x(shackPDESWECart2D->h0*i_u) - ops->diff_f_y(shackPDESWECart2D->h0*i_v);
		}
	}

	return true;

}

}}}
