/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_na_erk.hpp>


PDEAdvectionCart2DTS_na_erk::PDEAdvectionCart2DTS_na_erk()
{
}


PDEAdvectionCart2DTS_na_erk::~PDEAdvectionCart2DTS_na_erk()
{
}


/*
 * Setup
 */
bool PDEAdvectionCart2DTS_na_erk::setup(sweet::Data::Cart2D::Operators *io_ops)
{
	PDEAdvectionCart2DTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDEAdvTimeDisc->timestepping_order;
	return true;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDEAdvectionCart2DTS_na_erk::euler_timestep_update(
		const sweet::Data::Cart2D::DataSpectral &i_phi,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_phi_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

		double i_simulation_timestamp
)
{
	/**
	 * We simply compute
	 * 	-DIV(rho*U) = -rho DIV(U) - U.GRAD(rho) = - U.GRAD(rho)
	 * which is the Lagrangian contribution only.
	 *
	 * This is the case because the velocity field is divergence free!!!
	 */

	if (shackPDEAdvBenchmark->getExternalForcesCallback != nullptr)
	{
		sweet::Data::Cart2D::DataSpectral u(i_phi.cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral v(i_phi.cart2DDataConfig);

		shackPDEAdvBenchmark->getExternalForcesCallback(
				1,
				shackTimestepControl->currentSimulationTime,
				&u,
				shackPDEAdvBenchmark
			);
		shackPDEAdvBenchmark->getExternalForcesCallback(
				2,
				shackTimestepControl->currentSimulationTime,
				&v,
				shackPDEAdvBenchmark
			);

		o_phi_t = -ops->diff_c_x(i_phi*u) - ops->diff_c_y(i_phi*v);
	}
	else
	{
		o_phi_t = -ops->diff_c_x(i_phi*i_u) - ops->diff_c_y(i_phi*i_v);
	}

	o_u_t.spectral_setZero();
	o_v_t.spectral_setZero();
}


void PDEAdvectionCart2DTS_na_erk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_phi,		//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,		//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	SWEET_ASSERT(i_dt > 0);

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDEAdvectionCart2DTS_na_erk::euler_timestep_update,	//!< pointer to function to compute euler time step updates
			io_phi, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);

	if (shackPDEAdvBenchmark->getExternalForcesCallback != nullptr)
	{
		// this is just called for cosmetic reasons to update the velocity field
		shackPDEAdvBenchmark->getExternalForcesCallback(
				1,
				shackTimestepControl->currentSimulationTime+i_dt,
				&io_u,
				shackPDEAdvBenchmark
			);
		shackPDEAdvBenchmark->getExternalForcesCallback(
				2,
				shackTimestepControl->currentSimulationTime+i_dt,
				&io_v,
				shackPDEAdvBenchmark
			);
	}
}


