/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_erk.hpp"




bool PDESWESphere2DTS_lg_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphere2DTS_lg_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order	//!< order of RK time stepping method
)
{
	ops = io_ops;

	timestepping_order = i_order;
	return true;
}




void PDESWESphere2DTS_lg_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphere2DTS_lg_erk::euler_timestep_update,	//!< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}

/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_lg_erk::euler_timestep_update(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vort,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

	o_phi_t = -gh*i_div;
	o_div_t = -ops->laplace(i_phi);
	o_vort_t.spectral_setZero();
}



PDESWESphere2DTS_lg_erk::PDESWESphere2DTS_lg_erk()
{
}


PDESWESphere2DTS_lg_erk::~PDESWESphere2DTS_lg_erk()
{
}

