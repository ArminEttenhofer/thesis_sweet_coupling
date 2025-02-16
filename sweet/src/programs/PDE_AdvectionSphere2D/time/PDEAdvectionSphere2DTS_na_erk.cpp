/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_na_erk.hpp>




bool PDEAdvectionSphere2DTS_na_erk::testImplementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	return i_timestepping_method == "na_erk";
}

std::string PDEAdvectionSphere2DTS_na_erk::getStringId()
{
	return "na_erk";
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDEAdvectionSphere2DTS_na_erk::euler_timestep_update(
		const sweet::Data::Sphere2D::DataSpectral &i_prognostic_field,	//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_u,
		sweet::Data::Sphere2D::DataGrid &io_v,

		sweet::Data::Sphere2D::DataSpectral &o_prognostic_field,	//!< time updates

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
	sweet::Data::Sphere2D::DataSpectral phi = i_prognostic_field;

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (shackPDEAdvBenchmark->getVelocities)
	{
		shackPDEAdvBenchmark->getVelocities(io_u, io_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);

		sweet::Data::Sphere2D::DataSpectral vrt(phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div(phi.sphere2DDataConfig);
		ops->uv_2_vrtdiv(io_u, io_v, vrt, div);
	}

	sweet::Data::Sphere2D::DataGrid phig = phi.toGrid();

	sweet::Data::Sphere2D::DataGrid tmpg1 = io_u*phig;
	sweet::Data::Sphere2D::DataGrid tmpg2 = io_v*phig;

	o_prognostic_field = -ops->uv_2_div(tmpg1, tmpg2);
}



void PDEAdvectionSphere2DTS_na_erk::runTimestep(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_u,
		sweet::Data::Sphere2D::DataGrid &io_v,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	for (std::size_t i = 0; i < io_prognostic_fields.size(); i++)
	{
		// standard time stepping
		timestepping_rk.runTimestep_na(
				this,
				&PDEAdvectionSphere2DTS_na_erk::euler_timestep_update,	//!< pointer to function to compute euler time step updates
				io_prognostic_fields[i], io_u, io_v,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
}


bool PDEAdvectionSphere2DTS_na_erk::setup(
	sweet::Data::Sphere2D::Operators *io_ops
)
{
	PDEAdvectionSphere2DTS_BaseInterface::setup(io_ops);

	SWEET_ASSERT(shackPDEAdvectionTimeDisc != nullptr);
	timestepping_order = shackPDEAdvectionTimeDisc->timestepping_order;

	if (timestepping_order < 1 || timestepping_order > 4)
		return error.set("Invalid time stepping order");

	return true;
}


void PDEAdvectionSphere2DTS_na_erk::printImplementedTimesteppingMethods(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << " + PDEAdvectionSphere2DTS_na_erk:" << std::endl;
	o_ostream << i_prefix << "    * 'na_erk'" << std::endl;
}


PDEAdvectionSphere2DTS_na_erk::~PDEAdvectionSphere2DTS_na_erk()
{
}

