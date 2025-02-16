/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_nair_lauritzen_sl.hpp>
#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_na_trajectories.hpp>
#include "../PDEAdvectionSphere2DBenchmarksCombined.hpp"



bool PDEAdvectionSphere2DTS_na_trajectories::testImplementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "na_trajectories";
}

std::string PDEAdvectionSphere2DTS_na_trajectories::getStringId()
{
	return "na_trajectories";
}


void PDEAdvectionSphere2DTS_na_trajectories::printImplementedTimesteppingMethods(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << " + Sphere2DAdvection_TS_na_trajectories:" << std::endl;
	o_ostream << i_prefix << "    * 'na_trajectories'" << std::endl;
}



void PDEAdvectionSphere2DTS_na_trajectories::runTimestep(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &io_U_phi,		//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_U_u,
		sweet::Data::Sphere2D::DataGrid &io_U_v,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	for (std::size_t i = 0; i < io_U_phi.size(); i++)
	{
		run_timestep_1(
				io_U_phi[i],
				io_U_u,
				io_U_v,
				i_fixed_dt,
				i_simulation_timestamp
			);
	}
}

void PDEAdvectionSphere2DTS_na_trajectories::run_timestep_1(
		sweet::Data::Sphere2D::DataSpectral &io_U_phi,		//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_U_u,
		sweet::Data::Sphere2D::DataGrid &io_U_v,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi.sphere2DDataConfig;
	std::size_t N = sphere2DDataConfig->grid_number_elements;


	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (shackPDEAdvBenchmark->getVelocities != nullptr)
	{
		shackPDEAdvBenchmark->getVelocities(io_U_u, io_U_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);
	}


	sweet::Data::Vector::Vector<double> pos_lon_D(N), pos_lat_D(N);

	if (shackPDEAdvBenchmark->callback_slComputeDeparture3rdOrder != nullptr)
	{
		// compute 2nd order accurate departure points
		shackPDEAdvBenchmark->callback_slComputeDeparture3rdOrder(
				shackPDEAdvBenchmark->slComputeDeparture3rdOrderUserData,
				semiLagrangian.pos_lon_A,
				semiLagrangian.pos_lat_A,
				pos_lon_D,
				pos_lat_D,
				i_fixed_dt,
				i_simulation_timestamp+i_fixed_dt	// arrival time: current time + dt
			);
	}


	// sample phi at departure points
	sweet::Data::Sphere2D::DataGrid U_phi_phys_D =
		sphere2DSampler.bicubic_scalar_ret_phys(
			io_U_phi.toGrid(),
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			shackSemiLagrangian->semi_lagrangian_interpolation_limiter
		);

	io_U_phi = U_phi_phys_D;


	// sample phi at departure points

	U_phi_phys_D =
	sphere2DSampler.bicubic_scalar_ret_phys(
			io_U_phi.getSphere2DDataGrid(),
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			shackSemiLagrangian->semi_lagrangian_interpolation_limiter
		);


}



bool PDEAdvectionSphere2DTS_na_trajectories::setup(
	sweet::Data::Sphere2D::Operators *io_ops
)
{
	PDEAdvectionSphere2DTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDEAdvectionTimeDisc->timestepping_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		error.set("Only 1st and 2nd order for SL integration supported");

	semiLagrangian.setup(io_ops->sphere2DDataConfig, shackSemiLagrangian, timestepping_order);

	sphere2DSampler.setup(io_ops->sphere2DDataConfig);
	return true;
}


PDEAdvectionSphere2DTS_na_trajectories::~PDEAdvectionSphere2DTS_na_trajectories()
{
}

