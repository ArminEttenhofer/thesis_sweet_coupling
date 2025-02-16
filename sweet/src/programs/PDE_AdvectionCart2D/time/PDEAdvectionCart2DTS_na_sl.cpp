/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_na_sl.hpp>


PDEAdvectionCart2DTS_na_sl::PDEAdvectionCart2DTS_na_sl()
{
}

PDEAdvectionCart2DTS_na_sl::~PDEAdvectionCart2DTS_na_sl()
{
}


/*
 * Setup
 */
bool PDEAdvectionCart2DTS_na_sl::setup(sweet::Data::Cart2D::Operators *io_ops)
{
	PDEAdvectionCart2DTS_BaseInterface::setup(io_ops);

	timestepping_order = shackPDEAdvTimeDisc->timestepping_order;


	prog_u_prev.setup(ops->cart2DDataConfig);
	prog_v_prev.setup(ops->cart2DDataConfig);

	const sweet::Data::Cart2D::Config *cart2DDataConfig = ops->cart2DDataConfig;

	posx_a.setup(cart2DDataConfig->grid_number_elements);
	posy_a.setup(cart2DDataConfig->grid_number_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position
	posx_a.update_lambda_array_indices(
		[&](int idx, double &io_data)
		{
			int i = idx % cart2DDataConfig->grid_res[0];

			io_data = (double)i*(double)shackCart2DDataOps->cart2d_domain_size[0]/(double)cart2DDataConfig->grid_res[0];

			SWEET_ASSERT(io_data >= 0.0);
			SWEET_ASSERT(io_data <= shackCart2DDataOps->cart2d_domain_size[0]);
		}
	);
	posy_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
		{
			//int i = idx % cart2DDataConfig->grid_data_size[0];
			int j = idx / (double)cart2DDataConfig->grid_res[0];

			io_data = (double)j*(double)shackCart2DDataOps->cart2d_domain_size[1]/(double)cart2DDataConfig->grid_res[1];

			SWEET_ASSERT(io_data >= 0.0);
			SWEET_ASSERT(io_data <= shackCart2DDataOps->cart2d_domain_size[1]);
		}
	);

	// TODO: Use semiLagrangian.sampler2D
	sampler2D.setup(shackCart2DDataOps->cart2d_domain_size, cart2DDataConfig);

	//PXT- This just calls sampler2D.setup, so any reason for having it?
	semiLagrangian.setup(shackCart2DDataOps->cart2d_domain_size, cart2DDataConfig);

	return true;
}



void PDEAdvectionCart2DTS_na_sl::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_phi,		//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,		//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	SWEET_ASSERT(i_dt > 0);

	if (shackPDEAdvBenchmark->getExternalForcesCallback != nullptr)
	{
		shackPDEAdvBenchmark->getExternalForcesCallback(
				1,
				i_simulation_timestamp,
				&io_u,
				shackPDEAdvBenchmark
			);
		shackPDEAdvBenchmark->getExternalForcesCallback(
				2,
				i_simulation_timestamp,
				&io_v,
				shackPDEAdvBenchmark
			);
	}

	if (i_simulation_timestamp == 0)
	{
		prog_u_prev = io_u;
		prog_v_prev = io_v;
	}


	// OUTPUT: position of departure points at t
	sweet::Data::Vector::Vector<double> posx_d(io_phi.cart2DDataConfig->grid_number_elements);
	sweet::Data::Vector::Vector<double> posy_d(io_phi.cart2DDataConfig->grid_number_elements);

	semiLagrangian.semi_lag_departure_points_settls(
			prog_u_prev.toGrid(), prog_v_prev.toGrid(),
			io_u.toGrid(), io_v.toGrid(),
			posx_a, posy_a,
			i_dt,
			posx_d, posy_d,
			shackCart2DDataOps->cart2d_domain_size,
			nullptr,
			timestepping_order,

			shackPDEAdvTimeDisc->semi_lagrangian_max_iterations,
			shackPDEAdvTimeDisc->semi_lagrangian_convergence_threshold
	);

	prog_u_prev = io_u;
	prog_v_prev = io_v;

	sweet::Data::Cart2D::DataSpectral new_prog_phi(io_phi.cart2DDataConfig);

	if (timestepping_order == 1)
	{
		sampler2D.bilinear_scalar(
				io_phi,
				posx_d,
				posy_d,
				new_prog_phi
		);
	}
	else if (timestepping_order == 2)
	{
		sampler2D.bicubic_scalar(
				io_phi,
				posx_d,
				posy_d,
				new_prog_phi
		);
	}
	else
	{
		SWEETErrorFatal("Timestepping order not available");
	}

	io_phi = new_prog_phi;
}



