/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_na_sl.hpp>


bool PDEAdvectionSphere2DTS_na_sl::testImplementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "na_sl";
}

std::string PDEAdvectionSphere2DTS_na_sl::getStringId()
{
	return "na_sl";
}



bool PDEAdvectionSphere2DTS_na_sl::setup(
	sweet::Data::Sphere2D::Operators *io_ops
)
{
	PDEAdvectionSphere2DTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDEAdvectionTimeDisc->timestepping_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		error.set("Only 1st and 2nd order for SL integration supported");

	SWEET_ASSERT(shackSemiLagrangian != nullptr);

	semiLagrangian.setup(
			ops->sphere2DDataConfig,
			shackSemiLagrangian,
			timestepping_order
		);

	sphere2DSampler.setup(io_ops->sphere2DDataConfig);
	return true;
}


void PDEAdvectionSphere2DTS_na_sl::printImplementedTimesteppingMethods(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << " + Sphere2DAdvection_TS_na_sl:" << std::endl;
	o_ostream << "    * 'na_sl'" << std::endl;
}




/*
 * SL treatment of 3D Vector in Cartesian space
 */
void PDEAdvectionSphere2DTS_na_sl::interpolate_departure_point_vec_3d(
		const sweet::Data::Sphere2D::DataSpectral &i_vec0,
		const sweet::Data::Sphere2D::DataSpectral &i_vec1,
		const sweet::Data::Sphere2D::DataSpectral &i_vec2,

		const sweet::Data::Vector::Vector<double> &i_pos_lon_D,
		const sweet::Data::Vector::Vector<double> &i_pos_lat_D,

		sweet::Data::Sphere2D::DataSpectral &o_vec0,
		sweet::Data::Sphere2D::DataSpectral &o_vec1,
		sweet::Data::Sphere2D::DataSpectral &o_vec2
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = i_vec0.sphere2DDataConfig;

	o_vec0.setup_if_required(i_vec0.sphere2DDataConfig);
	o_vec1.setup_if_required(i_vec1.sphere2DDataConfig);
	o_vec2.setup_if_required(i_vec2.sphere2DDataConfig);

	/*
	 * First we sample the field at the departure point
	 */

	bool velocity_sampling = false;

	sweet::Data::Sphere2D::DataGrid u_tmp_D;
	sweet::Data::Sphere2D::DataGrid v_tmp_D;
	sweet::Data::Sphere2D::DataGrid w_tmp_D;

	if (timestepping_order == 1 && false)
	{
		u_tmp_D = sphere2DSampler.bilinear_scalar_ret_phys(
				o_vec0.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				velocity_sampling,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points
			);

		v_tmp_D = sphere2DSampler.bilinear_scalar_ret_phys(
				o_vec1.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				velocity_sampling,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points
			);

		w_tmp_D = sphere2DSampler.bilinear_scalar_ret_phys(
				o_vec2.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				velocity_sampling,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points
			);
	}
	else
	{
		u_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				o_vec0.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				velocity_sampling,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSemiLagrangian->semi_lagrangian_interpolation_limiter
			);

		v_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				o_vec1.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				velocity_sampling,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSemiLagrangian->semi_lagrangian_interpolation_limiter
			);

		w_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				o_vec2.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				velocity_sampling,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSemiLagrangian->semi_lagrangian_interpolation_limiter
			);
	}

	/*
	 * Convert to scalar data arrays
	 */
	sweet::Data::Vector::Vector<double> V_u_D = sweet::Data::Sphere2D::Convert::DataGrid_2_Vector::convert(u_tmp_D);
	sweet::Data::Vector::Vector<double> V_v_D = sweet::Data::Sphere2D::Convert::DataGrid_2_Vector::convert(v_tmp_D);
	sweet::Data::Vector::Vector<double> V_w_D = sweet::Data::Sphere2D::Convert::DataGrid_2_Vector::convert(w_tmp_D);

#if 1
	/*
	 * Here we have the velocity at the departure points.
	 *
	 * Now we need to take the rotation of this vector into account!
	 */

	/*
	 * Compute departure position
	 */
	sweet::Data::Vector::Vector<double> P_x_D, P_y_D, P_z_D;
	sweet::LibMath::VectorMath::point_latlon_2_cartesian__array(i_pos_lon_D, i_pos_lat_D, P_x_D, P_y_D, P_z_D);

	sweet::Data::Vector::Vector<double> &P_x_A = semiLagrangian.pos_x_A;
	sweet::Data::Vector::Vector<double> &P_y_A = semiLagrangian.pos_y_A;
	sweet::Data::Vector::Vector<double> &P_z_A = semiLagrangian.pos_z_A;

	/*
	 * Compute rotation angle based on departure and arrival position
	 */
	sweet::Data::Vector::Vector<double> rotation_angle_ =
			sweet::LibMath::VectorMath::dot_prod(
				P_x_D, P_y_D, P_z_D,
				P_x_A, P_y_A, P_z_A
			);

	// Can be slightly larger than 1 due to round-off issues, leading to NaN, hence this hack
	rotation_angle_ = sweet::LibMath::VectorMath::min(rotation_angle_, 1.0);

	/*
	 * Compute rotation angle
	 */
	sweet::Data::Vector::Vector<double> rotation_angle = sweet::LibMath::VectorMath::arccos(rotation_angle_);

	/*
	 * Compute Rotation axis and normalize
	 */
	sweet::Data::Vector::Vector<double> rot_x, rot_y, rot_z;
	sweet::LibMath::VectorMath::cross_prod(
			P_x_D, P_y_D, P_z_D,
			P_x_A, P_y_A, P_z_A,
			rot_x, rot_y, rot_z
		);
	sweet::LibMath::VectorMath::normalize_with_threshold(rot_x, rot_y, rot_z);

	/*
	 * Rotate vector (using transpose of rotation matrix without translation!)
	 */
	sweet::Data::Vector::Vector<double> V_u_A, V_v_A, V_w_A;
	sweet::LibMath::VectorMath::point_rotate_3d_normalized_rotation_axis__array(
			V_u_D, V_v_D, V_w_D,
			rotation_angle,
			rot_x, rot_y, rot_z,
			V_u_A, V_v_A, V_w_A
		);

	o_vec0 = sweet::Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_u_A, sphere2DDataConfig);
	o_vec1 = sweet::Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_v_A, sphere2DDataConfig);
	o_vec2 = sweet::Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_w_A, sphere2DDataConfig);

#else

	o_vec0 = Convert_Vector_2_Sphere2DDataGrid::convert(V_u_D, sphere2DDataConfig);
	o_vec1 = Convert_Vector_2_Sphere2DDataGrid::convert(V_v_D, sphere2DDataConfig);
	o_vec2 = Convert_Vector_2_Sphere2DDataGrid::convert(V_w_D, sphere2DDataConfig);

#endif
}




/*
 * SL treatment of 2D Vector in UV space
 */
void PDEAdvectionSphere2DTS_na_sl::interpolate_departure_point_vec_uv(
		const sweet::Data::Sphere2D::DataGrid &i_u,
		const sweet::Data::Sphere2D::DataGrid &i_v,

		const sweet::Data::Vector::Vector<double> &i_pos_lon_D,
		const sweet::Data::Vector::Vector<double> &i_pos_lat_D,

		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = i_u.sphere2DDataConfig;

	o_u.setup_if_required(i_u.sphere2DDataConfig);
	o_v.setup_if_required(i_v.sphere2DDataConfig);


	/*********************************************************************
	 * Step 1)
	 * Sample velocity at departure points and convert to Cartesian space
	 *********************************************************************/

	/*
	 * First we sample the field at the departure point
	 */

	sweet::Data::Sphere2D::DataGrid u_tmp_D;
	sweet::Data::Sphere2D::DataGrid v_tmp_D;

	if (timestepping_order == 1 && false)
	{
		u_tmp_D = sphere2DSampler.bilinear_scalar_ret_phys(
				i_u,
				i_pos_lon_D, i_pos_lat_D,
				true,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points
			);

		v_tmp_D = sphere2DSampler.bilinear_scalar_ret_phys(
				i_v,
				i_pos_lon_D, i_pos_lat_D,
				true,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points
			);
	}
	else
	{
		u_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				i_u,
				i_pos_lon_D, i_pos_lat_D,
				true,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSemiLagrangian->semi_lagrangian_interpolation_limiter
			);

		v_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				i_v,
				i_pos_lon_D, i_pos_lat_D,
				true,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSemiLagrangian->semi_lagrangian_interpolation_limiter
			);
	}

	/*
	 * Convert to Cartesian space
	 */
	sweet::Data::Vector::Vector<double> V_lon_D = sweet::Data::Sphere2D::Convert::DataGrid_2_Vector::convert(u_tmp_D);
	sweet::Data::Vector::Vector<double> V_lat_D = sweet::Data::Sphere2D::Convert::DataGrid_2_Vector::convert(v_tmp_D);

	/*
	 * Convert from UV Velocity space to 3D Cartesian space
	 */
	sweet::Data::Vector::Vector<double> V_x_D, V_y_D, V_z_D;
	sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
			i_pos_lon_D,
			i_pos_lat_D,
			V_lon_D,
			V_lat_D,
			V_x_D,
			V_y_D,
			V_z_D
		);

#if 1
	/*********************************************************************
	 * Step 2)
	 * Prepare rotation
	 *********************************************************************
	 * Here we have the velocity at the departure points.
	 *
	 * Now we need to take the rotation of this vector into account!
	 */

	/*
	 * Compute departure position
	 */
	sweet::Data::Vector::Vector<double> P_x_D, P_y_D, P_z_D;
	sweet::LibMath::VectorMath::point_latlon_2_cartesian__array(i_pos_lon_D, i_pos_lat_D, P_x_D, P_y_D, P_z_D);

	sweet::Data::Vector::Vector<double> &P_x_A = semiLagrangian.pos_x_A;
	sweet::Data::Vector::Vector<double> &P_y_A = semiLagrangian.pos_y_A;
	sweet::Data::Vector::Vector<double> &P_z_A = semiLagrangian.pos_z_A;

	/*
	 * Compute rotation angle based on departure and arrival position
	 */
	sweet::Data::Vector::Vector<double> rotation_angle_ =
			sweet::LibMath::VectorMath::dot_prod(
				P_x_D, P_y_D, P_z_D,
				P_x_A, P_y_A, P_z_A
			);

	// Can be slightly larger than 1 due to round-off issues, leading to NaN, hence this hack
	rotation_angle_ = sweet::LibMath::VectorMath::min(rotation_angle_, 1.0);

	/*
	 * Compute rotation angle
	 */
	sweet::Data::Vector::Vector<double> rotation_angle = sweet::LibMath::VectorMath::arccos(rotation_angle_);

	/*
	 * Compute Rotation axis and normalize
	 */
	sweet::Data::Vector::Vector<double> rot_x, rot_y, rot_z;
	sweet::LibMath::VectorMath::cross_prod(
			P_x_D, P_y_D, P_z_D,
			P_x_A, P_y_A, P_z_A,
			rot_x, rot_y, rot_z
		);
	sweet::LibMath::VectorMath::normalize_with_threshold(rot_x, rot_y, rot_z);


	/*
	 * Rotate vector
	 */
	sweet::Data::Vector::Vector<double> V_x_A, V_y_A, V_z_A;
	sweet::LibMath::VectorMath::point_rotate_3d_normalized_rotation_axis__array(
			V_x_D, V_y_D, V_z_D,
			rotation_angle,
			rot_x, rot_y, rot_z,
			V_x_A, V_y_A, V_z_A
		);

#else

	sweet::Data::Vector::Vector<double> V_x_A = V_x_D;
	sweet::Data::Vector::Vector<double> V_y_A = V_y_D;
	sweet::Data::Vector::Vector<double> V_z_A = V_z_D;

#endif

	/*********************************************************************
	 * Step 3)
	 * Convert Velocity from Cartesian to u/v space
	 *********************************************************************/

	/*
	 * Return velocity in lat/lon space
	 */
	sweet::Data::Vector::Vector<double> V_lon_A, V_lat_A;
	sweet::LibMath::VectorMath::velocity_cartesian_2_latlon__array(
			semiLagrangian.pos_lon_A,
			semiLagrangian.pos_lat_A,
			V_x_A,
			V_y_A,
			V_z_A,
			V_lon_A, V_lat_A
	);

	o_u = sweet::Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_lon_A, sphere2DDataConfig);
	o_v = sweet::Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_lat_A, sphere2DDataConfig);
}



void PDEAdvectionSphere2DTS_na_sl::run_timestep_1(
		sweet::Data::Sphere2D::DataSpectral &io_U_phi,		//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_U_u,
		sweet::Data::Sphere2D::DataGrid &io_U_v,

		double i_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi.sphere2DDataConfig;

	if (i_simulation_timestamp == 0)
	{
		U_u_prev = io_U_u;
		U_v_prev = io_U_v;
	}

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */

	if (shackPDEAdvBenchmark->getVelocities)
	{
		shackPDEAdvBenchmark->getVelocities(io_U_u, io_U_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);
		shackPDEAdvBenchmark->getVelocities(U_u_prev, U_v_prev, i_simulation_timestamp - i_dt, shackPDEAdvBenchmark->getVelocitiesUserData);
	}

	// OUTPUT: position of departure points at t
	sweet::Data::Vector::Vector<double> pos_lon_D(sphere2DDataConfig->grid_number_elements);
	sweet::Data::Vector::Vector<double> pos_lat_D(sphere2DDataConfig->grid_number_elements);

	double dt_div_radius = shackTimestepControl->currentTimestepSize / shackSphere2DDataOps->sphere2d_radius;

	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*io_U_u, dt_div_radius*io_U_v,

			pos_lon_D, pos_lat_D
	);

	U_u_prev = io_U_u;
	U_v_prev = io_U_v;

	sweet::Data::Sphere2D::DataGrid new_prog_phi_phys;

	if (timestepping_order == 1 && false)
	{
		new_prog_phi_phys =
			sphere2DSampler.bilinear_scalar_ret_phys(
				io_U_phi.toGrid(),
				pos_lon_D,
				pos_lat_D,
				false,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points
		);
	}
	else
	{
		new_prog_phi_phys =
			sphere2DSampler.bicubic_scalar_ret_phys(
				io_U_phi.toGrid(),
				pos_lon_D,
				pos_lat_D,
				false,
				shackSemiLagrangian->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSemiLagrangian->semi_lagrangian_interpolation_limiter
		);
	}

	io_U_phi = new_prog_phi_phys;
}


void PDEAdvectionSphere2DTS_na_sl::run_timestep_2(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_U_u,
		sweet::Data::Sphere2D::DataGrid &io_U_v,

		double i_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_prognostic_fields[0].sphere2DDataConfig;

	if (i_simulation_timestamp == 0)
	{
		U_u_prev = io_U_u;
		U_v_prev = io_U_v;
	}

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (shackPDEAdvBenchmark->getVelocities)
	{
		shackPDEAdvBenchmark->getVelocities(io_U_u, io_U_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);
		shackPDEAdvBenchmark->getVelocities(U_u_prev, U_v_prev, i_simulation_timestamp - i_dt, shackPDEAdvBenchmark->getVelocitiesUserData);
	}

	// OUTPUT: position of departure points at t
	sweet::Data::Vector::Vector<double> pos_lon_d(sphere2DDataConfig->grid_number_elements);
	sweet::Data::Vector::Vector<double> pos_lat_d(sphere2DDataConfig->grid_number_elements);

	double dt_div_radius = shackTimestepControl->currentTimestepSize / shackSphere2DDataOps->sphere2d_radius;

	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*io_U_u, dt_div_radius*io_U_v,

			pos_lon_d, pos_lat_d
	);


	sweet::Data::Sphere2D::DataGrid u, v;
	ops->vrtdiv_2_uv(
			io_prognostic_fields[0], io_prognostic_fields[1],
			u, v
		);

	sweet::Data::Sphere2D::DataGrid new_u(sphere2DDataConfig), new_v(sphere2DDataConfig);
	interpolate_departure_point_vec_uv(
			u, v,

			pos_lon_d,
			pos_lat_d,

			new_u, new_v
	);

	ops->uv_2_vrtdiv(
			new_u, new_v,
			io_prognostic_fields[0], io_prognostic_fields[1]
		);

}



void PDEAdvectionSphere2DTS_na_sl::run_timestep_3(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_U_u,
		sweet::Data::Sphere2D::DataGrid &io_U_v,

		double i_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_prognostic_fields[0].sphere2DDataConfig;

	if (i_simulation_timestamp == 0)
	{
		U_u_prev = io_U_u;
		U_v_prev = io_U_v;
	}

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (shackPDEAdvBenchmark->getVelocities)
	{
		shackPDEAdvBenchmark->getVelocities(io_U_u, io_U_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);
		shackPDEAdvBenchmark->getVelocities(U_u_prev, U_v_prev, i_simulation_timestamp - i_dt, shackPDEAdvBenchmark->getVelocitiesUserData);
	}

	// OUTPUT: position of departure points at t
	sweet::Data::Vector::Vector<double> pos_lon_d(sphere2DDataConfig->grid_number_elements);
	sweet::Data::Vector::Vector<double> pos_lat_d(sphere2DDataConfig->grid_number_elements);

	double dt_div_radius = shackTimestepControl->currentTimestepSize / shackSphere2DDataOps->sphere2d_radius;

	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*io_U_u, dt_div_radius*io_U_v,

			pos_lon_d, pos_lat_d
	);



	interpolate_departure_point_vec_3d(
			io_prognostic_fields[0],
			io_prognostic_fields[1],
			io_prognostic_fields[2],

			pos_lon_d,
			pos_lat_d,

			io_prognostic_fields[0],
			io_prognostic_fields[1],
			io_prognostic_fields[2]
	);
}



void PDEAdvectionSphere2DTS_na_sl::runTimestep(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
		sweet::Data::Sphere2D::DataGrid &io_U_u,
		sweet::Data::Sphere2D::DataGrid &io_U_v,

		double i_dt,
		double i_simulation_timestamp
)
{
	if (io_prognostic_fields.size() == 1)
	{
		run_timestep_1(
				io_prognostic_fields[0],
				io_U_u,
				io_U_v,
				i_dt,
				i_simulation_timestamp
			);
		return;
	}

	if (io_prognostic_fields.size() == 2)
	{
		run_timestep_2(
				io_prognostic_fields,
				io_U_u,
				io_U_v,
				i_dt,
				i_simulation_timestamp
			);
		return;
	}


	if (io_prognostic_fields.size() == 3)
	{
		run_timestep_3(
				io_prognostic_fields,
				io_U_u,
				io_U_v,
				i_dt,
				i_simulation_timestamp
			);
		return;
	}

	SWEETErrorFatal("Should never happen");
}


PDEAdvectionSphere2DTS_na_sl::~PDEAdvectionSphere2DTS_na_sl()
{
}

