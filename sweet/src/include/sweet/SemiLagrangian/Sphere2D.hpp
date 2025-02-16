/*
 * Sphere2DDataSemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Updated to sphere2D on 28th March 2018
 */
#ifndef INCLUDE_SWEET_SEMILAGRANGIAN_SPHERE2D_HPP
#define INCLUDE_SWEET_SEMILAGRANGIAN_SPHERE2D_HPP

#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Vector.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <sweet/Data/Sphere2D/Operators_Sampler_Sphere2DDataGrid.hpp>
#include <sweet/Data/Vector/Convert/Vector_2_Sphere2D_DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <limits>
#include <sweet/LibMath/VectorMath.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>


namespace sweet {
namespace SemiLagrangian {

class Sphere2D
{
	const Data::Sphere2D::Config *sphere2DDataConfig;
	sweet::SemiLagrangian::Shack *shackSLData;

	Data::Sphere2D::DataGrid sl_coriolis;

	// Arrival points
public:
	Data::Vector::Vector<double> pos_lat_A;
	Data::Vector::Vector<double> pos_lon_A;

	Data::Vector::Vector<double> pos_x_A;
	Data::Vector::Vector<double> pos_y_A;
	Data::Vector::Vector<double> pos_z_A;

	Data::Sphere2D::Operators_Sampler_DataGrid sphere2DSampler;


	int timestepping_order;
	int semi_lagrangian_max_iterations;
	double semi_lagrangian_convergence_threshold;
	bool semi_lagrangian_approximate_sphere2d_geometry;
	bool semi_lagrangian_interpolation_limiter;

	enum EnumTrajectories
	{
		E_TRAJECTORY_METHOD_CANONICAL,
		E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE,
		E_TRAJECTORY_METHOD_SETTLS_HORTAL
	};


	EnumTrajectories trajectory_method;


	Sphere2D()	:
		sphere2DDataConfig(nullptr),
		shackSLData(nullptr)
	{
	}


	void setup(
			const Data::Sphere2D::Config *i_sphere2DDataConfig,
			Shack *i_shackSLData,
			int i_timestepping_order
	)
	{
		shackSLData = i_shackSLData;
		sphere2DDataConfig = i_sphere2DDataConfig;
		timestepping_order = i_timestepping_order;

		sphere2DSampler.setup(sphere2DDataConfig);

		if (shackSLData->semi_lagrangian_departure_point_method == "settls")
		{
			trajectory_method = E_TRAJECTORY_METHOD_SETTLS_HORTAL;
		}
		else if (shackSLData->semi_lagrangian_departure_point_method == "canonical")
		{
			trajectory_method = E_TRAJECTORY_METHOD_CANONICAL;
		}
		else if (shackSLData->semi_lagrangian_departure_point_method == "midpoint")
		{
			trajectory_method = E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE;
		}
		else
		{
			SWEETErrorFatal(std::string("Trajectory method '")+shackSLData->semi_lagrangian_departure_point_method+"' not supported");
		}


		semi_lagrangian_max_iterations = shackSLData->semi_lagrangian_max_iterations;
		semi_lagrangian_convergence_threshold = shackSLData->semi_lagrangian_convergence_threshold;
		semi_lagrangian_approximate_sphere2d_geometry = shackSLData->semi_lagrangian_approximate_sphere2d_geometry;
		semi_lagrangian_interpolation_limiter = shackSLData->semi_lagrangian_interpolation_limiter;



		pos_lon_A.setup(sphere2DDataConfig->grid_number_elements);
		pos_lon_A.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % sphere2DDataConfig->grid_num_lon;

				io_data = 2.0*M_PI*(double)i/(double)sphere2DDataConfig->grid_num_lon;
				SWEET_ASSERT(io_data >= 0);
				SWEET_ASSERT(io_data < 2.0*M_PI);
			}
		);

		pos_lat_A.setup(sphere2DDataConfig->grid_number_elements);
		pos_lat_A.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				int j = idx / sphere2DDataConfig->grid_num_lon;

				io_data = sphere2DDataConfig->lat[j];

				SWEET_ASSERT(io_data >= -M_PI*0.5);
				SWEET_ASSERT(io_data <= M_PI*0.5);
			}
		);

		pos_x_A.setup(i_sphere2DDataConfig->grid_number_elements);
		pos_y_A.setup(i_sphere2DDataConfig->grid_number_elements);
		pos_z_A.setup(i_sphere2DDataConfig->grid_number_elements);

		sweet::LibMath::VectorMath::point_latlon_2_cartesian__array(
				pos_lon_A, pos_lat_A,
				pos_x_A, pos_y_A, pos_z_A
			);


#if 0
		//Sphere2DData_Grid u_lon_
		if (simVars.sim.sphere2d_use_fsphere2D)
		{
			Error.set("Not supported");
		}

		sl_coriolis.setup_if_required(sphere2DDataConfig);
		sl_coriolis.grid_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu*2.0*simVars.sim.sphere2d_rotating_coriolis_omega*simVars.sim.sphere2d_radius;
			}
		);
#endif

	}



	/**
	 * Do 1st order accurate advection on the sphere2D for given
	 *
	 *  - starting point,
	 *  - velocity*dt
	 *
	 * in Cartesian space 
	 */
	inline static
	void doAdvectionOnSphere2D(
		const Data::Vector::Vector<double> &i_pos_x,
		const Data::Vector::Vector<double> &i_pos_y,
		const Data::Vector::Vector<double> &i_pos_z,

		const Data::Vector::Vector<double> &i_dt_velocity_x,
		const Data::Vector::Vector<double> &i_dt_velocity_y,
		const Data::Vector::Vector<double> &i_dt_velocity_z,

		Data::Vector::Vector<double> &o_pos_x,
		Data::Vector::Vector<double> &o_pos_y,
		Data::Vector::Vector<double> &o_pos_z,

		double i_approximate_sphere2d_geometry
	)
	{
		if (i_approximate_sphere2d_geometry)
		{
//			std::cout << "i_approximate_sphere2d_geometry" << std::endl;
			/*
			 * This just uses an approximation of the sphere2D geometry.
			 */

			/*
			 * Step 1) Apply the velocity vector in Cartesian space, ignoring the sphere2D's curvature
			 */
			o_pos_x = i_pos_x + i_dt_velocity_x;
			o_pos_y = i_pos_y + i_dt_velocity_y;
			o_pos_z = i_pos_z + i_dt_velocity_z;

			/*
			 * Step 2) Normalize the resulting position
			 */
			sweet::LibMath::VectorMath::normalize(o_pos_x, o_pos_y, o_pos_z);
			return;
		}


		/*
		 * This version implements an accurate geometry.
		 */
		{
			/*
			 * Step 1) Compute rotation axis with cross product
			 */
			Data::Vector::Vector<double> rotation_axis_x(i_pos_x.numberOfElements);
			Data::Vector::Vector<double> rotation_axis_y(i_pos_x.numberOfElements);
			Data::Vector::Vector<double> rotation_axis_z(i_pos_x.numberOfElements);

			sweet::LibMath::VectorMath::cross_prod(
					i_pos_x, i_pos_y, i_pos_z,
					i_dt_velocity_x, i_dt_velocity_y, i_dt_velocity_z,
					rotation_axis_x, rotation_axis_y, rotation_axis_z
				);
#if 1
			/*
			 * Normalize rotation axis since it's likely not normalized yet
			 */

			static int asdf = 1;

			if (asdf == 1)
			{
				std::cout << "TODO: REPLACE ME" << std::endl;
				asdf = 0;
			}

			/*
			 * TODO: replace this!
			 *
			 * Use a formulation without normalization by using the
			 * angular vector / velocity e.g. by using quaternions
			 */
			sweet::LibMath::VectorMath::normalize_with_threshold(
					rotation_axis_x,
					rotation_axis_y,
					rotation_axis_z
				);

			/*
			 * Compute rotation angle
			 *
			 * No rescaling by 1/(2pi) since the change of angle is directly
			 * given by the magnitude of the angular velocity by its definition
			 */
			Data::Vector::Vector<double> angle = sweet::LibMath::VectorMath::length(i_dt_velocity_x, i_dt_velocity_y, i_dt_velocity_z);
#else
			// doesn't work
			Data::Vector::Vector<double> angle = i_dt_velocity_x;
			angle.grid_set_all(1.0);
#endif
			/*
			 * Rotate
			 */
			sweet::LibMath::VectorMath::point_rotate_3d_normalized_rotation_axis__array(
					i_pos_x,
					i_pos_y,
					i_pos_z,
					angle,
					rotation_axis_x,
					rotation_axis_y,
					rotation_axis_z,
					o_pos_x,
					o_pos_y,
					o_pos_z
				);
		}
	}



	/*!
	 * Compute SL departure points on unit sphere2D for given dt*(u,v) velocities
	 *
	 * All this is for the unit sphere2D and unit time!
	 *
	 * Hence, the velocities need to be rescaled by "dt/radius"
	 */
	void semi_lag_departure_points_settls_specialized(
		const Data::Sphere2D::DataGrid &i_u_lon_prev,	//!< Velocities at time t-1
		const Data::Sphere2D::DataGrid &i_v_lat_prev,

		const Data::Sphere2D::DataGrid &i_dt_u_lon, 		//!< Velocities at time t
		const Data::Sphere2D::DataGrid &i_dt_v_lat,

		Data::Vector::Vector<double> &o_pos_lon_D, 	//!< OUTPUT: Position of departure points x / y
		Data::Vector::Vector<double> &o_pos_lat_D
	)
	{
		o_pos_lon_D.setup_if_required(pos_lon_A);
		o_pos_lat_D.setup_if_required(pos_lon_A);

		std::size_t num_elements = o_pos_lon_D.numberOfElements;

		if (timestepping_order == 1)
		{
			/*
			 * Compute Cartesian velocity
			 */
			Data::Vector::Vector<double> u_lon_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(i_dt_u_lon);
			Data::Vector::Vector<double> v_lat_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(i_dt_v_lat);

			Data::Vector::Vector<double> vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
			sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
					pos_lon_A, pos_lat_A,
					u_lon_array, v_lat_array,
					vel_x_A, vel_y_A, vel_z_A
				);

			/*
			 * Do advection in Cartesian space
			 */

			Data::Vector::Vector<double> new_pos_x_d(num_elements), new_pos_y_d(num_elements), new_pos_z_d(num_elements);
			doAdvectionOnSphere2D(
				pos_x_A, pos_y_A, pos_z_A,
				-vel_x_A, -vel_y_A, -vel_z_A,
				new_pos_x_d, new_pos_y_d, new_pos_z_d,

				semi_lagrangian_approximate_sphere2d_geometry
			);

			/*
			 * Departure point to lat/lon coordinate
			 */
			sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(
					new_pos_x_d, new_pos_y_d, new_pos_z_d,
					o_pos_lon_D, o_pos_lat_D
				);
			return;
		}


		if (timestepping_order == 2)
		{
			if (trajectory_method == E_TRAJECTORY_METHOD_CANONICAL)
			{
				/*
				 * Standard iterative method
				 *
				 * See also Michail Diamantarkis paper, p. 185
				 */

				/*
				 * Prepare
				 */
				const Data::Sphere2D::DataGrid &vel_lon = i_dt_u_lon;
				const Data::Sphere2D::DataGrid &vel_lat = i_dt_v_lat;

				Data::Vector::Vector<double> vel_lon_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(vel_lon);
				Data::Vector::Vector<double> vel_lat_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(vel_lat);

				/*
				 * Polar => Cartesian velocities
				 */
				Data::Vector::Vector<double> vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Step 1)
				 */
				// Departure points for iterations
				Data::Vector::Vector<double> pos_x_D(sphere2DDataConfig->grid_number_elements);
				Data::Vector::Vector<double> pos_y_D(sphere2DDataConfig->grid_number_elements);
				Data::Vector::Vector<double> pos_z_D(sphere2DDataConfig->grid_number_elements);

				doAdvectionOnSphere2D(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-vel_x_A,
					-vel_y_A,
					-vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					semi_lagrangian_approximate_sphere2d_geometry
				);

				/*
				 * 2 iterations to get midpoint
				 */
				double diff = -1;
				for (int i = 0; i < semi_lagrangian_max_iterations; i++)
				{
					/*
					 * Step 2a
					 */
					Data::Vector::Vector<double> pos_x_mid = 0.5*(pos_x_A + pos_x_D);
					Data::Vector::Vector<double> pos_y_mid = 0.5*(pos_y_A + pos_y_D);
					Data::Vector::Vector<double> pos_z_mid = 0.5*(pos_z_A + pos_z_D);

					Data::Vector::Vector<double> pos_lon_mid(sphere2DDataConfig->grid_number_elements);
					Data::Vector::Vector<double> pos_lat_mid(sphere2DDataConfig->grid_number_elements);

					sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(
							pos_x_mid, pos_y_mid, pos_z_mid,
							pos_lon_mid, pos_lat_mid
						);

					Data::Vector::Vector<double> vel_u_mid = sphere2DSampler.bilinear_scalar(vel_lon, pos_lon_mid, pos_lat_mid, true, shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points);
					Data::Vector::Vector<double> vel_v_mid = sphere2DSampler.bilinear_scalar(vel_lat, pos_lon_mid, pos_lat_mid, true, shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					Data::Vector::Vector<double> vel_x_mid(num_elements), vel_y_mid(num_elements), vel_z_mid(num_elements);
					sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
							pos_lon_mid, pos_lat_mid,
							vel_u_mid, vel_v_mid,
							vel_x_mid, vel_y_mid, vel_z_mid
						);

					// convert final points from Cartesian space to angular space
					Data::Vector::Vector<double> new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					/*
					 * Step 2b
					 */
					doAdvectionOnSphere2D(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-vel_x_mid,
						-vel_y_mid,
						-vel_z_mid,

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						semi_lagrangian_approximate_sphere2d_geometry
					);


					if (semi_lagrangian_convergence_threshold > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
										(pos_y_D-new_pos_y_D).reduce_maxAbs() +
										(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < semi_lagrangian_convergence_threshold)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (semi_lagrangian_convergence_threshold > 0)
				{
					if (diff > semi_lagrangian_convergence_threshold)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << semi_lagrangian_convergence_threshold << std::endl;
					}
				}

				// convert final points from Cartesian space to angular space
				sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else if (trajectory_method == E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE)
			{

				/*
				 * Ritchies midpoint rule
				 */
				Data::Vector::Vector<double> vel_lon_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(i_dt_u_lon);
				Data::Vector::Vector<double> vel_lat_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(i_dt_v_lat);

				// Polar => Cartesian velocities
				Data::Vector::Vector<double> vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Setup iterations
				 */
				// Departure points for iterations
				Data::Vector::Vector<double> pos_x_D(sphere2DDataConfig->grid_number_elements);
				Data::Vector::Vector<double> pos_y_D(sphere2DDataConfig->grid_number_elements);
				Data::Vector::Vector<double> pos_z_D(sphere2DDataConfig->grid_number_elements);

				doAdvectionOnSphere2D(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-0.5*vel_x_A,
					-0.5*vel_y_A,
					-0.5*vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					semi_lagrangian_approximate_sphere2d_geometry
				);


				/*
				 * 2 iterations to get midpoint
				 */
				double diff = -1;
				for (int i = 0; i < semi_lagrangian_max_iterations; i++)
				{
					sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(
							pos_x_D, pos_y_D, pos_z_D,
							o_pos_lon_D, o_pos_lat_D
						);

					Data::Vector::Vector<double> u_D = sphere2DSampler.bilinear_scalar(i_dt_u_lon, o_pos_lon_D, o_pos_lat_D, true, shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points);
					Data::Vector::Vector<double> v_D = sphere2DSampler.bilinear_scalar(i_dt_v_lat, o_pos_lon_D, o_pos_lat_D, true, shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					Data::Vector::Vector<double> vel_x_D(num_elements), vel_y_D(num_elements), vel_z_D(num_elements);
					sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
							o_pos_lon_D, o_pos_lat_D,
							u_D, v_D,
							vel_x_D, vel_y_D, vel_z_D
						);

					// convert final points from Cartesian space to angular space
					Data::Vector::Vector<double> new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					doAdvectionOnSphere2D(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-0.5*vel_x_D,
						-0.5*vel_y_D,
						-0.5*vel_z_D,

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						semi_lagrangian_approximate_sphere2d_geometry
					);


					if (semi_lagrangian_convergence_threshold > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
										(pos_y_D-new_pos_y_D).reduce_maxAbs() +
										(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < semi_lagrangian_convergence_threshold)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (semi_lagrangian_convergence_threshold > 0)
				{
					if (diff > semi_lagrangian_convergence_threshold)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << semi_lagrangian_convergence_threshold << std::endl;
					}
				}

				/*
				 * Given the midpoint at pos_?_D, we compute the full time step
				 */

				Data::Vector::Vector<double> dot2 = 2.0*(pos_x_D*pos_x_A + pos_y_D*pos_y_A + pos_z_D*pos_z_A);

				pos_x_D = dot2*pos_x_D - pos_x_A;
				pos_y_D = dot2*pos_y_D - pos_y_A;
				pos_z_D = dot2*pos_z_D - pos_z_A;


				sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);

			}
			else if (trajectory_method == E_TRAJECTORY_METHOD_SETTLS_HORTAL)
			{
				/**
				 * See SETTLS paper
				 * Hortal, M. (2002). The development and testing of a new two-time-level semi-Lagrangian scheme (SETTLS) in the ECMWF forecast model. Q. J. R. Meteorol. Soc., 2, 1671–1687.
				 *
				 * We use the SETTLS formulation also to compute the departure points.
				 */

				// Extrapolate velocities at departure points
				Data::Sphere2D::DataGrid u_extrapol = 2.0*i_dt_u_lon - i_u_lon_prev;
				Data::Sphere2D::DataGrid v_extrapol = 2.0*i_dt_v_lat - i_v_lat_prev;

				// Convert velocities along lon/lat to scalar data array
				Data::Vector::Vector<double> vel_lon_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(i_dt_u_lon);
				Data::Vector::Vector<double> vel_lat_array = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(i_dt_v_lat);

				// Polar => Cartesian velocities
				Data::Vector::Vector<double> vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Setup iterations
				 */

				// Departure points for iterations
				Data::Vector::Vector<double> pos_x_D(sphere2DDataConfig->grid_number_elements);
				Data::Vector::Vector<double> pos_y_D(sphere2DDataConfig->grid_number_elements);
				Data::Vector::Vector<double> pos_z_D(sphere2DDataConfig->grid_number_elements);

				doAdvectionOnSphere2D(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-vel_x_A,
					-vel_y_A,
					-vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					semi_lagrangian_approximate_sphere2d_geometry
				);

				double diff = -1;
				for (int iters = 0; iters < semi_lagrangian_max_iterations; iters++)
				{
					sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);

					/*
					 * WARNING: Never convert this to vort/div space!!!
					 * This creates some artificial waves
					 */
					Data::Vector::Vector<double> vel_lon_extrapol_D = sphere2DSampler.bilinear_scalar(
							u_extrapol,
							o_pos_lon_D, o_pos_lat_D,
							true,
							shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points
						);
					Data::Vector::Vector<double> vel_lat_extrapol_D = sphere2DSampler.bilinear_scalar(
							v_extrapol,
							o_pos_lon_D, o_pos_lat_D,
							true,
							shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points
						);

					// convert extrapolated velocities to Cartesian velocities
					Data::Vector::Vector<double> vel_x_extrapol_D(num_elements), vel_y_extrapol_D(num_elements), vel_z_extrapol_D(num_elements);

					// polar => Cartesian coordinates
					sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
							o_pos_lon_D, o_pos_lat_D,
							vel_lon_extrapol_D, vel_lat_extrapol_D,
							vel_x_extrapol_D, vel_y_extrapol_D, vel_z_extrapol_D
						);

					Data::Vector::Vector<double> new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					doAdvectionOnSphere2D(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-0.5*(vel_x_extrapol_D + vel_x_A),
						-0.5*(vel_y_extrapol_D + vel_y_A),
						-0.5*(vel_z_extrapol_D + vel_z_A),

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						semi_lagrangian_approximate_sphere2d_geometry
					);

					if (semi_lagrangian_convergence_threshold > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
								(pos_y_D-new_pos_y_D).reduce_maxAbs() +
								(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < semi_lagrangian_convergence_threshold)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (semi_lagrangian_convergence_threshold > 0)
				{
					if (diff > semi_lagrangian_convergence_threshold)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << semi_lagrangian_convergence_threshold << std::endl;
					}
				}

				// convert final points from Cartesian space to angular space
				sweet::LibMath::VectorMath::point_cartesian_2_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else
			{
				SWEETErrorFatal("Unknown departure point calculation method");
			}

			return;
		}

		SWEETErrorFatal("Only 1st and 2nd order time integration supported");
	}



	/*
	 * Interpolation of prognostic fields at departure points.
	 *
	 * We assume the velocity U-V to be the SL advected field!
	 */
	void apply_sl_timeintegration_vd(
			const Data::Sphere2D::Operators *i_ops,

			const Data::Sphere2D::DataSpectral &i_phi,
			const Data::Sphere2D::DataSpectral &i_vrt,
			const Data::Sphere2D::DataSpectral &i_div,

			const Data::Vector::Vector<double> &i_pos_lon_d,
			const Data::Vector::Vector<double> &i_pos_lat_d,

			Data::Sphere2D::DataSpectral &o_phi,
			Data::Sphere2D::DataSpectral &o_vrt,
			Data::Sphere2D::DataSpectral &o_div
	)
	{
		o_phi.setup_if_required(i_phi.sphere2DDataConfig);
		o_vrt.setup_if_required(i_phi.sphere2DDataConfig);
		o_div.setup_if_required(i_phi.sphere2DDataConfig);

		o_phi = sphere2DSampler.bicubic_scalar_ret_phys(
				i_phi.toGrid(),
				i_pos_lon_d, i_pos_lat_d,
				false,
				shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSLData->semi_lagrangian_interpolation_limiter
			);

		o_vrt = sphere2DSampler.bicubic_scalar_ret_phys(
				i_vrt.toGrid(),
				i_pos_lon_d, i_pos_lat_d,
				false,
				shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSLData->semi_lagrangian_interpolation_limiter
			);

		o_div = sphere2DSampler.bicubic_scalar_ret_phys(
				i_div.toGrid(),
				i_pos_lon_d, i_pos_lat_d,
				false,
				shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSLData->semi_lagrangian_interpolation_limiter
			);
	}


	/*
	 * Interpolation of prognostic fields at departure points.
	 *
	 * We assume the velocity U-V to be the SL advected field!
	 */
	void apply_sl_timeintegration_uv(
			const Data::Sphere2D::Operators *i_ops,

			const Data::Sphere2D::DataSpectral &i_phi,	// input field
			const Data::Sphere2D::DataSpectral &i_vrt,
			const Data::Sphere2D::DataSpectral &i_div,

			const Data::Vector::Vector<double> &i_pos_lon_D,	// position
			const Data::Vector::Vector<double> &i_pos_lat_D,

			Data::Sphere2D::DataSpectral &o_phi,	// output field
			Data::Sphere2D::DataSpectral &o_vrt,
			Data::Sphere2D::DataSpectral &o_div
	)
	{
		o_phi.setup_if_required(i_phi.sphere2DDataConfig);
		o_vrt.setup_if_required(i_phi.sphere2DDataConfig);
		o_div.setup_if_required(i_phi.sphere2DDataConfig);


		/*************************************************************************
		 * Phi
		 *************************************************************************
		 */
		o_phi = sphere2DSampler.bicubic_scalar_ret_phys(
				i_phi.toGrid(),
				i_pos_lon_D, i_pos_lat_D,
				false,
				shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSLData->semi_lagrangian_interpolation_limiter
			);


	#if 1

		/*************************************************************************
		 * Prepare rotation system for handling velocities
		 *************************************************************************
		 */

		Data::Vector::Vector<double> P_x_D, P_y_D, P_z_D;
		sweet::LibMath::VectorMath::point_latlon_2_cartesian__array(i_pos_lon_D, i_pos_lat_D, P_x_D, P_y_D, P_z_D);

		Data::Vector::Vector<double> &P_x_A = pos_x_A,
						&P_y_A = pos_y_A,
						&P_z_A = pos_z_A;

		/*
		 * Compute rotation angle
		 */

		Data::Vector::Vector<double> rotation_angle_ =
				sweet::LibMath::VectorMath::dot_prod(
					P_x_D, P_y_D, P_z_D,
					P_x_A, P_y_A, P_z_A
				);

		// Can be slightly larger than 1, leading to NaN, hence this hack
		rotation_angle_ = sweet::LibMath::VectorMath::min(rotation_angle_, 1.0);

		Data::Vector::Vector<double> rotation_angle = sweet::LibMath::VectorMath::arccos(rotation_angle_);


		/*
		 * Compute Rotation axis
		 */
		Data::Vector::Vector<double> rot_x, rot_y, rot_z;
		sweet::LibMath::VectorMath::cross_prod(
				P_x_D, P_y_D, P_z_D,
				P_x_A, P_y_A, P_z_A,
				rot_x, rot_y, rot_z
			);

		sweet::LibMath::VectorMath::normalize_with_threshold(rot_x, rot_y, rot_z);



		/*************************************************************************
		 * Velocity
		 *************************************************************************
		 */

		Data::Sphere2D::DataGrid u_tmp, v_tmp;
		i_ops->vrtdiv_2_uv(i_vrt, i_div, u_tmp, v_tmp);

		Data::Sphere2D::DataGrid u_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				u_tmp,
				i_pos_lon_D, i_pos_lat_D,
				true,
				shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSLData->semi_lagrangian_interpolation_limiter
			);

		Data::Sphere2D::DataGrid v_tmp_D = sphere2DSampler.bicubic_scalar_ret_phys(
				v_tmp,
				i_pos_lon_D, i_pos_lat_D,
				true,
				shackSLData->semi_lagrangian_sampler_use_pole_pseudo_points,
				shackSLData->semi_lagrangian_interpolation_limiter
			);

		/*
		 * Convert to Cartesian space
		 */

		Data::Vector::Vector<double> V_lon_D = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(u_tmp_D);
		Data::Vector::Vector<double> V_lat_D = Data::Sphere2D::Convert::DataGrid_2_Vector::convert(v_tmp_D);

		Data::Vector::Vector<double> V_x_D, V_y_D, V_z_D;
		sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
				i_pos_lon_D,
				i_pos_lat_D,
				V_lon_D,
				V_lat_D,
				V_x_D,
				V_y_D,
				V_z_D
			);

		/*
		 * Rotate to velocity vector
		 */
		Data::Vector::Vector<double> V_x_A, V_y_A, V_z_A;
		sweet::LibMath::VectorMath::vector_rotate_3d_normalized_rotation_axis__array(
				V_x_D, V_y_D, V_z_D,
				rotation_angle,
				rot_x, rot_y, rot_z,
				V_x_A, V_y_A, V_z_A
			);


		/*
		 * Return velocity in lat/lon space
		 */
		Data::Vector::Vector<double> V_lon_A, V_lat_A;
		sweet::LibMath::VectorMath::velocity_cartesian_2_latlon__array(
				pos_lon_A,
				pos_lat_A,
				V_x_A,
				V_y_A,
				V_z_A,
				V_lon_A, V_lat_A
		);

		i_ops->uv_2_vrtdiv(
				Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_lon_A, i_vrt.sphere2DDataConfig),
				Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(V_lat_A, i_vrt.sphere2DDataConfig),
				o_vrt, o_div
			);

	#else

		op.uv_2_vrtdiv(u_tmp_D, v_tmp_D, o_vrt, o_div);

	#endif
	}

};

}}

#endif
