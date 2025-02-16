/*
 * sweet::LibMath::VectorMath.hpp
 *
 *  Created on: Mar 27, 2020
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Error/Assert.hpp>

#ifndef INCLUDE_SWEET_LIBMATH_VECTORMATH_HPP
#define INCLUDE_SWEET_LIBMATH_VECTORMATH_HPP

namespace sweet {
namespace LibMath {


class VectorMath
{
public:
	inline
	static
	void normalize_vec3(
		const double i_vec[3],
		double o_vec[3]
	)
	{
		double s = 1.0/std::sqrt(i_vec[0]*i_vec[0] + i_vec[1]*i_vec[1] + i_vec[2]*i_vec[2]);

		o_vec[0] = i_vec[0] * s;
		o_vec[1] = i_vec[1] * s;
		o_vec[2] = i_vec[2] * s;
	}


	inline
	static
	void normalize_vec3(
		double i_vec[3]
	)
	{
		double s = 1.0/std::sqrt(i_vec[0]*i_vec[0] + i_vec[1]*i_vec[1] + i_vec[2]*i_vec[2]);

		i_vec[0] *= s;
		i_vec[1] *= s;
		i_vec[2] *= s;
	}


	inline
	static
	void normalize_vec3(
		double &io_vec_x,
		double &io_vec_y,
		double &io_vec_z
	)
	{
		double s = 1.0/std::sqrt(io_vec_x*io_vec_x + io_vec_y*io_vec_y + io_vec_z*io_vec_z);

		io_vec_x *= s;
		io_vec_y *= s;
		io_vec_z *= s;
	}


	inline
	static
	void normalize_vec3(
		const double &i_vec_x,
		const double &i_vec_y,
		const double &i_vec_z,
		double &o_vec_x,
		double &o_vec_y,
		double &o_vec_z
	)
	{
		double s = 1.0/std::sqrt(i_vec_x*i_vec_x + i_vec_y*i_vec_y + i_vec_z*i_vec_z);

		o_vec_x = i_vec_x*s;
		o_vec_y = i_vec_y*s;
		o_vec_z = i_vec_z*s;
	}


	inline
	static
	void normalize(
			sweet::Data::Vector::Vector<double> &io_v0,
			sweet::Data::Vector::Vector<double> &io_v1,
			sweet::Data::Vector::Vector<double> &io_v2
	)
	{
		sweet::Data::Vector::Vector<double> norm = (io_v0*io_v0 + io_v1*io_v1 + io_v2*io_v2).inv_sqrt();

		io_v0 *= norm;
		io_v1 *= norm;
		io_v2 *= norm;

	}



	/*
	 * Normalize by using threshold
	 *
	 * If the L2^2 value is smaller than the threshold, the norm is set to 1 instead of infinity
	 */
	inline
	static
	void normalize_with_threshold(
			sweet::Data::Vector::Vector<double> &io_v0,
			sweet::Data::Vector::Vector<double> &io_v1,
			sweet::Data::Vector::Vector<double> &io_v2,
			double i_threshold = 1e-20	// Default threshold for unit sphere2D and assuming this to be for the rotation axis
	)
	{
		double threshold2 = i_threshold*i_threshold;

		sweet::Data::Vector::Vector<double> length2 = io_v0*io_v0 + io_v1*io_v1 + io_v2*io_v2;

		length2.update_lambda_array_indices(
				[&](int i, double &io_data)
				{
					if (io_data < threshold2)
						io_data = threshold2;
				}
		);
		sweet::Data::Vector::Vector<double> norm = length2.inv_sqrt();

		io_v0 *= norm;
		io_v1 *= norm;
		io_v2 *= norm;

	}



	inline
	static
	sweet::Data::Vector::Vector<double> length(
			const sweet::Data::Vector::Vector<double> &i_v0,
			const sweet::Data::Vector::Vector<double> &i_v1,
			const sweet::Data::Vector::Vector<double> &i_v2
	)
	{
		return (i_v0*i_v0 + i_v1*i_v1 + i_v2*i_v2).sqrt();
	}



	static
	void cross_prod_v3(double w[3], double v[3], double ret[3])
	{
		ret[0] = w[1]*v[2] - w[2]*v[1];
		ret[1] = w[2]*v[0] - w[0]*v[2];
		ret[2] = w[0]*v[1] - w[1]*v[0];
	}


	static
	void cross_prod(
			const sweet::Data::Vector::Vector<double> &w0,
			const sweet::Data::Vector::Vector<double> &w1,
			const sweet::Data::Vector::Vector<double> &w2,

			const sweet::Data::Vector::Vector<double> &v0,
			const sweet::Data::Vector::Vector<double> &v1,
			const sweet::Data::Vector::Vector<double> &v2,

			sweet::Data::Vector::Vector<double> &o_ret0,
			sweet::Data::Vector::Vector<double> &o_ret1,
			sweet::Data::Vector::Vector<double> &o_ret2
	)
	{
		o_ret0.setup_if_required(w0);
		o_ret1.setup_if_required(w0);
		o_ret2.setup_if_required(w0);


		o_ret0 = w1*v2 - w2*v1;
		o_ret1 = w2*v0 - w0*v2;
		o_ret2 = w0*v1 - w1*v0;
	}


	static
	sweet::Data::Vector::Vector<double> dot_prod(
			const sweet::Data::Vector::Vector<double> &w0,
			const sweet::Data::Vector::Vector<double> &w1,
			const sweet::Data::Vector::Vector<double> &w2,

			const sweet::Data::Vector::Vector<double> &v0,
			const sweet::Data::Vector::Vector<double> &v1,
			const sweet::Data::Vector::Vector<double> &v2
	)
	{
		return w0*v0 + w1*v1 + w2*v2;
	}



	/*
	 * Make sure that lat/lon is within boundaries
	 */
	inline
	static
	void point_latlon_normalize__scalar(
			double &io_lon,
			double &io_lat
	)
	{
		if (io_lat > M_PI*0.5)
		{
			//io_lat = M_PI*0.5 - (io_lat-M_PI*0.5);
			io_lat = M_PI - io_lat;
			io_lon += M_PI;
			// longitude will be fixed if required below
		}
		else if (io_lat < -M_PI*0.5)
		{
			//io_lat = -M_PI*0.5 - (io_lat+M_PI*0.5);
			io_lat = -M_PI - io_lat;
			io_lon += M_PI;
			// longitude will be fixed if required below
		}

		if (io_lon < 0)
			io_lon += M_PI*2.0;
		else if (io_lon >= M_PI*2.0)
			io_lon -= M_PI*2.0;
	}



	/*
	 * Make sure that lat/lon is within boundaries
	 */
	inline
	static
	void point_latlon_normalize__array(
			const sweet::Data::Vector::Vector<double> &io_lon,
			const sweet::Data::Vector::Vector<double> &io_lat
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < io_lon.numberOfElements; i++)
		{
			point_latlon_normalize__scalar(
					io_lon.data[i],
					io_lat.data[i]
				);
		}
	}


	inline
	static
	void point_latlon_2_cartesian__scalar(
			const double i_lon,	//!< \in [0; 2pi]
			const double i_lat,	//!< \in [-pi/2; pi/2]
			double &o_ret_x,
			double &o_ret_y,
			double &o_ret_z
	)
	{
		SWEET_ASSERT(0 <= i_lon && i_lon <= 2.0*M_PI);
		SWEET_ASSERT(-0.5*M_PI <= i_lat && i_lat <= 0.5*M_PI);

		o_ret_x = std::cos(i_lon)*std::cos(i_lat);
		o_ret_y = std::sin(i_lon)*std::cos(i_lat);
		o_ret_z =                 std::sin(i_lat);
	}


	inline
	static
	void point_latlon_2_cartesian__array(
			const sweet::Data::Vector::Vector<double> &i_lon,	//!< \in [0; 2pi]
			const sweet::Data::Vector::Vector<double> &i_lat,	//!< \in [-pi/2; pi/2]
			sweet::Data::Vector::Vector<double> &o_x,
			sweet::Data::Vector::Vector<double> &o_y,
			sweet::Data::Vector::Vector<double> &o_z
	)
	{
		o_x.setup_if_required(i_lon);
		o_y.setup_if_required(i_lon);
		o_z.setup_if_required(i_lon);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.numberOfElements; i++)
		{
			point_latlon_2_cartesian__scalar(
					i_lon.data[i],
					i_lat.data[i],
					o_x.data[i],
					o_y.data[i],
					o_z.data[i]
				);
		}
	}




	inline
	static
	void point_cartesian_2_latlon__scalar(
			const double i_x,
			const double i_y,
			const double i_z,
			double &o_lon,
			double &o_lat
	)
	{
		SWEET_ASSERT(-1.0 <= i_x && i_x <= 1.0);
		SWEET_ASSERT(-1.0 <= i_y && i_y <= 1.0);
		SWEET_ASSERT(-1.0 <= i_z && i_z <= 1.0);

		/*
		 * Make sure that coordinates are in valid range
		 */
#if 0
		i_x = std::min(1., i_x);
		i_x = std::max(-1., i_x);

		i_y = std::min(1., i_y);
		i_y = std::max(-1., i_y);

		i_z = std::min(1., i_z);
		i_z = std::max(-1., i_z);
#endif

#if 0
		o_lon = std::atan(i_y/i_x);

		if (i_x < 0)
			o_lon += M_PI;
		else if (i_y < 0)
			o_lon += M_PI*2.0;

#else
		/*
		 * Now compute the angles using atan2
		 */
		o_lon = std::atan2(i_y, i_x);

		// Make sure that the angles are within [0;2pi]
		if (o_lon < 0)
			o_lon += 2.0*M_PI;
#endif

		o_lat = std::acos(-i_z) - M_PI*0.5;

		SWEET_ASSERT(0 <= o_lon && o_lon <= 2.0*M_PI);
		SWEET_ASSERT(-0.5*M_PI <= o_lat && o_lat <= 0.5*M_PI);

	}



	static
	void point_cartesian_2_latlon__array(
			const sweet::Data::Vector::Vector<double> &i_x,
			const sweet::Data::Vector::Vector<double> &i_y,
			const sweet::Data::Vector::Vector<double> &i_z,
			sweet::Data::Vector::Vector<double> &o_lon,
			sweet::Data::Vector::Vector<double> &o_lat
	)
	{
		o_lon.setup_if_required(i_x);
		o_lat.setup_if_required(i_x);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_x.numberOfElements; i++)
		{
			point_cartesian_2_latlon__scalar(
					i_x.data[i],
					i_y.data[i],
					i_z.data[i],
					o_lon.data[i],
					o_lat.data[i]
				);
		}
	}






	inline
	static
	void velocity_cartesian_2_latlon__scalar(
			const double i_lon,
			const double i_lat,
			const double i_v_x,
			const double i_v_y,
			const double i_v_z,
			double &o_vel_lon,
			double &o_vel_lat
	)
	{
		SWEET_ASSERT(0 <= i_lon && i_lon <= 2.0*M_PI);
		SWEET_ASSERT(-0.5*M_PI <= i_lat && i_lat <= 0.5*M_PI);

		o_vel_lon =	- std::sin(i_lon)*i_v_x
					+ std::cos(i_lon)*i_v_y;

		o_vel_lat =	- std::cos(i_lon)*std::sin(i_lat)*i_v_x
					- std::sin(i_lon)*std::sin(i_lat)*i_v_y
					+ std::cos(i_lat)*i_v_z;
	}




	inline
	static
	void velocity_cartesian_2_latlon__array(
			const sweet::Data::Vector::Vector<double> &i_lon,
			const sweet::Data::Vector::Vector<double> &i_lat,
			const sweet::Data::Vector::Vector<double> &i_v_x,
			const sweet::Data::Vector::Vector<double> &i_v_y,
			const sweet::Data::Vector::Vector<double> &i_v_z,
			sweet::Data::Vector::Vector<double> &o_vel_lon,
			sweet::Data::Vector::Vector<double> &o_vel_lat
	)
	{
		o_vel_lon.setup_if_required(i_lon);
		o_vel_lat.setup_if_required(i_lon);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.numberOfElements; i++)
		{
			velocity_cartesian_2_latlon__scalar(
					i_lon[i],
					i_lat[i],
					i_v_x[i],
					i_v_y[i],
					i_v_z[i],
					o_vel_lon[i],
					o_vel_lat[i]
				);
		}
	}




	/*
	 * Convert velocity in u-v (lat/lon) space to Cartesian space
	 */
	inline
	static
	void velocity_latlon_2_cartesian__scalar(
			const double i_lon,
			const double i_lat,
			const double i_vel_lon,
			const double i_vel_lat,
			double &o_v_x,
			double &o_v_y,
			double &o_v_z
	)
	{
			o_v_x = -i_vel_lon*std::sin(i_lon) - i_vel_lat*std::cos(i_lon)*std::sin(i_lat);
			o_v_y = i_vel_lon*std::cos(i_lon) - i_vel_lat*std::sin(i_lon)*std::sin(i_lat);
			o_v_z = i_vel_lat*std::cos(i_lat);
	}



	inline
	static
	void velocity_latlon_2_cartesian__array(
			const sweet::Data::Vector::Vector<double> &i_lon,
			const sweet::Data::Vector::Vector<double> &i_lat,
			const sweet::Data::Vector::Vector<double> &i_vel_lon,
			const sweet::Data::Vector::Vector<double> &i_vel_lat,
			sweet::Data::Vector::Vector<double> &o_v_x,
			sweet::Data::Vector::Vector<double> &o_v_y,
			sweet::Data::Vector::Vector<double> &o_v_z
	)
	{
		o_v_x.setup_if_required(i_lon);
		o_v_y.setup_if_required(i_lon);
		o_v_z.setup_if_required(i_lon);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.numberOfElements; i++)
		{
			velocity_latlon_2_cartesian__scalar(
					i_lon[i],
					i_lat[i],
					i_vel_lon[i],
					i_vel_lat[i],
					o_v_x[i],
					o_v_y[i],
					o_v_z[i]
				);
		}
	}



	static
	sweet::Data::Vector::Vector<double> cos(const sweet::Data::Vector::Vector<double> &i_values)
	{
		sweet::Data::Vector::Vector<double> retval = i_values;
		retval.update_lambda_array_indices(
			[](int, double &io_data)
			{
				io_data = std::cos(io_data);
			}
		);

		return retval;
	}


	static
	sweet::Data::Vector::Vector<double> min(const sweet::Data::Vector::Vector<double> &i_values, double i_scalar)
	{
		sweet::Data::Vector::Vector<double> retval = i_values;
		retval.update_lambda_array_indices(
			[&](int, double &io_data)
			{
				io_data = std::min(io_data, i_scalar);
			}
		);

		return retval;
	}



	static
	sweet::Data::Vector::Vector<double> max(const sweet::Data::Vector::Vector<double> &i_values, double i_scalar)
	{
		sweet::Data::Vector::Vector<double> retval = i_values;
		retval.update_lambda_array_indices(
			[&](int, double &io_data)
			{
				io_data = std::max(io_data, i_scalar);
			}
		);

		return retval;
	}



	static
	sweet::Data::Vector::Vector<double> arccos(const sweet::Data::Vector::Vector<double> &i_values)
	{
		sweet::Data::Vector::Vector<double> retval = i_values;
		retval.update_lambda_array_indices(
			[](int, double &io_data)
			{
				io_data = std::acos(io_data);
			}
		);

		return retval;
	}



	static
	sweet::Data::Vector::Vector<double> sin(const sweet::Data::Vector::Vector<double> &i_values)
	{
		sweet::Data::Vector::Vector<double> retval = i_values;
		retval.update_lambda_array_indices(
			[](int, double &io_data)
			{
				io_data = std::sin(io_data);
			}
		);

		return retval;
	}



	/**
	 * Rotate a point around an axis for the given rotation angle.
	 *
	 * The rotation axis is assumed to be normalized
	 */
	static
	void point_rotate_3d_normalized_rotation_axis__scalar(
		const double i_pos_start_x,
		const double i_pos_start_y,
		const double i_pos_start_z,
		const double i_rotation_angle,
		const double i_rotation_axis_x,
		const double i_rotation_axis_y,
		const double i_rotation_axis_z,
		double &o_pos_new_x,
		double &o_pos_new_y,
		double &o_pos_new_z
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */
	   double costheta = std::cos(i_rotation_angle);
	   double sintheta = std::sin(i_rotation_angle);

	   o_pos_new_x = (costheta + (1 - costheta) * i_rotation_axis_x * i_rotation_axis_x) * i_pos_start_x;
	   o_pos_new_x += ((1 - costheta) * i_rotation_axis_x * i_rotation_axis_y - i_rotation_axis_z * sintheta) * i_pos_start_y;
	   o_pos_new_x += ((1 - costheta) * i_rotation_axis_x * i_rotation_axis_z + i_rotation_axis_y * sintheta) * i_pos_start_z;

	   o_pos_new_y = ((1 - costheta) * i_rotation_axis_x * i_rotation_axis_y + i_rotation_axis_z * sintheta) * i_pos_start_x;
	   o_pos_new_y += (costheta + (1 - costheta) * i_rotation_axis_y * i_rotation_axis_y) * i_pos_start_y;
	   o_pos_new_y += ((1 - costheta) * i_rotation_axis_y * i_rotation_axis_z - i_rotation_axis_x * sintheta) * i_pos_start_z;

	   o_pos_new_z = ((1 - costheta) * i_rotation_axis_x * i_rotation_axis_z - i_rotation_axis_y * sintheta) * i_pos_start_x;
	   o_pos_new_z += ((1 - costheta) * i_rotation_axis_y * i_rotation_axis_z + i_rotation_axis_x * sintheta) * i_pos_start_y;
	   o_pos_new_z += (costheta + (1 - costheta) * i_rotation_axis_z * i_rotation_axis_z) * i_pos_start_z;
	}



	/**
	 * Rotate a point around an axis for the given rotation angle.
	 *
	 * The rotation axis is *not* assumed to be normalized!
	 * This is done as part of this function call.
	 */
	static
	void point_rotate_3d__scalar(
		const double i_pos_start_x,
		const double i_pos_start_y,
		const double i_pos_start_z,
		const double i_rotation_angle,
		const double &i_rotation_axis_x,
		const double &i_rotation_axis_y,
		const double &i_rotation_axis_z,
		double &o_pos_new_x,
		double &o_pos_new_y,
		double &o_pos_new_z
	)
	{
		double rotation_axis_x, rotation_axis_y, rotation_axis_z;
		normalize_vec3(
				i_rotation_axis_x, i_rotation_axis_y, i_rotation_axis_z,
				rotation_axis_x, rotation_axis_y, rotation_axis_z
			);

		point_rotate_3d_normalized_rotation_axis__scalar(
			i_pos_start_x,
			i_pos_start_y,
			i_pos_start_z,
			i_rotation_angle,
			rotation_axis_x,
			rotation_axis_y,
			rotation_axis_z,
			o_pos_new_x,
			o_pos_new_y,
			o_pos_new_z
		);
	}



	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void point_rotate_3d__array(
			const sweet::Data::Vector::Vector<double> &i_pos_start_0,
			const sweet::Data::Vector::Vector<double> &i_pos_start_1,
			const sweet::Data::Vector::Vector<double> &i_pos_start_2,
			const sweet::Data::Vector::Vector<double> &i_rotation_angle,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_0,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_1,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_2,
			sweet::Data::Vector::Vector<double> &o_pos_new_0,
			sweet::Data::Vector::Vector<double> &o_pos_new_1,
			sweet::Data::Vector::Vector<double> &o_pos_new_2
	)
	{
		o_pos_new_0.setup_if_required(i_pos_start_0);
		o_pos_new_1.setup_if_required(i_pos_start_0);
		o_pos_new_2.setup_if_required(i_pos_start_0);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_pos_start_0.numberOfElements; i++)
		{
			point_rotate_3d__scalar(
					i_pos_start_0[i],
					i_pos_start_1[i],
					i_pos_start_2[i],
					i_rotation_angle[i],
					i_rotation_axis_0[i],
					i_rotation_axis_1[i],
					i_rotation_axis_2[i],
					o_pos_new_0[i],
					o_pos_new_1[i],
					o_pos_new_2[i]
				);
		}
	}



	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void point_rotate_3d_normalized_rotation_axis__array(
			const sweet::Data::Vector::Vector<double> &i_pos_start_0,
			const sweet::Data::Vector::Vector<double> &i_pos_start_1,
			const sweet::Data::Vector::Vector<double> &i_pos_start_2,
			const sweet::Data::Vector::Vector<double> &i_rotation_angle,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_0,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_1,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_2,
			sweet::Data::Vector::Vector<double> &o_pos_new_0,
			sweet::Data::Vector::Vector<double> &o_pos_new_1,
			sweet::Data::Vector::Vector<double> &o_pos_new_2
	)
	{
		o_pos_new_0.setup_if_required(i_pos_start_0);
		o_pos_new_1.setup_if_required(i_pos_start_0);
		o_pos_new_2.setup_if_required(i_pos_start_0);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_pos_start_0.numberOfElements; i++)
		{
			point_rotate_3d_normalized_rotation_axis__scalar(
					i_pos_start_0[i],
					i_pos_start_1[i],
					i_pos_start_2[i],
					i_rotation_angle[i],
					i_rotation_axis_0[i],
					i_rotation_axis_1[i],
					i_rotation_axis_2[i],
					o_pos_new_0[i],
					o_pos_new_1[i],
					o_pos_new_2[i]
				);
		}
	}



#if 0
	/**
	 * Rotate a point around an axis for the given rotation angle.
	 *
	 * The rotation axis is assumed to be already normalized
	 */
	static
	void rotate_3d_point_normalized_rotation_axis(
		double i_pos_start[3],
		double i_rotation_angle,
		double i_rotation_axis_normalized[3],
		double i_pos_new[3]
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */

	   double cos_theta = std::cos(i_rotation_angle);
	   double sin_theta = std::sin(i_rotation_angle);

	   i_pos_new[0] = (cos_theta + (1 - cos_theta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[0]) * i_pos_start[0];
	   i_pos_new[0] += ((1 - cos_theta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[1] - i_rotation_axis_normalized[2] * sin_theta) * i_pos_start[1];
	   i_pos_new[0] += ((1 - cos_theta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[2] + i_rotation_axis_normalized[1] * sin_theta) * i_pos_start[2];

	   i_pos_new[1] = ((1 - cos_theta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[1] + i_rotation_axis_normalized[2] * sin_theta) * i_pos_start[0];
	   i_pos_new[1] += (cos_theta + (1 - cos_theta) * i_rotation_axis_normalized[1] * i_rotation_axis_normalized[1]) * i_pos_start[1];
	   i_pos_new[1] += ((1 - cos_theta) * i_rotation_axis_normalized[1] * i_rotation_axis_normalized[2] - i_rotation_axis_normalized[0] * sin_theta) * i_pos_start[2];

	   i_pos_new[2] = ((1 - cos_theta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[2] - i_rotation_axis_normalized[1] * sin_theta) * i_pos_start[0];
	   i_pos_new[2] += ((1 - cos_theta) * i_rotation_axis_normalized[1] * i_rotation_axis_normalized[2] + i_rotation_axis_normalized[0] * sin_theta) * i_pos_start[1];
	   i_pos_new[2] += (cos_theta + (1 - cos_theta) * i_rotation_axis_normalized[2] * i_rotation_axis_normalized[2]) * i_pos_start[2];
	}
#endif

#if 0
	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void rotate_3d_point_normalized_rotation_axis(
			const sweet::Data::Vector::Vector<double> &i_pos_start_0,
			const sweet::Data::Vector::Vector<double> &i_pos_start_1,
			const sweet::Data::Vector::Vector<double> &i_pos_start_2,
			const sweet::Data::Vector::Vector<double> &i_rotation_angle,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_0,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_1,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_2,
			sweet::Data::Vector::Vector<double> &o_pos_new_0,
			sweet::Data::Vector::Vector<double> &o_pos_new_1,
			sweet::Data::Vector::Vector<double> &o_pos_new_2
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */

		sweet::Data::Vector::Vector<double> cos_theta = cos(i_rotation_angle);
		sweet::Data::Vector::Vector<double> sin_theta = sin(i_rotation_angle);

		o_pos_new_0 = (cos_theta + (1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_0) * i_pos_start_0;
		o_pos_new_0 += ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_1 - i_rotation_axis_2 * sin_theta) * i_pos_start_1;
		o_pos_new_0 += ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_2 + i_rotation_axis_1 * sin_theta) * i_pos_start_2;

		o_pos_new_1 = ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_1 + i_rotation_axis_2 * sin_theta) * i_pos_start_0;
		o_pos_new_1 += (cos_theta + (1.0 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_1) * i_pos_start_1;
		o_pos_new_1 += ((1.0 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_2 - i_rotation_axis_0 * sin_theta) * i_pos_start_2;

		o_pos_new_2 = ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_2 - i_rotation_axis_1 * sin_theta) * i_pos_start_0;
		o_pos_new_2 += ((1.0 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_2 + i_rotation_axis_0 * sin_theta) * i_pos_start_1;
		o_pos_new_2 += (cos_theta + (1.0 - cos_theta) * i_rotation_axis_2 * i_rotation_axis_2) * i_pos_start_2;
	}
#endif


#if 0
	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void vector_rotate_3d_normalized_rotation_axis__scalar(
			const double i_vec_start_0,
			const double i_vec_start_1,
			const double i_vec_start_2,
			const double i_rotation_angle,
			const double i_rotation_axis_0,
			const double i_rotation_axis_1,
			const double i_rotation_axis_2,
			double &o_vec_new_0,
			double &o_vec_new_1,
			double &o_vec_new_2
	)
	{
		double cos_theta = std::cos(i_rotation_angle);
		double sin_theta = std::sin(i_rotation_angle);

		/*
		 * Based on source code from website:
		 *
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 *
		 * Main modification: Transposed version added for vector rotation
		 */

		o_vec_new_0 = 	i_vec_start_0 * (cos_theta + (1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_0)
						+ i_vec_start_1 * ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_1 + i_rotation_axis_2 * sin_theta)
						+ i_vec_start_2 * ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_2 - i_rotation_axis_1 * sin_theta);

		o_vec_new_1 =	i_vec_start_0 * ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_1 - i_rotation_axis_2 * sin_theta)
						+ i_vec_start_1 * (cos_theta + (1 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_1)
						+ i_vec_start_2 * ((1.0 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_2 + i_rotation_axis_0 * sin_theta);

		o_vec_new_2 =	i_vec_start_0 * ((1.0 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_2 + i_rotation_axis_1 * sin_theta)
						+ i_vec_start_1 * ((1.0 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_2 - i_rotation_axis_0 * sin_theta)
						+ i_vec_start_2 * (cos_theta + (1.0 - cos_theta) * i_rotation_axis_2 * i_rotation_axis_2);
	}




	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void vector_rotate_3d__scalar(
			const double i_vec_start_0,
			const double i_vec_start_1,
			const double i_vec_start_2,
			const double i_rotation_angle,
			const double i_rotation_axis_0,
			const double i_rotation_axis_1,
			const double i_rotation_axis_2,
			double &o_vec_new_0,
			double &o_vec_new_1,
			double &o_vec_new_2
	)
	{
		double rotation_axis_0, rotation_axis_1, rotation_axis_2;

		normalize_vec3(
				i_rotation_axis_0,
				i_rotation_axis_1,
				i_rotation_axis_2,
				rotation_axis_0,
				rotation_axis_1,
				rotation_axis_2
			);

		vector_rotate_3d_normalized_rotation_axis__scalar(
				i_vec_start_0,
				i_vec_start_1,
				i_vec_start_2,
				i_rotation_angle,
				rotation_axis_0,
				rotation_axis_1,
				rotation_axis_2,
				o_vec_new_0,
				o_vec_new_1,
				o_vec_new_2
			);
	}




	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void vector_rotate_3d(
			const sweet::Data::Vector::Vector<double> &i_vec_start_0,
			const sweet::Data::Vector::Vector<double> &i_vec_start_1,
			const sweet::Data::Vector::Vector<double> &i_vec_start_2,
			const sweet::Data::Vector::Vector<double> &i_rotation_angle,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_0,	//!< WARNING: This axis needs to be normalized!!!
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_1,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_2,
			sweet::Data::Vector::Vector<double> &o_vec_new_0,
			sweet::Data::Vector::Vector<double> &o_vec_new_1,
			sweet::Data::Vector::Vector<double> &o_vec_new_2
	)
	{
		o_vec_new_0.setup_if_required(i_vec_start_0);
		o_vec_new_1.setup_if_required(i_vec_start_1);
		o_vec_new_2.setup_if_required(i_vec_start_2);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_vec_start_0.numberOfElements; i++)
		{
			vector_rotate_3d__scalar(
						i_vec_start_0[i],
						i_vec_start_1[i],
						i_vec_start_2[i],
						i_rotation_angle[i],
						i_rotation_axis_0[i],
						i_rotation_axis_1[i],
						i_rotation_axis_2[i],
						o_vec_new_0[i],
						o_vec_new_1[i],
						o_vec_new_2[i]
				);
		}
	}


	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void vector_rotate_3d_normalized_rotation_axis__array(
			const sweet::Data::Vector::Vector<double> &i_vec_start_0,
			const sweet::Data::Vector::Vector<double> &i_vec_start_1,
			const sweet::Data::Vector::Vector<double> &i_vec_start_2,
			const sweet::Data::Vector::Vector<double> &i_rotation_angle,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_0,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_1,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_2,
			sweet::Data::Vector::Vector<double> &o_vec_new_0,
			sweet::Data::Vector::Vector<double> &o_vec_new_1,
			sweet::Data::Vector::Vector<double> &o_vec_new_2
	)
	{
		o_vec_new_0.setup_if_required(i_vec_start_0);
		o_vec_new_1.setup_if_required(i_vec_start_1);
		o_vec_new_2.setup_if_required(i_vec_start_2);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_vec_start_0.numberOfElements; i++)
		{
			vector_rotate_3d_normalized_rotation_axis__scalar(
						i_vec_start_0[i],
						i_vec_start_1[i],
						i_vec_start_2[i],
						i_rotation_angle[i],
						i_rotation_axis_0[i],
						i_rotation_axis_1[i],
						i_rotation_axis_2[i],
						o_vec_new_0[i],
						o_vec_new_1[i],
						o_vec_new_2[i]
				);
		}
	}
#endif


	/*
	 * The transformation doesn't include any stretching, etc.
	 * Hence, we can directly use the rotation matrix also for vectors rather
	 * than using M^-1^T as it would be the case in e.g. OpenGL
	 */
	static
	void vector_rotate_3d_normalized_rotation_axis__array(
			const sweet::Data::Vector::Vector<double> &i_vec_start_0,
			const sweet::Data::Vector::Vector<double> &i_vec_start_1,
			const sweet::Data::Vector::Vector<double> &i_vec_start_2,
			const sweet::Data::Vector::Vector<double> &i_rotation_angle,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_0,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_1,
			const sweet::Data::Vector::Vector<double> &i_rotation_axis_2,
			sweet::Data::Vector::Vector<double> &o_vec_new_0,
			sweet::Data::Vector::Vector<double> &o_vec_new_1,
			sweet::Data::Vector::Vector<double> &o_vec_new_2
	)
	{
		point_rotate_3d_normalized_rotation_axis__array(
					i_vec_start_0,
					i_vec_start_1,
					i_vec_start_2,
					i_rotation_angle,
					i_rotation_axis_0,
					i_rotation_axis_1,
					i_rotation_axis_2,
					o_vec_new_0,
					o_vec_new_1,
					o_vec_new_2
		);
	}
};

}}

#endif
