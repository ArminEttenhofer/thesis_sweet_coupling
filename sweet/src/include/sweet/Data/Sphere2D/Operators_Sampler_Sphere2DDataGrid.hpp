/*
 * Sphere2DDataSampler2D.hpp
 *
 *  Created on: 29 Mar 2018
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_OPERATORS_SAMPLER_SPHERE2DDATAGRID_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_OPERATORS_SAMPLER_SPHERE2DDATAGRID_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/LibMath/Interpolation.hpp>



namespace sweet {
namespace Data {
namespace Sphere2D {



/*!
 * This is a sampler class which provides method to provide
 * interpolated sampled values on 2D physical sphere2D data which
 * is provided by Sphere2DDataGrid
 */
class Operators_Sampler_DataGrid
{
public:
	int res[2];						//! resolution of domain
	const Sphere2D::Config *sphere2DDataConfig;

	std::vector<double> sampling_data;

private:
//	double cached_scale_factor[2];			//! cached parameters for sampling

	// number of extended latitudinal points (num_lat + 4)
	int ext_lat_M;

	// lookup table with latitudinal angles extended by 2 at front and back
	std::vector<double> phi_lookup;

	// lookup table using pseudo points at poles
	std::vector<double> phi_lookup_pseudo_points;


public:
	Operators_Sampler_DataGrid(
		Config *i_sphere2DDataConfig
	)
	{
		SWEET_ASSERT(i_sphere2DDataConfig != nullptr);

		sphere2DDataConfig = i_sphere2DDataConfig;
		setup(sphere2DDataConfig);
	}


	Operators_Sampler_DataGrid()
	{
		sphere2DDataConfig = nullptr;

		res[0] = -1;
		res[1] = -1;
	}


	/*
	 * Handler to sampling data, but with assertions
	 */
	inline
	double sampling_data_(
			int j,		// row
			int i		// col
	)
	{
		int num_lon = sphere2DDataConfig->grid_num_lon;
		int num_lat = sphere2DDataConfig->grid_num_lat;
		SWEET_ASSERT((i >= 0 && i < num_lon));
		SWEET_ASSERT((j >= 0 && j < (num_lat+4)));

		return sampling_data[j*num_lon + i];
	}

public:
	void setup(
		const Sphere2D::Config *i_sphere2DDataConfig
	)
	{
		if (i_sphere2DDataConfig == nullptr)
			SWEETErrorFatal("Setup called twice for Sphere2DOperators_Sampler_Sphere2DDataGrid");

		sphere2DDataConfig = i_sphere2DDataConfig;

		res[0] = i_sphere2DDataConfig->grid_num_lon;
		res[1] = i_sphere2DDataConfig->grid_num_lat;


		/*
		 * Use extended lat/lon lookup table to avoid if brances
		 */
		ext_lat_M = sphere2DDataConfig->grid_num_lat+4;

		/*
		 * Regular phi lookup
		 */
		phi_lookup.resize(ext_lat_M);

		for (int i = 0; i < sphere2DDataConfig->grid_num_lat; i++)
			phi_lookup[i+2] = sphere2DDataConfig->lat[i];

		phi_lookup[0] = M_PI - phi_lookup[3];
		phi_lookup[1] = M_PI - phi_lookup[2];
		phi_lookup[ext_lat_M-1] = -M_PI - phi_lookup[ext_lat_M-4];
		phi_lookup[ext_lat_M-2] = -M_PI - phi_lookup[ext_lat_M-3];


		/*
		 * Phi lookup with pseudo points at the poles
		 */
		phi_lookup_pseudo_points.resize(ext_lat_M);

		for (int i = 0; i < sphere2DDataConfig->grid_num_lat; i++)
			phi_lookup_pseudo_points[i+2] = sphere2DDataConfig->lat[i];

		phi_lookup_pseudo_points[0] = M_PI - phi_lookup_pseudo_points[2];
		phi_lookup_pseudo_points[1] = M_PI*0.5;	// northpole
		phi_lookup_pseudo_points[ext_lat_M-1] = -M_PI - phi_lookup_pseudo_points[ext_lat_M-3];
		phi_lookup_pseudo_points[ext_lat_M-2] = -M_PI*0.5;	// southpole



		/*
		 * Test for proper cubic interpolation
		 */
		double a[4] = {1.0, 2.0, 3.0, 4.0};
		double x0 = 1.3;
		auto f = [&](double x) -> double
		{
			return a[0] + x*a[1] + x*x*a[2] + x*x*x*a[3];
		};

		double x[4] = {0.1, 0.2, 0.4, 0.8};
		double y[4] = {f(x[0]), f(x[1]), f(x[2]), f(x[3])};

		double retval = sweet::LibMath::interpolation_lagrange_nonequidistant<4>(x, y, x0);

		if (std::abs(retval - f(x0)) > 1e-10)
			SWEETErrorFatal("Cubic interpolation Buggy!!!");
	}



	void updateSamplingData(
			const Sphere2D::DataGrid &i_data,
			int i_order,
			bool i_velocity_sampling,// = false,
			bool i_pole_pseudo_points// = false
	)
	{
#if SWEET_DEBUG
		if (i_pole_pseudo_points && i_velocity_sampling)
			SWEETErrorFatal("Not supported yet (or it doesn't make sense)!");
#endif

		int num_lon = sphere2DDataConfig->grid_num_lon;
		int num_lat = sphere2DDataConfig->grid_num_lat;
		int num_lat_ext = num_lat+4;

		// resize to store 4 additional rows (2 for north and 2 for south pole)
		sampling_data.resize(num_lon*num_lat_ext);


		int num_lon_d2 = sphere2DDataConfig->grid_num_lon/2;

		SWEET_ASSERT((num_lon & 1) == 0);

		/*
		 * Copy data at the center first!
		 */
		for (int i = 0; i < num_lon*num_lat; i++)
			sampling_data[2*num_lon + i] = i_data.grid_space_data[i];

		if (!i_velocity_sampling)
		{
			// first block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[0*num_lon + i] = i_data.grid_space_data[1*num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[0*num_lon + num_lon_d2 + i] = i_data.grid_space_data[1*num_lon + i];

			// second block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[1*num_lon + i] = i_data.grid_space_data[0*num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[1*num_lon + num_lon_d2 + i] = i_data.grid_space_data[0*num_lon + i];

			if (i_pole_pseudo_points)
			{
				// compute average
				double a = 0;

				if (1)
				{
					for (int i = 0; i < num_lon; i++)
					{
						double q[4];
						for (int j = 0; j < 4; j++)
							q[j] = sampling_data[j*num_lon + i];

						SWEET_ASSERT(phi_lookup_pseudo_points[1] == M_PI/2);
						a += sweet::LibMath::interpolation_lagrange_nonequidistant<4>(&phi_lookup[0], q, phi_lookup_pseudo_points[1]);
					}
				}
				else
				{
					SWEETErrorFatal("This order is not supported");
				}

				// average all points
				// TODO: Avg over half of points should be enough
				a /= num_lon;

				// copy 2nd last row to last row
				for (int i = 0; i < num_lon; i++)
					sampling_data[0*num_lon + i] = sampling_data[1*num_lon + i];

				// fill in avg values at pole
				for (int i = 0; i < num_lon; i++)
					sampling_data[1*num_lon + i] = a;
			}
		}
		else
		{
			// first block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[i] = -i_data.grid_space_data[num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon_d2 + i] = -i_data.grid_space_data[num_lon + i];

			// second block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon + i] = -i_data.grid_space_data[num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon + num_lon_d2 + i] = -i_data.grid_space_data[i];
		}


		if (!i_velocity_sampling)
		{
			// last block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[(num_lat+2)*num_lon + i] = i_data.grid_space_data[(num_lat-1)*num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[(num_lat+2)*num_lon + num_lon_d2 + i] = i_data.grid_space_data[(num_lat-1)*num_lon + i];

			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[(num_lat+3)*num_lon + i] = i_data.grid_space_data[(num_lat-2)*num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[(num_lat+3)*num_lon + num_lon_d2 + i] = i_data.grid_space_data[(num_lat-2)*num_lon + i];

			if (i_pole_pseudo_points)
			{
				// compute average
				double a = 0;

				if (1)
				{
					for (int i = 0; i < num_lon; i++)
					{
						double q[4];
						for (int j = 0; j < 4; j++)
							q[j] = sampling_data[(num_lat_ext-4+j)*num_lon + i];

						SWEET_ASSERT(phi_lookup_pseudo_points[num_lat_ext-2] == -M_PI/2);
						a += sweet::LibMath::interpolation_lagrange_nonequidistant<4>(&phi_lookup[num_lat_ext-4], q, -M_PI/2);
					}
				}
				else
				{
					SWEETErrorFatal("This order is not supported");
				}

				// average all points
				a /= num_lon;

				// copy 2nd last row to last row
				for (int i = 0; i < num_lon; i++)
					sampling_data[(num_lat_ext-1)*num_lon + i] = sampling_data[(num_lat_ext-2)*num_lon + i];

				for (int i = 0; i < num_lon; i++)
					sampling_data[(num_lat_ext-2)*num_lon + i] = a;
			}
		}
		else
		{
			// last block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+2) + i] = -i_data.grid_space_data[num_lon*(num_lat-1) + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+2) + num_lon_d2 + i] = -i_data.grid_space_data[num_lon*(num_lat-1) + i];

			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+3) + i] = -i_data.grid_space_data[num_lon*(num_lat-2) + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+3) + num_lon_d2 + i] = -i_data.grid_space_data[num_lon*(num_lat-2) + i];
		}
	}



public:
	/**
	 * wrap the position i in a periodic domain of size i_res
	 */
	template <typename T>
	static
	inline
	T wrapPeriodic(T i, T i_res)
	{
#if 1

		i = fmod(i, i_res);

		if (i < 0)
			i += i_res;

		// if i = -1e-14, then there's the special case if i==i_res
		if (i >= i_res)
			i = 0;

#else

		while (i < 0)
			i += i_res;

		while (i >= i_res)
			i -= i_res;
#endif

		SWEET_ASSERT(i >= 0 && i < i_res);

		return i;
	}



public:
	void bicubic_scalar(
			const Sphere2D::DataGrid &i_data,		//!< sampling data

			const Vector::Vector<double> &i_pos_lon,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_lat,		//!< y positions of interpolation points

			double *o_data,							//!< output values
			bool i_velocity_sampling,
			bool i_pole_pseudo_points,				//!< reconstruct pole points
			bool i_limiter							//!< Use limiter for interpolation to avoid unphysical local extrema
	)
	{
		SWEET_ASSERT(res[0] > 0);
		SWEET_ASSERT(i_pos_lon.numberOfElements == i_pos_lat.numberOfElements);
		SWEET_ASSERT((sphere2DDataConfig->grid_num_lon & 1) == 0);

		updateSamplingData(i_data, 3, i_velocity_sampling, i_pole_pseudo_points);

		std::vector<double> &phi_lookup_ = (i_pole_pseudo_points ? phi_lookup_pseudo_points : phi_lookup);

		// longitude angle delta
		double dlon = (double)i_data.sphere2DDataConfig->grid_num_lon / (2.0*M_PI);

		double L = -(-M_PI*0.5 - M_PI/ext_lat_M*1.5);

		// total size of lat field (M_PI + extension)
		// divided by number of cells
		//double s = (M_PI+M_PI/ext_lat_M*3)/(double)(ext_lat_M-1);
		double inv_s = (double)(ext_lat_M-1)/(M_PI+M_PI/ext_lat_M*3);

		// iterate over all positions in parallel
#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (std::size_t pos_idx = 0; pos_idx < i_pos_lon.numberOfElements; pos_idx++)
		{
			// load scalars
			double pos_lat = i_pos_lat.data[pos_idx];
			double pos_lon = i_pos_lon.data[pos_idx];

			/*
			 * We first ensure that the lon position is within [-0.5*pi; 0.5*pi]
			 * Otherwise, we correct it
			 */
#if 1
			if (pos_lat > M_PI*0.5)
			{
				std::cout << "TODO: Validate this since this special case is not yet tested" << std::endl;
				pos_lat = M_PI-pos_lat;
				pos_lon += M_PI;
			}

			if (pos_lat < -M_PI*0.5)
			{
				std::cout << "TODO: Validate this since this special case is not yet tested" << std::endl;
				pos_lat = -M_PI-pos_lat;
				pos_lon += M_PI;
			}
#endif

			/*
			 * Compute position in array
			 */
			double pos_array_x = wrapPeriodic(pos_lon*dlon, (double)res[0]);

			// compute position relative in cell \in [0;1]
			double cell_rel_x = pos_array_x - std::floor(pos_array_x);
			SWEET_ASSERT(cell_rel_x >= 0);
			SWEET_ASSERT(cell_rel_x <= 1);

			// compute array index
			int array_idx_x = std::floor(pos_array_x);
			SWEET_ASSERT(array_idx_x >= 0);
			SWEET_ASSERT(array_idx_x < sphere2DDataConfig->grid_num_lon);

			/*
			 * Compute Y array position
			 */
			/*
			 * Estimate array index for latitude.
			 * This should be off by only one cell.
			 */
			int est_lat_idx = (L - pos_lat)*inv_s;

#if SWEET_DEBUG
			if (!(est_lat_idx >= 0))
			{
				std::cout << "est_lat_idx: " << est_lat_idx << std::endl;
				std::cout << "L: " << L << std::endl;
				std::cout << "phi: " << pos_lat << std::endl;
				std::cout << "inv_s: " << inv_s << std::endl;
				SWEETErrorFatal("est_lat_idx");
			}
#endif
			SWEET_ASSERT(est_lat_idx >= 1);
			SWEET_ASSERT(est_lat_idx < ext_lat_M-1);

			if (phi_lookup_[est_lat_idx] < pos_lat)
				est_lat_idx--;
			else if (phi_lookup_[est_lat_idx+1] > pos_lat)
				est_lat_idx++;

			int array_idx_y = est_lat_idx;
			SWEET_ASSERT(array_idx_y >= 1);
			SWEET_ASSERT(array_idx_y < ext_lat_M-1);


			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */

			// precompute x-position indices since they are reused 4 times
			int idx_i[4];
			{
				idx_i[0] = wrapPeriodic(array_idx_x-1, res[0]);
				idx_i[1] = wrapPeriodic(array_idx_x+0, res[0]);
				idx_i[2] = wrapPeriodic(array_idx_x+1, res[0]);
				idx_i[3] = wrapPeriodic(array_idx_x+2, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = array_idx_y-1;

			double q[4];
			for (int kj = 0; kj < 4; kj++)
			{
				double p[4];

				p[0] = sampling_data_(idx_j, idx_i[0]);
				p[1] = sampling_data_(idx_j, idx_i[1]);
				p[2] = sampling_data_(idx_j, idx_i[2]);
				p[3] = sampling_data_(idx_j, idx_i[3]);

				//
				// 0....1..x..2....3
				//
				SWEET_ASSERT(cell_rel_x+1.0 >= 1.0 && cell_rel_x+1.0 <= 2.0);
				q[kj] = sweet::LibMath::interpolation_lagrange_equidistant<4>(p, cell_rel_x+1.0);

				if (i_limiter)
				{
					double max = std::max(p[1], p[2]);
					double min = std::min(p[1], p[2]);

					q[kj] = std::min(q[kj], max);
					q[kj] = std::max(q[kj], min);
				}

				idx_j++;
			}

			SWEET_ASSERT((phi_lookup_[array_idx_y] >= pos_lat) && (phi_lookup_[array_idx_y+1] <= pos_lat));
			double value = sweet::LibMath::interpolation_lagrange_nonequidistant<4>(&phi_lookup_[array_idx_y-1], q, pos_lat);

			if (i_limiter)
			{
				double max = std::max(q[1], q[2]);
				double min = std::min(q[1], q[2]);

				value = std::min(value, max);
				value = std::max(value, min);
			}

			o_data[pos_idx] = value;
		}
	}



public:
	void bicubic_scalar_new(
			const Sphere2D::DataGrid &i_data,			//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			Sphere2D::DataGrid &o_data,		//!< output values
			bool i_velocity_sampling,
			bool i_pole_pseudo_points,			//!< Use limiter for interpolation to avoid unphysical local extrema
			bool i_limiter
	)
	{
		o_data.setup_if_required(i_data.sphere2DDataConfig);

		SWEET_ASSERT(i_data.sphere2DDataConfig->grid_number_elements == o_data.sphere2DDataConfig->grid_number_elements);
		SWEET_ASSERT(res[0] > 0);
		SWEET_ASSERT(i_pos_x.numberOfElements == i_pos_y.numberOfElements);
		SWEET_ASSERT(i_pos_x.numberOfElements == (std::size_t)o_data.sphere2DDataConfig->grid_number_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.grid_space_data,  i_velocity_sampling, i_pole_pseudo_points, i_limiter);
	}


public:
	Sphere2D::DataGrid bicubic_scalar_ret_phys(
			const Sphere2D::DataGrid &i_data,			//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			bool i_velocity_sampling,
			bool i_pole_pseudo_points,			//!< Use limiter for interpolation to avoid unphysical local extrema
			bool i_limiter
	)
	{
		Sphere2D::DataGrid o_data(i_data.sphere2DDataConfig);

		SWEET_ASSERT(i_data.sphere2DDataConfig->grid_number_elements == o_data.sphere2DDataConfig->grid_number_elements);
		SWEET_ASSERT(res[0] > 0);
		SWEET_ASSERT(i_pos_x.numberOfElements == i_pos_y.numberOfElements);
		SWEET_ASSERT(i_pos_x.numberOfElements == (std::size_t)o_data.sphere2DDataConfig->grid_number_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.grid_space_data,  i_velocity_sampling, i_pole_pseudo_points, i_limiter);
		return o_data;
	}



public:
	void bilinear_scalar(
			const Sphere2D::DataGrid &i_data,	//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			double *o_data,						//!< output values
			bool i_velocity_sampling,			//!< swap sign for velocities,
			bool i_pole_pseudo_points
	)
	{
		SWEET_ASSERT(res[0] > 0);
		SWEET_ASSERT(i_pos_x.numberOfElements == i_pos_y.numberOfElements);
		SWEET_ASSERT((sphere2DDataConfig->grid_num_lon & 1) == 0);

		// copy the data to an internal buffer including halo layers
		updateSamplingData(i_data, 1, i_velocity_sampling, i_pole_pseudo_points);

		//int num_lon = sphere2DDataConfig->grid_num_lon;
		int num_lat = sphere2DDataConfig->grid_num_lat;

		// longitude spacing
		double s_lon = (double)i_data.sphere2DDataConfig->grid_num_lon / (2.0*M_PI);

		double L = -(-M_PI*0.5 - M_PI/ext_lat_M*1.5);
		// total size of lat field (M_PI + extension)
		// divided by number of cells
		//double s = (M_PI+M_PI/ext_lat_M*3)/(double)(ext_lat_M-1);
		double inv_s = (double)(ext_lat_M-1)/(M_PI+M_PI/ext_lat_M*3);

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t pos_idx = 0; pos_idx < i_pos_x.numberOfElements; pos_idx++)
		{
			/*
			 * Compute X information
			 */
			double array_x = wrapPeriodic(i_pos_x.data[pos_idx]*s_lon, (double)res[0]);

			// compute position relative in cell \in [0;1]
			double cell_rel_x = array_x - std::floor(array_x);
			SWEET_ASSERT(cell_rel_x >= 0);
			SWEET_ASSERT(cell_rel_x <= 1);

			// compute array index
			int array_idx_x = std::floor(array_x);
			SWEET_ASSERT(array_idx_x >= 0);
			SWEET_ASSERT(array_idx_x < sphere2DDataConfig->grid_num_lon);

			/*
			 * Compute Y information
			 *
			 * This is done via a lookup into phi_lookup since these
			 * coordinates are not equidistantly spaced, but close to it
			 */
			// estimate array index for latitude
			double phi = i_pos_y.data[pos_idx];
			int est_lat_idx = (L - phi)*inv_s;
#if SWEET_DEBUG
			if (!(est_lat_idx >= 0))
			{
				std::cout << "pos_idx: " << pos_idx << std::endl;
				std::cout << "est_lat_idx: " << est_lat_idx << std::endl;
				std::cout << "L: " << L << std::endl;
				std::cout << "phi: " << phi << std::endl;
				std::cout << "inv_s: " << inv_s << std::endl;
				SWEETErrorFatal("est_lat_idx");
			}
#endif
			SWEET_ASSERT(est_lat_idx >= 0);
			SWEET_ASSERT(est_lat_idx < ext_lat_M-1);

			if (phi_lookup[est_lat_idx] < phi)
				est_lat_idx--;
			else if (phi_lookup[est_lat_idx+1] > phi)
				est_lat_idx++;

			int array_idx_y = est_lat_idx;
			SWEET_ASSERT(array_idx_y >= 0);
			SWEET_ASSERT(array_idx_y < ext_lat_M);

			SWEET_ASSERT(phi_lookup[array_idx_y] >= phi);
			SWEET_ASSERT(phi_lookup[array_idx_y+1] <= phi);


			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */

			// precompute x-position indices since they are reused 4 times
			int idx_i[2];
			{
				idx_i[0] = wrapPeriodic(array_idx_x, res[0]);
				idx_i[1] = wrapPeriodic(array_idx_x+1, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = array_idx_y;

			double q[2];
			for (int kj = 0; kj < 2; kj++)
			{
				SWEET_ASSERT(idx_j >= 0);
				SWEET_ASSERT(idx_j < num_lat+4);
				double p[2];

				p[0] = sampling_data_(idx_j, idx_i[0]);
				p[1] = sampling_data_(idx_j, idx_i[1]);

				SWEET_ASSERT(0.0 <= cell_rel_x && cell_rel_x <= 1.0);
				q[kj] = sweet::LibMath::interpolation_lagrange_equidistant<2>(p, cell_rel_x);

				idx_j++;
			}

			// interpolation in y direction
			double value = sweet::LibMath::interpolation_lagrange_nonequidistant<2>(&phi_lookup[array_idx_y], q, phi);

			o_data[pos_idx] = value;
		}
	}

public:
	void bilinear_scalar(
			const Sphere2D::DataGrid &i_data,	//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			Vector::Vector<double> &o_data,			//!< output values
			bool i_velocity_sampling,			//!< swap sign for velocities
			bool i_pole_pseudo_points
	)
	{
		o_data.setup_if_required(i_pos_x.numberOfElements);

		SWEET_ASSERT(i_pos_x.numberOfElements == i_pos_y.numberOfElements);
		SWEET_ASSERT(i_pos_x.numberOfElements == (std::size_t)o_data.numberOfElements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.data, i_velocity_sampling, i_pole_pseudo_points);
	}


public:
	void bilinear_scalar(
			const Sphere2D::DataGrid &i_data,	//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			Sphere2D::DataGrid &o_data,		//!< output values
			bool i_velocity_sampling,			//!< swap sign for velocities
			bool i_pole_pseudo_points
	)
	{
		o_data.setup_if_required(i_data.sphere2DDataConfig);

		SWEET_ASSERT(i_pos_x.numberOfElements == i_pos_y.numberOfElements);
		SWEET_ASSERT(i_pos_x.numberOfElements == (std::size_t)o_data.sphere2DDataConfig->grid_number_elements);
		SWEET_ASSERT(o_data.grid_space_data != nullptr);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.grid_space_data, i_velocity_sampling, i_pole_pseudo_points);
	}


public:
	const Vector::Vector<double> bilinear_scalar(
			const Sphere2D::DataGrid &i_data,	//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			bool i_velocity_sampling,			//!< swap sign for velocities in halo regions
			bool i_pole_pseudo_points
	)
	{
		Vector::Vector<double> out(i_pos_x.numberOfElements);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, out, i_velocity_sampling, i_pole_pseudo_points);
		return out;
	}

public:
	const Sphere2D::DataGrid bilinear_scalar_ret_phys(
			const Sphere2D::DataGrid &i_data,	//!< sampling data

			const Vector::Vector<double> &i_pos_x,		//!< x positions of interpolation points
			const Vector::Vector<double> &i_pos_y,		//!< y positions of interpolation points

			bool i_velocity_sampling,			//!< swap sign for velocities in halo regions
			bool i_pole_pseudo_points
	)
	{
		Sphere2D::DataGrid o_data(i_data.sphere2DDataConfig);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.grid_space_data, i_velocity_sampling, i_pole_pseudo_points);

		return o_data;
	}


};

}}}

#endif
