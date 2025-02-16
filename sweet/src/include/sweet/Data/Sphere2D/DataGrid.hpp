/*
 * Sphere2DDataGrid.hpp
 *
 *  Created on: 01 Jan 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_DATAGRID_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <complex>
#include <functional>
#include <array>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <utility>
#include <functional>
#include <vector>
#include <iterator>

#include <cmath>
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Fatal.hpp>



namespace sweet {
namespace Data {
namespace Sphere2D {

class DataSpectral;

class DataGrid
{
public:
	const Config *sphere2DDataConfig;

public:
	double *grid_space_data;


	void swap(
			DataGrid &i_sphere2DData
	);

public:
	DataGrid(
			const Config &i_sphere2DDataConfig
	);

public:
	DataGrid(
			const Config *i_sphere2DDataConfig
	);


public:
	DataGrid(
			const Config *i_sphere2DDataConfig,
			double i_value
	);

public:
	DataGrid();

public:
	DataGrid(
			const DataGrid &i_sph_data
	);

public:
	DataGrid(
			DataGrid &&i_sph_data
	);


	/*!
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const Config *i_sphere2DDataConfig
	) const;


public:
	DataGrid& operator=(
			const DataGrid &i_sph_data
	);


public:
	DataGrid& operator=(
			DataGrid &&i_sph_data
	);


	DataGrid operator+(
			const DataGrid &i_sph_data
	) const;


	DataGrid& operator+=(
			const DataGrid &i_sph_data
	);


	DataGrid& operator+=(
			double i_scalar
	);


	DataGrid& operator-=(
			const DataGrid &i_sph_data
	);


	DataGrid& operator-=(
			double i_scalar
	);


	DataGrid operator-(
			const DataGrid &i_sph_data
	)	const;


	DataGrid operator-();

	DataGrid operator*(
			const DataGrid &i_sph_data
	)	const;


	DataGrid operator/(
			const DataGrid &i_sph_data
	)	const;


	DataGrid operator*(
			const double i_value
	)	const;


	const DataGrid& operator*=(
			const double i_value
	)	const;


	DataGrid operator/(
			double i_value
	)	const;


	DataGrid operator+(
			double i_value
	)	const;



	DataGrid operator-(
			double i_value
	)	const;


	DataGrid operator_scalar_sub_this(
			double i_value
	)	const;


public:
	bool setup(
		const Config *i_sphere2DDataConfig
	);


public:
	bool setup(
		const Config &i_sphere2DDataConfig
	);



private:
	bool alloc_data();


public:
	void setup_if_required(
			const Config *i_sphere2DDataConfig
	);

public:
	~DataGrid();

public:
	void clear();


	/*!
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
	void grid_update_lambda(
			std::function<void(double,double,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	void grid_update_lambda_array(
			std::function<void(int,int,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	void grid_update_lambda_array_idx(
			std::function<void(int,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude sin(phi) \in [-1;1])
	 */
	void grid_update_lambda_gaussian_grid(
			std::function<void(double,double,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	/*!
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters:
	 *   (longitude \in [0;2*pi], Cogaussian latitude cos(phi) \in [0;1])
	 */
	void grid_update_lambda_cogaussian_grid(
			std::function<void(double,double,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);

	void grid_update_lambda_sinphi_grid(
			std::function<void(double,double,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);

	void grid_update_lambda_cosphi_grid(
			std::function<void(double,double,double&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	/*!
	 * Set all values to zero
	 */
	void grid_setZero();


	/*!
	 * Set all values to a specific value
	 */
	void grid_setValue(
			double i_value
	);


	/*!
	 * Set all values to a specific value
	 */
	__attribute__((deprecated))
	void grid_setValue(
			int i_lon_idx,
			int i_lat_idx,
			double i_value
	);

    double grid_getValue(
            int i_lon_idx,
            int i_lat_idx
    );


	/*!
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double grid_reduce_max(
			const DataGrid &i_sph_data
	);

	double grid_reduce_rms();

	double grid_reduce_norm1();

	double grid_reduce_norm2();


	double grid_reduce_sum()	const;


	/*!
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double grid_reduce_sum_quad()	const;


	/*!
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double grid_reduce_max_abs(
			const DataGrid &i_sph_data
	)	const;

	/*!
	 * Return the maximum absolute value
	 */
	double grid_reduce_max_abs()	const;


	/*!
	 * Return the minimum value
	 */
	double grid_reduce_min()	const;

	/*!
	 * Return the minimum value
	 */
	double grid_reduce_max()	const;


	bool grid_isAnyNaNorInf()	const;



	void grid_print(
			int i_precision = -1
	)	const;


	void grid_file_write(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const;


	void grid_file_write_lon_pi_shifted(
			const char *i_filename,
			const std::string &i_title = "",
			int i_precision = 20
	)	const;

	bool file_read_csv_grid(
			const char *i_filename,		//!< Name of file to load data from
			bool i_binary_data = false	//!< load as binary data (disabled per default)
	);

	void file_write_raw(
			const std::string &i_filename
	)	const;


	void file_read_raw(
			const std::string &i_filename
	)	const;


	void print_debug(
			const char *name
	)	const;

	void print()	const;
};




}}}



/*!
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
sweet::Data::Sphere2D::DataGrid operator*(
		const double i_value,
		const sweet::Data::Sphere2D::DataGrid &i_array_data
)
{
	return i_array_data*i_value;
}



/*!
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 */
inline
static
sweet::Data::Sphere2D::DataGrid operator-(
		const double i_value,
		const sweet::Data::Sphere2D::DataGrid &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}



/*!
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 */
inline
static
sweet::Data::Sphere2D::DataGrid operator+(
		const double i_value,
		const sweet::Data::Sphere2D::DataGrid &i_array_data
)
{
	return i_array_data+i_value;
}

#endif
