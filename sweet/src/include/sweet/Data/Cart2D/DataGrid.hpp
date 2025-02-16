/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_CART2D_DATAGRID_HPP


#include <complex>
#include <string>
#include <iostream>
#include <functional>
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid_Kernels.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>


namespace sweet {
namespace Data {
namespace Cart2D {


class DataSpectral;

/*!
 * \brief Data storage class for real valued 2D periodic data on a Cartesian grid
 */
class DataGrid
	:	private DataGrid_Kernels
{

public:
	const Config *cart2DDataConfig;

public:
	double *grid_space_data;

	public:
	DataGrid();

DataGrid(
			const Config *i_cart2DDataConfig
	);

public:
	DataGrid(
			const Config &i_cart2DDataConfig
	);

public:
	DataGrid(
			const DataGrid &i_cart2d_data
	);

	DataGrid(
			DataGrid &&i_cart2d_data
	);


	DataGrid(
			const Config *i_cart2DDataConfig,
			double i_value
	);

public:
		/*!
	 * dummy initialization by handing over an unused integer
	 */
public:
	DataGrid(int i);

public:
	public:
		/*!
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	void swap(
			DataGrid &i_cart2DData
	);

void check(
			const Config *i_cart2DDataConfig
	) const;



public:
	DataGrid& operator=(
			const DataGrid &i_cart2d_data
	);


public:
	DataGrid& operator=(
			DataGrid &&i_cart2d_data
	);


public:
	/*!
	 * assignment operator
	 */
	DataGrid &operator=(double i_value);


public:
	/*!
	 * assignment operator
	 */
	DataGrid &operator=(int i_value);


public:
	void setup(
			const Config &i_cart2DDataConfig
	);

public:
	void setup(
			const Config *i_cart2DDataConfig
	);


private:
	void alloc_data();

public:
	void setup_if_required(
			const Config *i_cart2DDataConfig
	);

public:
	~DataGrid();

public:
	void clear();

public:
	template <int S>
	void kernel_stencil_setup(
			const double i_kernel_array[S][S],
			double i_scale = 1.0
	)
	{
		((DataGrid_Kernels&)*this).kernel_stencil_setup(
				i_kernel_array,
				i_scale,

				cart2DDataConfig,
				grid_space_data
		);
	}



public:
	void grid_update_lambda_array_indices(
			std::function<void(int,int,double&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool anti_aliasing = true
	);


	void grid_update_lambda_array_idx(
			std::function<void(int,double&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool anti_aliasing = true
	);


	void grid_update_lambda_unit_coordinates_corner_centered(
			std::function<void(double,double,double&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool anti_aliasing = true
	);

	void grid_update_lambda_unit_coordinates_cell_centered(
			std::function<void(double,double,double&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool anti_aliasing = true
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
	 * Set given point a specific value
	 */
	void grid_setValue(
			int i_y_idx,
			int i_x_idx,
			double i_value
	);

	double grid_get(
			int i_y_idx,
			int i_x_idx
	)	const;


	DataGrid operator+(
			const DataGrid &i_cart2d_data
	)	const;


	DataGrid& operator+=(
			const DataGrid &i_cart2d_data
	);


	DataGrid& operator+=(
			double i_scalar
	);


	DataGrid& operator-=(
			const DataGrid &i_cart2d_data
	);


	DataGrid& operator-=(
			double i_scalar
	);


	DataGrid operator-(
			const DataGrid &i_cart2d_data
	)	const;


	DataGrid operator-();



	DataGrid operator*(
			const DataGrid &i_cart2d_data
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

	DataGrid& operator/=(
			double i_scalar
	);

	DataGrid operator+(
			double i_value
	)	const;

	DataGrid operator-(
			double i_value
	)	const;

	DataGrid operator_scalar_sub_this(
			double i_value
	)	const;


	/*!
	 * Apply a linear operator given by this class to the input data array.
	 */
	DataGrid operator()(
			const DataGrid &i_array_data
	)	const;


	DataGrid operator/(
			const DataGrid &i_cart2d_data
	)	const;


	/*!
	 * Simply apply a dealiasing to physical data
	 */
	void dealiasing(
		DataGrid& io_data
	) const;


	/*!
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double grid_reduce_max(
			const DataGrid &i_cart2d_data
	);


	double grid_reduce_rms();


	double grid_reduce_sum()	const;


	/*!
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double grid_reduce_sum_quad()	const;

	/*!
	 * reduce to root mean square
	 */
	double grid_reduce_rms_quad();

	/*!
	 * return the sum of the absolute values
	 */
	double grid_reduce_norm1()	const;


	/*!
	 * return the sqrt of the sum of the squared values
	 */
	double grid_reduce_norm2()	const;


	/*!
	 * return the sqrt of the sum of the squared values, use quad precision for reduction
	 */
	double grid_reduce_norm2_quad()	const;




	/*!
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double grid_reduce_max_abs(
			const DataGrid &i_cart2d_data
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

	/*!
	 * return true, if any value is infinity
	 */
	bool grid_reduce_boolean_all_finite() const;

	void grid_print(
			int i_precision = -1
	)	const;

	/*!
	 * print spectral data and zero out values which are numerically close to zero
	 */
	void print_physicalData_zeroNumZero(double i_zero_threshold = 1e-13)	const;


	void grid_file_write(
			const std::string &i_filename,
			double cart2d_domain_size[2],
			const char *i_title = "",
			int i_precision = 20
	)	const;

	bool file_grid_saveData_ascii(
			const char *i_filename,		//!< Name of file to store data to
			char i_separator = '\t',	//!< separator to use for each line
			int i_precision = 16,		//!< number of floating point digits
			int dimension = 2			//!< store 1D or 2D
	)	const;


	/*!
	 * Load data from ASCII file.
	 *
	 * Read csv files from reference data in parareal, in order to compute parareal erros online
	 * instead of producing csv files during the parareal simulation.
	 *
	 * \return true if data was successfully read
	 */

	///__attribute__((deprecated))
	bool file_read_csv_grid(
			const char *i_filename		//!< Name of file to load data from
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

#include <sweet/Data/Cart2D/DataGrid_Operators_inc.hpp>

}}}

#include <sweet/Data/Cart2D/DataGrid_Operators_inc.hpp>

#endif
