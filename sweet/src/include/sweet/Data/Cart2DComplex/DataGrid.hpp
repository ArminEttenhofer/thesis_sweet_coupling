/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATAGRID_HPP

#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
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

#include <cmath>
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Fatal.hpp>


namespace sweet {
namespace Data {
namespace Cart2DComplex {


/*!
 * \brief Data storage class for complex-valued 2D periodic data on a Cartesian grid
 */
class DataGrid
{

public:
	const Cart2D::Config *cart2DDataConfig;

public:
	std::complex<double> *grid_space_data;


	void swap(
			DataGrid &i_cart2DData
	);

public:
	DataGrid(
			const Cart2D::Config *i_cart2DDataConfig
	);


public:
	DataGrid();


public:
	DataGrid(
			const DataGrid &i_cart2d_data
	);


public:
	DataGrid(
			const Cart2D::DataGrid &i_cart2d_data
	);


	/*
	 * load real and imaginary data from physical arrays
	 */
	void loadRealImag(
			const Cart2D::DataGrid &i_re,
			const Cart2D::DataGrid &i_im
	);


	/*!
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const Cart2D::Config *i_cart2DDataConfig
	)	const;



public:
	DataGrid& operator=(
			const DataGrid &i_cart2d_data
	);


public:
	DataGrid& operator=(
			DataGrid &&i_cart2d_data
	);


public:
	DataGrid& operator=(
			const Cart2D::DataGrid &i_cart2d_data
	);


public:
	void setup(
			const Cart2D::Config *i_cart2DDataConfig
	);


public:
	void setup_if_required(
			const Cart2D::Config *i_cart2DDataConfig
	);


public:
	~DataGrid();

public:
	DataGrid operator+(
			const DataGrid &i_cart2d_data
	)	const;


	DataGrid& operator+=(
			const DataGrid &i_cart2d_data
	);


	DataGrid& operator-=(
			const DataGrid &i_cart2d_data
	);



	DataGrid operator-(
			const DataGrid &i_cart2d_data
	)	const;


	DataGrid operator-();

	DataGrid operator*(
			const DataGrid &i_cart2d_data
	)	const;


	DataGrid operator/(
			const DataGrid &i_cart2d_data
	)	const;


	DataGrid operator*(
			const double i_value
	)	const;


	DataGrid operator*(
			const std::complex<double> &i_value
	)	const;



	const DataGrid& operator*=(
			const double i_value
	)	const;

	DataGrid operator/(
			double i_value
	)	const;



	DataGrid operator+(
			const std::complex<double> &i_value
	)	const;


	DataGrid operator-(
			const std::complex<double> &i_value
	)	const;


	void grid_update_lambda_array_indices(
			std::function<void(int,int,std::complex<double>&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool i_anti_aliasing = true
	);

	void grid_update_lambda_array_idx(
			std::function<void(int,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	void grid_update_lambda_unit_coordinates_corner_centered(
			std::function<void(double,double,std::complex<double>&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool i_anti_aliasing = true
	);

	void grid_update_lambda_unit_coordinates_cell_centered(
			std::function<void(double,double,std::complex<double>&)> i_lambda,	//!< lambda function to return value for lat/mu
			bool i_anti_aliasing = true
	);

	/*
	 * Set all values to zero
	 */
	void grid_setZero();


	/*
	 * Set all values to a specific value
	 */
	void grid_setValue(
			std::complex<double> &i_value
	);

	/*
	 * Set all values to a specific value
	 */
	void grid_setValue(
			double i_value_real,
			double i_value_imag
	);


	/*
	 * Set all values to a specific value
	 */
	void grid_setValue(
			int i_y_idx,
			int i_x_idx,
			std::complex<double> &i_value
	);


	/*
	 * Set all values to a specific value
	 */
	void grid_setValue(
			int i_y_idx,
			int i_x_idx,
			double i_value_real,
			double i_value_imag
	);

	/*!
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double grid_reduce_max(
			const DataGrid &i_cart2d_data
	);

	double grid_reduce_rms();


	/*!
	 * Return the maximum error norm
	 */
	double grid_reduce_max_abs()	const;

	/*!
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double grid_reduce_sum_re_quad()	const;

	/*!
	 * reduce to root mean square
	 */
	double grid_reduce_rms_quad();

	/*!
	 * return the sqrt of the sum of the squared values, use quad precision for reduction
	 */
	double grid_reduce_norm2_quad()	const;
};

}}}

#include <sweet/Data/Cart2DComplex/DataGrid_Operators_inc.hpp>


std::ostream& operator<<(
		const std::ostream &o_ostream,
		const sweet::Data::Cart2DComplex::DataGrid &i_dataArray
);


#endif
