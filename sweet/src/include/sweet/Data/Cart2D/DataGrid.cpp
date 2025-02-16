/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#include <complex>
#include <cfloat>
#include <cstddef>
#include <algorithm>
#include <memory>
#include <string.h>
#include <string>
#include <iostream>
#include <utility>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <functional>
#include <utility>
#include <cmath>
#include <iterator>

#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid_Kernels.hpp>
#include <sweet/Memory/parmemcpy.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Fatal.hpp>

#include "DataGrid.hpp"

namespace sweet {
namespace Data {
namespace Cart2D {


void DataGrid::swap(
		DataGrid &i_cart2DData
)
{
	SWEET_ASSERT(cart2DDataConfig == i_cart2DData.cart2DDataConfig);

	std::swap(grid_space_data, i_cart2DData.grid_space_data);
}


DataGrid::DataGrid(
		const Config *i_cart2DDataConfig
)	:
	// important: set this to nullptr, since a check for this will be performed by setup(...)
	cart2DDataConfig(i_cart2DDataConfig),
	grid_space_data(nullptr)
{
	alloc_data();
}



DataGrid::DataGrid(
		const Config &i_cart2DDataConfig
)	:
	// important: set this to nullptr, since a check for this will be performed by setup(...)
	cart2DDataConfig(&i_cart2DDataConfig),
	grid_space_data(nullptr)
{
	alloc_data();
}



DataGrid::DataGrid(
		const Config *i_cart2DDataConfig,
		double i_value
)	:
	//! important: set this to nullptr, since a check for this will be performed by setup(...)
	cart2DDataConfig(i_cart2DDataConfig),
	grid_space_data(nullptr)
{
	alloc_data();
	grid_setValue(i_value);
}



DataGrid::DataGrid()	:
	cart2DDataConfig(nullptr),
	grid_space_data(nullptr)
{
}

DataGrid::DataGrid(int i)	:
	cart2DDataConfig(nullptr),
	grid_space_data(nullptr)
{
}



DataGrid::DataGrid(
		const DataGrid &i_cart2d_data
)	:
	cart2DDataConfig(i_cart2d_data.cart2DDataConfig),
	grid_space_data(nullptr)

{
	if (i_cart2d_data.cart2DDataConfig != nullptr)
		alloc_data();

	operator=(i_cart2d_data);
}



DataGrid::DataGrid(
		DataGrid &&i_cart2d_data
)	:
	cart2DDataConfig(i_cart2d_data.cart2DDataConfig),
	grid_space_data(nullptr)
{
	if (i_cart2d_data.cart2DDataConfig == nullptr)
		return;

	std::swap(grid_space_data, i_cart2d_data.grid_space_data);
}



/*!
 * Run validation checks to make sure that the physical and spectral spaces match in size
 */

DataGrid::~DataGrid()
{
	clear();
}



void DataGrid::setup(
		const Config *i_cart2DDataConfig
)
{
	if (cart2DDataConfig != nullptr)
		SWEETErrorFatal("Setup called twice!");

	cart2DDataConfig = i_cart2DDataConfig;
	alloc_data();
}



void DataGrid::setup(
		const Config &i_cart2DDataConfig
)
{
	setup(&i_cart2DDataConfig);
}



void DataGrid::setup_if_required(
		const Config *i_cart2DDataConfig
)
{
	if (cart2DDataConfig != nullptr)
		return;

	setup(i_cart2DDataConfig);
}



DataGrid& DataGrid::operator=(
		const DataGrid &i_cart2d_data
)
{
	if (i_cart2d_data.cart2DDataConfig == nullptr)
		return *this;

	if (cart2DDataConfig == nullptr)
		setup(i_cart2d_data.cart2DDataConfig);

	sweet::Memory::parmemcpy(grid_space_data, i_cart2d_data.grid_space_data, sizeof(double)*cart2DDataConfig->grid_number_elements);

	dealiasing(*this);

	return *this;
}



void DataGrid::check(
		const Config *i_cart2DDataConfig
)	const
{
	SWEET_ASSERT(cart2DDataConfig->grid_res[0] == i_cart2DDataConfig->grid_res[0]);
	SWEET_ASSERT(cart2DDataConfig->grid_res[1] == i_cart2DDataConfig->grid_res[1]);
}




DataGrid& DataGrid::operator=(
		DataGrid &&i_cart2d_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_cart2d_data.cart2DDataConfig);

	std::swap(grid_space_data, i_cart2d_data.grid_space_data);

	dealiasing(*this);

	return *this;
}



DataGrid &DataGrid::operator=(double i_value)
{
	grid_setValue(i_value);

	return *this;
}



DataGrid &DataGrid::operator=(int i_value)
{
	grid_setValue(i_value);

	return *this;
}







void DataGrid::alloc_data()
{
	SWEET_ASSERT(grid_space_data == nullptr);
	grid_space_data = sweet::Memory::MemBlockAlloc::alloc<double>(cart2DDataConfig->grid_number_elements * sizeof(double));
}



void DataGrid::clear()
{
	if (grid_space_data != nullptr)
	{
		sweet::Memory::MemBlockAlloc::free(grid_space_data, cart2DDataConfig->grid_number_elements * sizeof(double));
		grid_space_data = nullptr;
	}

	cart2DDataConfig = nullptr;
}



DataGrid DataGrid::operator+(
		double i_value
)	const
{
	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = grid_space_data[idx]+i_value;

	return out;
}

DataGrid DataGrid::operator+(
		const DataGrid &i_cart2d_data
)	const
{
	check(i_cart2d_data.cart2DDataConfig);

	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = grid_space_data[idx] + i_cart2d_data.grid_space_data[idx];

#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
	dealiasing(out);
#endif

	return out;
}


DataGrid& DataGrid::operator+=(
		const DataGrid &i_cart2d_data
)
{
	check(i_cart2d_data.cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		grid_space_data[idx] += i_cart2d_data.grid_space_data[idx];

	return *this;
}


DataGrid& DataGrid::operator+=(
		double i_scalar
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		grid_space_data[idx] += i_scalar;

	return *this;
}


DataGrid DataGrid::operator-(
		double i_value
)	const
{
	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = grid_space_data[idx]-i_value;

	return out;
}

DataGrid& DataGrid::operator-=(
		const DataGrid &i_cart2d_data
)
{
	check(i_cart2d_data.cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		grid_space_data[idx] -= i_cart2d_data.grid_space_data[idx];

	return *this;
}


DataGrid& DataGrid::operator-=(
		double i_scalar
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		grid_space_data[idx] -= i_scalar;

	return *this;
}


DataGrid DataGrid::operator-(
		const DataGrid &i_cart2d_data
)	const
{
	check(i_cart2d_data.cart2DDataConfig);

	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = grid_space_data[idx] - i_cart2d_data.grid_space_data[idx];

#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
	dealiasing(out);
#endif

	return out;
}


DataGrid DataGrid::operator-()
{
	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = -grid_space_data[idx];

	return out;
}



DataGrid DataGrid::operator*(
		const DataGrid &i_cart2d_data
)	const
{
	check(i_cart2d_data.cart2DDataConfig);


	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		out.grid_space_data[i] = grid_space_data[i]*i_cart2d_data.grid_space_data[i];


#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
	dealiasing(out);
#endif

	return out;

}


DataGrid DataGrid::operator*(
		const double i_value
)	const
{
	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		out.grid_space_data[i] = grid_space_data[i]*i_value;

	return out;
}

const DataGrid& DataGrid::operator*=(
		const double i_value
)	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		grid_space_data[idx] *= i_value;

	return *this;
}

DataGrid DataGrid::operator/(
		double i_value
)	const
{
	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = grid_space_data[idx]/i_value;

	return out;
}

DataGrid& DataGrid::operator/=(
		double i_scalar
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		grid_space_data[idx] /= i_scalar;

	return *this;
}

DataGrid DataGrid::operator()(
		const DataGrid &i_array_data
)	const
{
	DataGrid out(cart2DDataConfig);

	DataGrid &rw_array_data = (DataGrid&)i_array_data;

	kernel_apply(
			cart2DDataConfig->grid_data_size[0],
			cart2DDataConfig->grid_data_size[1],
			rw_array_data.grid_space_data,

			out.grid_space_data
	);


	return out;
}


DataGrid DataGrid::operator_scalar_sub_this(
		double i_value
)	const
{
	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		out.grid_space_data[idx] = i_value - grid_space_data[idx];

	return out;
}


DataGrid DataGrid::operator/(
		const DataGrid &i_cart2d_data
)	const
{
	check(i_cart2d_data.cart2DDataConfig);

	check(i_cart2d_data.cart2DDataConfig);

	DataGrid out(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		out.grid_space_data[i] = grid_space_data[i]/i_cart2d_data.grid_space_data[i];

	dealiasing(out);

	return out;
}


void DataGrid::dealiasing(
	DataGrid& io_data
) const
{
	// create spectral data container
	std::complex<double> *spectral_space_data = nullptr;
	spectral_space_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(cart2DDataConfig->spectral_array_data_number_of_elements * sizeof(std::complex<double>));

	// FFT
	cart2DDataConfig->fft_grid_2_spectral_OUTOFPLACE(io_data.grid_space_data, spectral_space_data);

	// Dealiasing
	//SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int k = 0; k < 2; k++)
	{
		if (k == 0)
		{
			/*
			 * First process part between top and bottom spectral data blocks
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
			for (std::size_t jj = cart2DDataConfig->spectral_data_iteration_ranges[0][1][1]; jj < cart2DDataConfig->spectral_data_iteration_ranges[1][1][0]; jj++)
				for (std::size_t ii = cart2DDataConfig->spectral_data_iteration_ranges[0][0][0]; ii < cart2DDataConfig->spectral_data_iteration_ranges[0][0][1]; ii++)
				{
					//spectral_space_data[jj*cart2DDataConfig->spectral_data_size[0]+ii] = 0;
					spectral_space_data[cart2DDataConfig->getArrayIndexByModes(jj, ii)] = 0;
				}
		}
		else
		{
			/*
			 * Then process the aliasing block on the right side
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
			for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_data_size[1]; jj++)
				for (std::size_t ii = cart2DDataConfig->spectral_data_iteration_ranges[0][0][1]; ii < cart2DDataConfig->spectral_data_size[0]; ii++)
				{
					//spectral_space_data[jj*cart2DDataConfig->spectral_data_size[0]+ii] = 0;
					spectral_space_data[cart2DDataConfig->getArrayIndexByModes(jj, ii)] = 0;
				}
		}
	}

	// IFFT
	cart2DDataConfig->fft_spectral_2_grid_INPLACE(spectral_space_data, io_data.grid_space_data);

	// Free spectral data
	sweet::Memory::MemBlockAlloc::free(spectral_space_data, cart2DDataConfig->spectral_array_data_number_of_elements * sizeof(std::complex<double>));
	spectral_space_data = nullptr;
}


double DataGrid::grid_get(
		int i_y_idx,
		int i_x_idx
)	const
{
	return grid_space_data[i_y_idx*cart2DDataConfig->grid_res[0] + i_x_idx];
}


void DataGrid::grid_update_lambda_unit_coordinates_cell_centered(
		std::function<void(double,double,double&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool anti_aliasing
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t j = 0; j < cart2DDataConfig->grid_data_size[1]; j++)
	{
		for (std::size_t i = 0; i < cart2DDataConfig->grid_data_size[0]; i++)
		{
			std::size_t idx = j*cart2DDataConfig->grid_data_size[0]+i;

			i_lambda(
					((double)i+0.5)/(double)cart2DDataConfig->grid_res[0],
					((double)j+0.5)/(double)cart2DDataConfig->grid_res[1],
					grid_space_data[idx]
			);
		}
	}
}


void DataGrid::grid_update_lambda_unit_coordinates_corner_centered(
		std::function<void(double,double,double&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool anti_aliasing
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t j = 0; j < cart2DDataConfig->grid_data_size[1]; j++)
	{
		for (std::size_t i = 0; i < cart2DDataConfig->grid_data_size[0]; i++)
		{
			std::size_t idx = j*cart2DDataConfig->grid_data_size[0]+i;

			i_lambda(
					(double)i/(double)cart2DDataConfig->grid_res[0],
					(double)j/(double)cart2DDataConfig->grid_res[1],
					grid_space_data[idx]
			);
		}
	}
}

void DataGrid::grid_setValue(
		double i_value
)
{

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t j = 0; j < cart2DDataConfig->grid_data_size[1]; j++)
	{
		for (std::size_t i = 0; i < cart2DDataConfig->grid_data_size[0]; i++)
		{
			std::size_t idx = j*cart2DDataConfig->grid_data_size[0]+i;
			grid_space_data[idx] = i_value;
		}
	}
}


void DataGrid::grid_update_lambda_array_idx(
		std::function<void(int,double&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool anti_aliasing
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR

	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		i_lambda(i, grid_space_data[i]);
	}
}


void DataGrid::grid_update_lambda_array_indices(
		std::function<void(int,int,double&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool anti_aliasing
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t j = 0; j < cart2DDataConfig->grid_data_size[1]; j++)
	{
		for (std::size_t i = 0; i < cart2DDataConfig->grid_data_size[0]; i++)
		{
			std::size_t idx = j*cart2DDataConfig->grid_data_size[0]+i;

			i_lambda(i, j, grid_space_data[idx]);
		}
	}
}


__attribute__((deprecated))
void DataGrid::grid_setValue(
		int i_y_idx,
		int i_x_idx,
		double i_value
)
{
	grid_space_data[i_y_idx*cart2DDataConfig->grid_res[0] + i_x_idx] = i_value;
}

void DataGrid::grid_setZero()
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t j = 0; j < cart2DDataConfig->grid_data_size[1]; j++)
	{
		for (std::size_t i = 0; i < cart2DDataConfig->grid_data_size[0]; i++)
		{
			std::size_t idx = j*cart2DDataConfig->grid_data_size[0]+i;

			grid_space_data[idx] = 0;
		}
	}
}


double DataGrid::grid_reduce_max(
		const DataGrid &i_cart2d_data
)
{
	check(i_cart2d_data.cart2DDataConfig);

	double error = -1;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
	{
		error = std::max(
					std::abs(
							grid_space_data[j] - i_cart2d_data.grid_space_data[j]
						),
						error
					);
	}
	return error;
}


double DataGrid::grid_reduce_rms()
{
	double error = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
	{
		double &d = grid_space_data[j];
		error += d*d;
	}

	return std::sqrt(error / (double)cart2DDataConfig->grid_number_elements);
}



double DataGrid::grid_reduce_sum()	const
{
	double sum = 0;
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
		sum += grid_space_data[j];

	return sum;
}


double DataGrid::grid_reduce_sum_quad()	const
{
	double sum = 0;
	double c = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		double value = grid_space_data[i];

		// Use Kahan summation
		double y = value - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	return sum;
}



double DataGrid::grid_reduce_rms_quad()
{

	double sum = 0;
	double c = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		double value = grid_space_data[i]*grid_space_data[i];

		// Use Kahan summation
		double y = value - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	sum = std::sqrt(sum/(double)(cart2DDataConfig->grid_number_elements));

	return sum;
}


double DataGrid::grid_reduce_norm1()	const
{
	double sum = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		sum += std::abs(grid_space_data[i]);


	return sum;
}


double DataGrid::grid_reduce_norm2()	const
{
	double sum = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		sum += grid_space_data[i]*grid_space_data[i];


	return std::sqrt(sum);
}


double DataGrid::grid_reduce_norm2_quad()	const
{
	double sum = 0.0;
	double c = 0.0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		double value = grid_space_data[i]*grid_space_data[i];

		// Use Kahan summation
		double y = value - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	return std::sqrt(sum);
}


double DataGrid::grid_reduce_max_abs(
		const DataGrid &i_cart2d_data
)	const
{
	check(i_cart2d_data.cart2DDataConfig);

	double error = -1;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
	{
		error = std::max(
					std::abs(
							grid_space_data[j] - i_cart2d_data.grid_space_data[j]
						),
						error	// leave the error variable as the 2nd parameter. In case of NaN of the 1st parameter, std::max returns NaN
					);
	}
	return error;
}


double DataGrid::grid_reduce_max_abs()	const
{
	double error = -1;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
	{
		error = std::max(
					std::abs(grid_space_data[j]),
					error		// leave the error variable as the 2nd parameter. In case of NaN of the 1st parameter, std::max returns NaN
			);
	}
	return error;
}


double DataGrid::grid_reduce_min()	const
{
	double error = std::numeric_limits<double>::infinity();

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MIN(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
		error = std::min(grid_space_data[j], error);

	return error;
}


double DataGrid::grid_reduce_max()	const
{
	double error = -std::numeric_limits<double>::infinity();

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
		error = std::max(grid_space_data[j], error);

	return error;
}

bool DataGrid::grid_isAnyNaNorInf()	const
{
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		if (std::isnan(grid_space_data[i]) || std::isinf(grid_space_data[i]) != 0)
			return true;
	}

	return false;
}



bool DataGrid::grid_reduce_boolean_all_finite() const
{

	bool isallfinite = true;

#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(&&:isallfinite)
#endif
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		isallfinite = isallfinite && std::isfinite(grid_space_data[i]);

	return isallfinite;
}


void DataGrid::grid_file_write(
		const std::string &i_filename,
		double cart2d_domain_size[2],
		const char *i_title,
		int i_precision
)	const
{
	std::ofstream file(i_filename, std::ios_base::trunc);

	if (i_precision >= 0)
		file << std::setprecision(i_precision);

	file << "#TI " << i_title << std::endl;
	file << "#TX" << std::endl;
	file << "#TY" << std::endl;

	//file << "lat\\lon\t";
	// Use 0 to make it processable by python
	file << "0\t";

	for (int i = 0; i < (int)cart2DDataConfig->grid_res[0]; i++)
	{
//			double lon_degree = ((double)i/(double)cart2DDataConfig->spat_num_lon)*2.0*M_PI;
		double x = ((double)i/(double)cart2DDataConfig->grid_res[0])*cart2d_domain_size[0]; // ????

		file << x;
		if (i < (int)cart2DDataConfig->grid_res[0]-1)
			file << "\t";
	}
	file << std::endl;

	for (int j = (int)cart2DDataConfig->grid_res[1]-1; j >= 0; j--)
	{
		double y = ((double)j/(double)cart2DDataConfig->grid_res[1])*cart2d_domain_size[1]; // ????

		file << y << "\t";

		for (int i = 0; i < (int)cart2DDataConfig->grid_res[0]; i++)
		{
			file << grid_space_data[j*cart2DDataConfig->grid_res[0]+i];
			if (i < (int)cart2DDataConfig->grid_res[0]-1)
				file << "\t";
		}
		file << std::endl;
	}
	file.close();
}



bool DataGrid::file_grid_saveData_ascii(
		const char *i_filename,		//!< Name of file to store data to
		char i_separator,			//!< separator to use for each line
		int i_precision,			//!< number of floating point digits
		int dimension				//!< store 1D or 2D
)	const
{

	std::ofstream file(i_filename, std::ios_base::trunc);
	file << std::setprecision(i_precision);

	file << "#SWEET" << std::endl;
	file << "#FORMAT ASCII" << std::endl;
	file << "#PRIMITIVE PLANE" << std::endl;

	std::size_t resx = cart2DDataConfig->grid_res[0];
	std::size_t resy = cart2DDataConfig->grid_res[1];

	file << "#SPACE PHYSICAL" << std::endl;
	file << "#RESX " << resx << std::endl;
	file << "#RESY " << resy << std::endl;

	std::size_t ymin = 0;
	if (dimension == 2)
		ymin = 0;
	else
		ymin = cart2DDataConfig->grid_res[1]-1;

	for (int y = (int) resy-1; y >= (int) ymin; y--)
	{
		for (std::size_t x = 0; x < resx; x++)
		{
			file << grid_get(y, x);

			if (x < cart2DDataConfig->grid_res[0]-1)
				file << i_separator;
			else
				file << std::endl;
		}
	}


	return true;
}


bool DataGrid::file_read_csv_grid(
		const char *i_filename		//!< Name of file to load data from
)
{

	Config cart2DDataConfig_ref;

	std::cout << "loading DATA from " << i_filename << std::endl;

	std::ifstream file(i_filename);


	int resx_ref = -1;
	int resy_ref = -1;
	for (int i = 0; i < 6; i++)
	{
		std::string line;
		std::getline(file, line);
		std::istringstream iss(line);
		std::vector<std::string> str_vector((std::istream_iterator<std::string>(iss)),
			std::istream_iterator<std::string>());

		if (i == 0)
		{
			SWEET_ASSERT(str_vector.size() == 1);
			SWEET_ASSERT(str_vector[0] == "#SWEET");
		}
		else if (i == 1)
		{
			SWEET_ASSERT(str_vector.size() == 2);
			SWEET_ASSERT(str_vector[0] == "#FORMAT");
			SWEET_ASSERT(str_vector[1] == "ASCII");
		}
		else if (i == 2)
		{
			SWEET_ASSERT(str_vector.size() == 2);
			SWEET_ASSERT(str_vector[0] == "#PRIMITIVE");
			SWEET_ASSERT(str_vector[1] == "PLANE");
		}
		else if (i == 3)
		{
			SWEET_ASSERT(str_vector.size() == 2);
			SWEET_ASSERT(str_vector[0] == "#SPACE");
			SWEET_ASSERT(str_vector[1] == "PHYSICAL");
		}
		else if (i == 4)
		{
			SWEET_ASSERT(str_vector.size() == 2);
			SWEET_ASSERT(str_vector[0] == "#RESX");
			resx_ref = stoi(str_vector[1]);
		}
		else if (i == 5)
		{
			SWEET_ASSERT(str_vector.size() == 2);
			SWEET_ASSERT(str_vector[0] == "#RESY");
			resy_ref = stoi(str_vector[1]);
		}
	}

	std::cout << "RESXY " << resx_ref << " " << resy_ref << std::endl;

	cart2DDataConfig_ref.setup(resx_ref, resy_ref, (int)((resx_ref * 2) / 3), (int)((resx_ref * 2) / 3), cart2DDataConfig->reuse_spectral_transformation_plans);
	*this = DataGrid(&cart2DDataConfig_ref);

	for (std::size_t row = 0; row < cart2DDataConfig_ref.grid_res[1]; row++)
	{
		std::string line;
		std::getline(file, line);
		if (!file.good())
		{
			std::cerr << "Failed to read data from file " << i_filename << " in line " << row << std::endl;
			return false;
		}

		std::size_t last_pos = 0;
		std::size_t col = 0;
		for (std::size_t pos = 0; pos < line.size()+1; pos++)
		{
			if (pos < line.size())
				if (line[pos] != '\t' && line[pos] != ' ')
					continue;

			std::string strvalue = line.substr(last_pos, pos-last_pos);

			double i_value = atof(strvalue.c_str());

			grid_setValue(cart2DDataConfig_ref.grid_res[1]-row-1, col, i_value);

			col++;
			last_pos = pos+1;
		}

		if (col < cart2DDataConfig_ref.grid_res[0])
		{
			std::cerr << "Failed to read data from file " << i_filename << " in line " << row << ", column " << col << std::endl;
			return false;
		}
	}

	std::cout << "DATA loaded OK" << std::endl;

	return true;
}


void DataGrid::file_write_raw(
		const std::string &i_filename
)	const
{
	std::fstream file(i_filename, std::ios::out | std::ios::binary);
	file.write((const char*)grid_space_data, sizeof(double)*cart2DDataConfig->grid_number_elements);
}



void DataGrid::file_read_raw(
		const std::string &i_filename
)	const
{
	std::fstream file(i_filename, std::ios::in | std::ios::binary);
	file.read((char*)grid_space_data, sizeof(double)*cart2DDataConfig->grid_number_elements);
}


void DataGrid::print_physicalData_zeroNumZero(
		double i_zero_threshold
)	const
{
	DataGrid &rw_array_data = (DataGrid&)*this;

	for (int y = cart2DDataConfig->grid_data_size[1]-1; y >= 0; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->grid_data_size[0]; x++)
		{
			double value = rw_array_data.grid_get(y, x);

			if (std::abs(value) < i_zero_threshold)
				value = 0.0;

			std::cout << value;

			if (x != cart2DDataConfig->grid_data_size[0]-1)
				std::cout << "\t";
		}
		std::cout << std::endl;
	}
}


void DataGrid::grid_print(
		int i_precision
)	const
{
	if (i_precision >= 0)
		std::cout << std::setprecision(i_precision);

	for (int j = (int)(cart2DDataConfig->grid_res[1]-1); j >= 0; j--)
	{
		for (int i = 0; i < (int)cart2DDataConfig->grid_res[0]; i++)
		{
			std::cout << grid_space_data[j*cart2DDataConfig->grid_res[0]+i];
			if (i < (int)cart2DDataConfig->grid_res[0]-1)
				std::cout << "\t";
		}
		std::cout << std::endl;
	}
}

void DataGrid::print_debug(
		const char *name
)	const
{
	std::cout << name << ":" << std::endl;
	std::cout << "                min: " << grid_reduce_min() << std::endl;
	std::cout << "                max: " << grid_reduce_max() << std::endl;
	std::cout << "                sum: " << grid_reduce_sum() << std::endl;
	std::cout << std::endl;
}


void DataGrid::print()	const
{
	for (std::size_t idx = 0; idx < cart2DDataConfig->grid_number_elements; idx++)
		std::cout << grid_space_data[idx] << "\t";
	std::cout << std::endl;
}

}}}
