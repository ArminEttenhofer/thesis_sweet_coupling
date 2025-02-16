/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#include "DataGrid.hpp"
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Memory/parmemcpy.hpp>
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
 * \brief Data stored in grid space
 */
void DataGrid::swap(
		DataGrid &i_cart2DData
)
{
	SWEET_ASSERT(cart2DDataConfig == i_cart2DData.cart2DDataConfig);

	std::swap(grid_space_data, i_cart2DData.grid_space_data);
}


DataGrid::DataGrid(
		const Cart2D::Config *i_cart2DDataConfig
)	:
	cart2DDataConfig(i_cart2DDataConfig),
	grid_space_data(nullptr)
{
	setup(i_cart2DDataConfig);
}



DataGrid::DataGrid()	:
	cart2DDataConfig(nullptr),
	grid_space_data(nullptr)
{
}



DataGrid::DataGrid(
		const DataGrid &i_data
)	:
	cart2DDataConfig(i_data.cart2DDataConfig),
	grid_space_data(nullptr)
{
	setup(i_data.cart2DDataConfig);

	operator=(i_data);
}



DataGrid::DataGrid(
		const Cart2D::DataGrid &i_data
)	:
	cart2DDataConfig(i_data.cart2DDataConfig),
	grid_space_data(nullptr)
{
	setup(i_data.cart2DDataConfig);

	operator=(i_data);
}


/*
 * load real and imaginary data from physical arrays
 */
/**
 * Run validation checks to make sure that the physical and spectral spaces match in size
 */

inline void DataGrid::check(
		const Cart2D::Config *i_cart2DDataConfig
)	const
{
	SWEET_ASSERT(cart2DDataConfig->grid_res[0] == i_cart2DDataConfig->grid_res[0]);
	SWEET_ASSERT(cart2DDataConfig->grid_res[1] == i_cart2DDataConfig->grid_res[1]);
}




void DataGrid::setup(
		const Cart2D::Config *i_cart2DDataConfig
)
{
	cart2DDataConfig = i_cart2DDataConfig;

	grid_space_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(cart2DDataConfig->grid_number_elements * sizeof(std::complex<double>));
}




void DataGrid::setup_if_required(
		const Cart2D::Config *i_cart2DDataConfig
)
{
	if (cart2DDataConfig != nullptr)
	{
		SWEET_ASSERT(grid_space_data != nullptr);
		return;
	}

	cart2DDataConfig = i_cart2DDataConfig;
	grid_space_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(cart2DDataConfig->grid_number_elements * sizeof(std::complex<double>));
}




DataGrid::~DataGrid()
{
	if (grid_space_data != nullptr)
		sweet::Memory::MemBlockAlloc::free(grid_space_data, cart2DDataConfig->grid_number_elements * sizeof(std::complex<double>));
}



void DataGrid::loadRealImag(
		const Cart2D::DataGrid &i_re,
		const Cart2D::DataGrid &i_im
)
{

	double *t = (double*)grid_space_data;
	double *in_re = (double*)i_re.grid_space_data;
	double *in_im = (double*)i_im.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
	{
		t[2*i] = in_re[i];
		t[2*i+1] = in_im[i];

	}
}


DataGrid& DataGrid::operator=(
		DataGrid &&i_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_data.cart2DDataConfig);

	std::swap(grid_space_data, i_data.grid_space_data);

	return *this;
}




DataGrid& DataGrid::operator=(
		const DataGrid &i_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_data.cart2DDataConfig);

	sweet::Memory::parmemcpy(grid_space_data, i_data.grid_space_data, sizeof(std::complex<double>)*cart2DDataConfig->grid_number_elements);

	return *this;
}



DataGrid& DataGrid::operator=(
		const Cart2D::DataGrid &i_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_data.cart2DDataConfig);

#if 1

	/*
	 * We implemented this double precision version
	 * since SIMDization with std::complex doesn't seem to work on LLVM
	 */
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		t[2*i] = in[i];
		t[2*i+1] = 0;
	}

#else

	for (std::size_t in = 0; in < cart2DDataConfig->grid_number_elements; in++)
		grid_space_data[in] = i_data.grid_space_data[in];

#endif

	return *this;
}



DataGrid DataGrid::operator+(
		const DataGrid &i_data
)	const
{
	check(i_data.cart2DDataConfig);

	DataGrid retval(cart2DDataConfig);

#if 1

	/*
	 * We implemented this double precision version
	 * since SIMDization with std::complex doesn't seem to work on LLVM
	 */
	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		out[i] = t[i] + in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i] + i_data.grid_space_data[i];

#endif

	return retval;
}



DataGrid& DataGrid::operator+=(
		const DataGrid &i_data
)
{
	check(i_data.cart2DDataConfig);

#if 1

	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		t[i] += in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] += i_data.grid_space_data[i];

#endif

	return *this;
}


DataGrid DataGrid::operator+(
		const std::complex<double> &i_value
)	const
{
	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i] + i_value.real();
		out[i+1] = t[i+1] + i_value.imag();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]+i_value;

#endif

	return retval;
}



DataGrid DataGrid::operator-(
		const std::complex<double> &i_value
)	const
{
	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i] - i_value.real();
		out[i+1] = t[i+1] - i_value.imag();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]-i_value;

#endif

	return retval;
}

DataGrid& DataGrid::operator-=(
		const DataGrid &i_data
)
{
	check(i_data.cart2DDataConfig);

#if 1

	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		t[i] -= in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] -= i_data.grid_space_data[i];

#endif
	return *this;
}



DataGrid DataGrid::operator-(
		const DataGrid &i_data
)	const
{
	check(i_data.cart2DDataConfig);

	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		out[i] = t[i] - in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i] - i_data.grid_space_data[i];

#endif

	return retval;
}



DataGrid DataGrid::operator-()
{
	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		out[i] = -t[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = -grid_space_data[i];

#endif

	return retval;
}



DataGrid DataGrid::operator*(
		const double i_value
)	const
{
	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		out[i] *= i_value;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]*i_value;

#endif

	return retval;
}



DataGrid DataGrid::operator*(
		const std::complex<double> &i_value
)	const
{
	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i]*i_value.real() - t[i+1]*i_value.imag();
		out[i+1] = t[i]*i_value.imag() + t[i+1]*i_value.real();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]*i_value;

#endif

	return retval;
}




DataGrid DataGrid::operator*(
		const DataGrid &i_data
)	const
{
	check(i_data.cart2DDataConfig);

	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i]*in[i] - t[i+1]*in[i+1];
		out[i+1] = t[i]*in[i+1] + t[i+1]*in[i];
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t in = 0; in < cart2DDataConfig->grid_number_elements; in++)
		retval.grid_space_data[in] = grid_space_data[in]*i_data.grid_space_data[in];

#endif

	return retval;
}



const DataGrid& DataGrid::operator*=(
		const double i_value
)	const
{
#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		t[i] *= i_value;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] *= i_value;

#endif

	return *this;
}



DataGrid DataGrid::operator/(
		const DataGrid &i_data
)	const
{
	check(i_data.cart2DDataConfig);

	DataGrid retval(cart2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	/*
	 * More tricky division
	 *
	 *   1          a-i*b          a-i*b
	 * ----- = --------------- = ---------
	 * a+i*b   (a+i*b)*(a-i*b)   a^2 + b^2
	 */
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		double &a = in[i];
		double &b = in[i+1];
		double denom = a*a + b*b;

		out[i] = (t[i]*a + t[i+1]*b)/denom;
		out[i+1] = (-t[i]*b + t[i+1]*a)/denom;
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		retval.grid_space_data[i] = grid_space_data[i]/i_data.grid_space_data[i];
	}

#endif

	return retval;
}



DataGrid DataGrid::operator/(
		double i_value
)	const
{
	DataGrid out(cart2DDataConfig);

#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		t[i] /= i_value;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		out.grid_space_data[i] = grid_space_data[i]/i_value;

#endif

	return out;
}




void DataGrid::grid_update_lambda_array_indices(
		std::function<void(int,int,std::complex<double>&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool i_anti_aliasing
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


void DataGrid::grid_update_lambda_array_idx(
		std::function<void(int,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		i_lambda(i, grid_space_data[i]);
	}
}


void DataGrid::grid_update_lambda_unit_coordinates_corner_centered(
		std::function<void(double,double,std::complex<double>&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool i_anti_aliasing
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

void DataGrid::grid_update_lambda_unit_coordinates_cell_centered(
		std::function<void(double,double,std::complex<double>&)> i_lambda,	//!< lambda function to return value for lat/mu
		bool i_anti_aliasing
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


/*
 * Set all values to zero
 */
void DataGrid::grid_setZero()
{
#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i++)
		t[i] = 0;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] = 0;

#endif
}


/*
 * Set all values to a specific value
 */
void DataGrid::grid_setValue(
		std::complex<double> &i_value
)
{
#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		t[i] = i_value.real();
		t[i+1] = i_value.imag();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] = i_value;
#endif
}

/*
 * Set all values to a specific value
 */
void DataGrid::grid_setValue(
		double i_value_real,
		double i_value_imag
)
{
#if 1
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->grid_number_elements; i+=2)
	{
		t[i] = i_value_real;
		t[i+1] = i_value_imag;
	}
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		grid_space_data[i].real(i_value_real);
		grid_space_data[i].imag(i_value_imag);
	}
#endif
}


/*
 * Set all values to a specific value
 */
void DataGrid::grid_setValue(
		int i_y_idx,
		int i_x_idx,
		std::complex<double> &i_value
)
{
	grid_space_data[i_y_idx*cart2DDataConfig->grid_res[0] + i_x_idx] = i_value;
}

/*
 * Set all values to a specific value
 */
void DataGrid::grid_setValue(
		int i_y_idx,
		int i_x_idx,
		double i_value_real,
		double i_value_imag
)
{
	grid_space_data[i_y_idx*cart2DDataConfig->grid_res[0] + i_x_idx].real(i_value_real);
	grid_space_data[i_y_idx*cart2DDataConfig->grid_res[0] + i_x_idx].imag(i_value_imag);
}


/**
 * Return the maximum error norm between this and the given data in physical space
 */
double DataGrid::grid_reduce_max(
		const DataGrid &i_data
)
{
	check(i_data.cart2DDataConfig);

	double error = -1;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
	{
		error = std::max(
					error,
					std::abs(
							grid_space_data[j] - i_data.grid_space_data[j]
						)
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
		std::complex<double> &d = grid_space_data[j];
		error += d.real()*d.real() + d.imag()*d.imag();
	}

	return std::sqrt(error / (double)cart2DDataConfig->grid_number_elements);
}


/**
 * Return the maximum error norm
 */
double DataGrid::grid_reduce_max_abs()	const
{
	double error = -1;
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(error)
	for (std::size_t j = 0; j < cart2DDataConfig->grid_number_elements; j++)
	{
		error = std::max(
					error,
					std::abs(grid_space_data[j].real())
					);

		error = std::max(
					error,
					std::abs(grid_space_data[j].imag())
					);
	}

	return error;
}

/**
 * return the maximum of all absolute values, use quad precision for reduction
 */
double DataGrid::grid_reduce_sum_re_quad()	const
{
	double sum = 0;
	double c = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		double value = grid_space_data[i].real();

		// Use Kahan summation
		double y = value - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	return sum;
}

/**
 * reduce to root mean square
 */
double DataGrid::grid_reduce_rms_quad()
{
	double sum = 0;
	double c = 0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		double radius2 = grid_space_data[i].real()*grid_space_data[i].real()+grid_space_data[i].imag()*grid_space_data[i].imag();
		double value = radius2;

		// Use Kahan summation
		double y = value - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	sum = std::sqrt(sum/double(cart2DDataConfig->grid_number_elements));
	return sum;
}

/**
 * return the sqrt of the sum of the squared values, use quad precision for reduction
 */
double DataGrid::grid_reduce_norm2_quad()	const
{
	double sum = 0.0;
	double c = 0.0;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		double value = grid_space_data[i].real()*grid_space_data[i].real() + grid_space_data[i].imag()*grid_space_data[i].imag();

		// Use Kahan summation
		double y = value - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	return std::sqrt(sum);
}


inline
std::ostream& operator<<(
		std::ostream &o_ostream,
		const DataGrid &i_dataArray
)
{
	DataGrid &rw_array_data = (DataGrid&)i_dataArray;

	for (int j = (int)rw_array_data.cart2DDataConfig->grid_data_size[1]-1; j >= 0; j--)
	{
		for (std::size_t i = 0; i < rw_array_data.cart2DDataConfig->grid_data_size[0]; i++)
		{
			std::cout << i_dataArray.grid_space_data[j*i_dataArray.cart2DDataConfig->grid_data_size[0]+i] << "\t";
		}
		std::cout << std::endl;
	}

	return o_ostream;
}

}}}

