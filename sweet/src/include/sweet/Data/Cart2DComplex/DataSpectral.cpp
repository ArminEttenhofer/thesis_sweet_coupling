/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#include <complex>
#include <cstddef>
#include <memory>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <utility>
#include <limits>
#include <fstream>
#include <iomanip>
#include <functional>
#include <cmath>

#include "DataSpectral.hpp"
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Memory/parmemcpy.hpp>



namespace sweet {
namespace Data {
namespace Cart2DComplex {


DataSpectral::DataSpectral()	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
}



DataSpectral::DataSpectral(
		const Cart2D::Config *i_cart2DDataConfig
)	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_cart2DDataConfig != 0);

	setup(i_cart2DDataConfig);
}



DataSpectral::DataSpectral(
		const Cart2D::Config &i_cart2DDataConfig
)	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	setup(&i_cart2DDataConfig);
}



DataSpectral::DataSpectral(
		const DataSpectral &i_cart2d_data
)	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	setup(i_cart2d_data.cart2DDataConfig);

	operator=(i_cart2d_data);
}



DataSpectral::DataSpectral(
		DataSpectral &&i_data
)	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	if (cart2DDataConfig == nullptr)
		setup(i_data.cart2DDataConfig);

	SWEET_ASSERT(i_data.spectral_space_data != nullptr);

	std::swap(spectral_space_data, i_data.spectral_space_data);
}





DataSpectral::DataSpectral(
		const DataGrid &i_cart2d_data
):
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	setup(i_cart2d_data.cart2DDataConfig);

	cart2DDataConfig->fft_complex_grid_2_spectral_OUTOFPLACE(i_cart2d_data.grid_space_data, spectral_space_data);
}


DataSpectral::~DataSpectral()
{
	if (spectral_space_data != nullptr)
	{
		sweet::Memory::MemBlockAlloc::free(spectral_space_data, cart2DDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
	}
}


DataGrid DataSpectral::toGrid()	const
{
	return getCart2DDataGridComplex();
}

DataGrid DataSpectral::getCart2DDataGridComplex()	const
{
	DataGrid out(cart2DDataConfig);

	/*
	 * WARNING:
	 * We have to use a temporary array here because of destructive FFT transformations
	 */
	spectral_zeroAliasingModes();
	DataSpectral tmp = *this;
	cart2DDataConfig->fft_complex_spectral_2_grid_OUTOFPLACE(tmp.spectral_space_data, out.grid_space_data);

	return out;
}


void DataSpectral::check_cart2DDataConfig_identical_res(
		 const Cart2D::Config *i_cart2DDataConfig
)	const
{
	SWEET_ASSERT(cart2DDataConfig->spectral_complex_data_size[0] == i_cart2DDataConfig->spectral_complex_data_size[0]);
	SWEET_ASSERT(cart2DDataConfig->spectral_complex_data_size[1] == i_cart2DDataConfig->spectral_complex_data_size[1]);
}


std::complex<double>& DataSpectral::operator[](std::size_t i)
{
	return spectral_space_data[i];
}

const std::complex<double>& DataSpectral::operator[](std::size_t i)	const
{
	return spectral_space_data[i];
}



DataSpectral& DataSpectral::operator=(
		const DataSpectral &i_cart2d_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_cart2d_data.cart2DDataConfig);

	SWEET_ASSERT(i_cart2d_data.spectral_space_data);
	sweet::Memory::parmemcpy(spectral_space_data, i_cart2d_data.spectral_space_data, sizeof(Tcomplex)*cart2DDataConfig->spectral_complex_array_data_number_of_elements);

	return *this;
}


DataSpectral& DataSpectral::operator=(
		DataSpectral &&i_cart2d_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_cart2d_data.cart2DDataConfig);

	SWEET_ASSERT(i_cart2d_data.spectral_space_data);
	std::swap(spectral_space_data, i_cart2d_data.spectral_space_data);

	return *this;
}


DataSpectral DataSpectral::operator+(
		double i_value
)	const
{
	DataSpectral out_cart2d_data(*this);

	out_cart2d_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_cart2d_data;
}


DataSpectral DataSpectral::operator+(
		const DataSpectral &i_cart2d_data
)	const
{
	check_cart2DDataConfig_identical_res(i_cart2d_data.cart2DDataConfig);

	DataSpectral out_cart2d_data(cart2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx] + i_cart2d_data.spectral_space_data[idx];

	out_cart2d_data.spectral_zeroAliasingModes();

	return out_cart2d_data;
}



DataSpectral& DataSpectral::operator+=(
		double i_value
)
{
	spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
	return *this;
}


DataSpectral& DataSpectral::operator+=(
		const DataSpectral &i_cart2d_data
)
{
	check_cart2DDataConfig_identical_res(i_cart2d_data.cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] += i_cart2d_data.spectral_space_data[idx];

	return *this;
}


DataSpectral DataSpectral::operator+(
		const std::complex<double> &i_value
)	const
{
	DataSpectral out_cart2d_data(*this);

	out_cart2d_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_cart2d_data;
}



DataSpectral DataSpectral::operator-(
		const std::complex<double> &i_value
)	const
{
	DataSpectral out_cart2d_data(*this);

	out_cart2d_data.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);

	return out_cart2d_data;
}


DataSpectral DataSpectral::operator-(
		const DataSpectral &i_cart2d_data
)	const
{
	check_cart2DDataConfig_identical_res(i_cart2d_data.cart2DDataConfig);

	DataSpectral out_cart2d_data(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx] - i_cart2d_data.spectral_space_data[idx];

	out_cart2d_data.spectral_zeroAliasingModes();

	return out_cart2d_data;
}


DataSpectral& DataSpectral::operator-=(
		const DataSpectral &i_cart2d_data
)
{
	check_cart2DDataConfig_identical_res(i_cart2d_data.cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] -= i_cart2d_data.spectral_space_data[idx];

	return *this;
}


DataSpectral DataSpectral::operator-()	const
{
	DataSpectral out_cart2d_data(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = -spectral_space_data[idx];

	return out_cart2d_data;
}


DataSpectral DataSpectral::operator*(
		const double i_value
)	const
{
	DataSpectral out_cart2d_data(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR

	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

	return out_cart2d_data;
}


const DataSpectral& DataSpectral::operator*=(
		const std::complex<double> &i_value
)	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] *= i_value;

	return *this;
}

DataSpectral DataSpectral::operator*(
		const DataSpectral &i_cart2d_data
)	const
{
	check_cart2DDataConfig_identical_res(i_cart2d_data.cart2DDataConfig);

	DataSpectral out_cart2d_data(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx] * i_cart2d_data.spectral_space_data[idx];

	out_cart2d_data.spectral_zeroAliasingModes();

	return out_cart2d_data;
}


DataSpectral DataSpectral::operator*(
		const std::complex<double> &i_value
)	const
{
	DataSpectral out_cart2d_data(cart2DDataConfig);


SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

	return out_cart2d_data;
}


DataSpectral DataSpectral::operator/(
		double i_value
)	const
{
	DataSpectral out_cart2d_data(cart2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

	return out_cart2d_data;
}

DataSpectral DataSpectral::operator/(
		const DataSpectral &i_cart2d_data
)	const
{
	check_cart2DDataConfig_identical_res(i_cart2d_data.cart2DDataConfig);

	DataSpectral out_cart2d_data(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
	{
		if (i_cart2d_data.spectral_space_data[idx] == std::complex<double>(0,0))
			out_cart2d_data.spectral_space_data[idx] = 0;
		else
			out_cart2d_data.spectral_space_data[idx] = spectral_space_data[idx] / i_cart2d_data.spectral_space_data[idx];
	}

	out_cart2d_data.spectral_zeroAliasingModes();


	return out_cart2d_data;
}


const DataSpectral& DataSpectral::operator/=(
		const std::complex<double> &i_value
)	const
{

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] /= i_value;

	return *this;
}


DataSpectral DataSpectral::operator()(
		const DataSpectral &i_array_data
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_array_data.spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
	{
		for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(jj, ii);
			double ar = t[2*idx];
			double ai = t[2*idx+1];
			double br = in[2*idx];
			double bi = in[2*idx+1];

			out[2*idx] = ar*br - ai*bi;
			out[2*idx+1] = ar*bi + ai*br;
		}
	}

#else
	CART2D_DATA_COMPLEX_SPECTRAL_FOR_IDX(
			retval.spectral_space_data[idx] = spectral_space_data[idx]*i_array_data.spectral_space_data[idx];
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


void DataSpectral::setup_if_required(
	const Cart2D::Config *i_cart2DDataConfig
)
{
	if (cart2DDataConfig != nullptr)
		return;

	setup(i_cart2DDataConfig);
}


void DataSpectral::setup(
		const Cart2D::Config *i_cart2dConfig
)
{
	// assure that the initialization is not done twice!
	SWEET_ASSERT(cart2DDataConfig == nullptr);

	cart2DDataConfig = i_cart2dConfig;

	spectral_space_data = sweet::Memory::MemBlockAlloc::alloc<Tcomplex>(cart2DDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
}


DataSpectral DataSpectral::spectral_solve_helmholtz(
		const std::complex<double> &i_a,
		const std::complex<double> &i_b,
		double r
)	const
{
	DataSpectral out(*this);

	const std::complex<double> a = i_a;
	const std::complex<double> b = i_b;

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			io_data /= (a + (-b*((double)(m * m + n * n))));
		}
	);

	return out;
}


void DataSpectral::spectral_update_lambda(
		std::function<void(int,int,Tcomplex&)> i_lambda
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t n = 0; n < cart2DDataConfig->spectral_complex_data_size[1]; n++)
	{
		for (std::size_t m = 0; m < cart2DDataConfig->spectral_complex_data_size[0]; m++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(n, m);
			i_lambda(n, m, spectral_space_data[idx]);
		}
	}
}

void DataSpectral::spectral_update_lambda_modes(
		std::function<void(int,int,std::complex<double>&)> i_lambda
)
{
	int modes_1 = cart2DDataConfig->spectral_complex_data_size[1];

	int half_modes_1 = modes_1/2;

#if 1
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
	{
		for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(jj, ii);

			int k0 = ii;

			int k1 = jj;
			if (k1 > half_modes_1)
				k1 = k1 - modes_1;

			i_lambda(k0, k1, spectral_space_data[idx]);
		}
	}
#else
	CART2D_DATA_COMPLEX_SPECTRAL_FOR_IDX(
		{
			int k0 = ii;

			int k1 = jj;
			if (k1 > half_modes_1)
				k1 = k1 - modes_1;

			i_lambda(k0, k1, spectral_space_data[idx]);
		}
	);
#endif
	spectral_zeroAliasingModes();
}


const std::complex<double>& DataSpectral::spectral_get(
		int in,
		int im
)	const
{
	SWEET_ASSERT(in >= 0);
	SWEET_ASSERT(in <= (int)cart2DDataConfig->spectral_complex_data_size[1]);
	SWEET_ASSERT(std::abs(im) <= (int)cart2DDataConfig->spectral_complex_data_size[0]);

	return spectral_space_data[cart2DDataConfig->getArrayIndexByModes_Complex(in, im)];
}


void DataSpectral::spectral_setZero()
{

	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_complex_data_size[1]; n++)
		{
			for (std::size_t m = 0; m < cart2DDataConfig->spectral_complex_data_size[0]; m++)
			{
				int idx = cart2DDataConfig->getArrayIndexByModes_Complex(n, m);
				spectral_space_data[idx] = 0;
			}
		}
}


void DataSpectral::spectral_set(
		int i_n,
		int i_m,
		double i_real,
		double i_imag
)	const
{
#if SWEET_DEBUG
	if (i_n < 0 ||  i_m < 0)
		SWEETErrorFatal("Out of boundary a");

	if (i_m >= (int)cart2DDataConfig->spectral_complex_data_size[0])
		SWEETErrorFatal("Out of boundary b");

	if (i_n >= (int)cart2DDataConfig->spectral_complex_data_size[1])
		SWEETErrorFatal("Out of boundary c");
#endif

	spectral_space_data[cart2DDataConfig->getArrayIndexByModes_Complex(i_n, i_m)].real(i_real);
	spectral_space_data[cart2DDataConfig->getArrayIndexByModes_Complex(i_n, i_m)].imag(i_imag);
}


void DataSpectral::spectral_set(
		int i_n,
		int i_m,
		std::complex<double> i_data
)	const
{
#if SWEET_DEBUG
	if (i_n < 0 ||  i_m < 0)
		SWEETErrorFatal("Out of boundary a");

	if (i_m >= (int)cart2DDataConfig->spectral_complex_data_size[0])
		SWEETErrorFatal("Out of boundary b");

	if (i_n >= (int)cart2DDataConfig->spectral_complex_data_size[1])
		SWEETErrorFatal("Out of boundary c");
#endif

	spectral_space_data[cart2DDataConfig->getArrayIndexByModes_Complex(i_n, i_m)].real(i_data.real());
	spectral_space_data[cart2DDataConfig->getArrayIndexByModes_Complex(i_n, i_m)].imag(i_data.imag());
}


DataSpectral DataSpectral::spectral_addScalarAll(
		const double &i_value
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;


	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
	{
		for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(jj, ii);

			out[2*idx] = t[2*idx] + i_value;
			out[2*idx+1] = t[2*idx+1] + i_value;
		}
	}

#else
	CART2D_DATA_COMPLEX_SPECTRAL_FOR_IDX(
			retval.spectral_space_data[idx] = spectral_space_data[idx] + i_value;
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::spectral_addScalarAll(
		const std::complex<double> &i_value
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
	{
		for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(jj, ii);

			out[2*idx] = t[2*idx] + i_value.real();
			out[2*idx+1] = t[2*idx+1] + i_value.imag();
		}
	}
#else
	CART2D_DATA_COMPLEX_SPECTRAL_FOR_IDX(
			retval.spectral_space_data[idx] = spectral_space_data[idx] + i_value;
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::spectral_div_element_wise(
		const DataSpectral &i_array_data
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_array_data.spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
	{
		for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(jj, ii);

			double ar = t[2*idx];
			double ai = t[2*idx+1];
			double br = in[2*idx];
			double bi = in[2*idx+1];

			double den = br*br+bi*bi;

			if (den == 0)
			{
				out[2*idx] = 0;
				out[2*idx+1] = 0;
			}
			else
			{
				out[2*idx] = (ar*br + ai*bi)/den;
				out[2*idx+1] = (ai*br - ar*bi)/den;
			}
		}
	}

#else
	CART2D_DATA_COMPLEX_SPECTRAL_FOR_IDX(
	{
		if (i_array_data.spectral_space_data[idx] == std::complex<double>(0,0))
			retval.spectral_space_data[idx] = 0;
		else
			retval.spectral_space_data[idx] = spectral_space_data[idx] / i_array_data.spectral_space_data[idx];
	}
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::spectral_invert()	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;


	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
	{
		for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
		{
			int idx = cart2DDataConfig->getArrayIndexByModes_Complex(jj, ii);

			// 1/(a+b*i) = (a-b*i)/(a*a + b*b)

			double re = t[2*idx];
			double im = t[2*idx+1];
			double den = re*re + im*im;

			out[2*idx] = re / den;
			out[2*idx+1] = im / den;
		}
	}
#else
	CART2D_DATA_COMPLEX_SPECTRAL_FOR_IDX(
			retval.spectral_space_data[idx] = 1.0/spectral_space_data[idx];
	);
#endif
	retval.spectral_zeroAliasingModes();

	return retval;
}


void DataSpectral::spectral_print(
		int i_precision
)	const
{
	std::cout << std::setprecision(i_precision);

	/*
	 * WARNING: This follows a different order contrast to how it is stored
	 */
	for (std::size_t m = 0; m < cart2DDataConfig->spectral_complex_data_size[0]; m++)
	{
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_complex_data_size[1]; n++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes_Complex(n, m);
			std::cout << spectral_space_data[idx] << "\t";
		}
		std::cout << std::endl;
	}
}


void DataSpectral::spectral_zeroAliasingModes()	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int k = 0; k < 2; k++)
	{
		if (k == 0)
		{
			/*
			 * First process part between top and bottom spectral data blocks
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
			for (std::size_t jj = cart2DDataConfig->spectral_complex_ranges[0][1][1]; jj < cart2DDataConfig->spectral_complex_ranges[1][1][0]; jj++)
				for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_complex_data_size[0]; ii++)
				{
					spectral_space_data[jj*cart2DDataConfig->spectral_complex_data_size[0]+ii] = 0;
				}
		}
		else
		{
			/*
			 * Then process the aliasing block on the right side
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
			for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_complex_data_size[1]; jj++)
				for (std::size_t ii = cart2DDataConfig->spectral_complex_ranges[0][0][1]; ii < cart2DDataConfig->spectral_complex_ranges[2][0][0]; ii++)
				{
					spectral_space_data[jj*cart2DDataConfig->spectral_complex_data_size[0]+ii] = 0;
				}
		}
	}
}


void DataSpectral::test_realphysical()	const
{

	for (int r = 0; r < 2; r++)
	{
		for (	std::size_t j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0];
				j < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1];
				j++
		) {
			for (	std::size_t i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]+1;
					i < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]-1;
					i++
			) {
				const std::complex<double> &data = spectral_get(j, i);
				std::complex<double> data2;

				if (j == 0)
					data2 = spectral_get(j, cart2DDataConfig->spectral_complex_data_size[0]-i);
				else
					data2 = spectral_get(cart2DDataConfig->spectral_complex_data_size[1]-j, cart2DDataConfig->spectral_complex_data_size[0]-i);

				data2.imag(-data2.imag());

				double error_real = std::abs(data.real() - data2.real());
				double error_imag = std::abs(data.imag() - data2.imag());
				if (error_real > 1e-6 || error_imag > 1e-6)
				{
					std::cout << std::setprecision(5);
					std::cout << "Mode " << j << ", " << i << " : " << error_real << " " << error_imag << std::endl;
					SWEETErrorFatal("Invalid symmetry detected");
				}
			}
		}
	}

}


void DataSpectral::loadCart2DDataGrid(
		const DataGrid &i_cart2DDataGrid
)
{

	cart2DDataConfig->fft_complex_grid_2_spectral_OUTOFPLACE(i_cart2DDataGrid.grid_space_data, spectral_space_data);
}


void DataSpectral::print_spectralData()	const
{
	DataSpectral &rw_array_data = (DataSpectral&)*this;

	//for (std::size_t y = cart2DDataConfig->spectral_complex_data_size[1]-1; y >= 0; y--) // https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
	for (std::size_t y = cart2DDataConfig->spectral_complex_data_size[1]-1; y < cart2DDataConfig->spectral_complex_data_size[1]; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_complex_data_size[0]; x++)
		{
			const std::complex<double> &value = rw_array_data.spectral_get(y, x);
			std::cout << "(" << value.real() << ", " << value.imag() << ")\t";
		}
		std::cout << std::endl;
	}
}


void DataSpectral::print_spectralData_zeroNumZero(
		double i_zero_threshold
)	const
{
	DataSpectral &rw_array_data = (DataSpectral&)*this;

	//for (std::size_t y = cart2DDataConfig->spectral_complex_data_size[1]-1; y >= 0; y--) // https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
	for (std::size_t y = cart2DDataConfig->spectral_complex_data_size[1]-1; y < cart2DDataConfig->spectral_complex_data_size[1]; y--)
	//for (int y = (int)cart2DDataConfig->spectral_complex_data_size[1]-1; y < cart2DDataConfig->spectral_complex_data_size[1]; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_complex_data_size[0]; x++)
		{
			const std::complex<double> &value = rw_array_data.spectral_get(y, x);

			double re = value.real();
			double im = value.imag();

			if (std::abs(re) < i_zero_threshold)	re = 0.0;
			if (std::abs(im) < i_zero_threshold)	im = 0.0;

			std::cout << "(" << re << ", " << im << ")\t";
		}
		std::cout << std::endl;
	}
}

}}}
