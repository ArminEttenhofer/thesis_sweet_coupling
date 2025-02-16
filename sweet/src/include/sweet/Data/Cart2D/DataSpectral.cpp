/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#include "DataSpectral.hpp"
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <complex>
#include <cfloat>
#include <functional>
#include <array>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <utility>

#include <cmath>
#include <sweet/Memory/parmemcpy.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Fatal.hpp>


#if !SWEET_USE_LIBFFT
#error "LIBFFT not activated, but spectral cart2d data compiled in"
#endif

#include <fftw3.h>

#if SWEET_THREADING_SPACE
#	include <omp.h>
#endif

namespace sweet {
namespace Data {
namespace Cart2D {


std::complex<double>& DataSpectral::operator[](std::size_t i)
{
	return spectral_space_data[i];
}

const std::complex<double>& DataSpectral::operator[](std::size_t i)	const
{
	return spectral_space_data[i];
}

void DataSpectral::swap(
		DataSpectral &i_cart2DData
)
{
	SWEET_ASSERT(cart2DDataConfig == i_cart2DData.cart2DDataConfig);

	std::swap(spectral_space_data, i_cart2DData.spectral_space_data);
}



DataSpectral::DataSpectral(
		const Config *i_cart2DDataConfig
)	:
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_cart2DDataConfig != 0);

	setup(i_cart2DDataConfig);
}



DataSpectral::DataSpectral(
		const Config *i_cart2DDataConfig,
		const std::complex<double> &i_value
)	:
	cart2DDataConfig(i_cart2DDataConfig),
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_cart2DDataConfig != 0);

	setup(i_cart2DDataConfig);
	spectral_setValue(i_value);
}




DataSpectral::DataSpectral(
		const Config *i_cart2DDataConfig,
		double &i_value
)	:
	cart2DDataConfig(i_cart2DDataConfig),
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_cart2DDataConfig != 0);

	setup(i_cart2DDataConfig);
	spectral_setValue(i_value);
}



/**
 * Without setup where we need to call setup(...) later on
 */

DataSpectral::DataSpectral()	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
}


DataSpectral::DataSpectral(int i)	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
}



DataSpectral::DataSpectral(
		const DataSpectral &i_cart2d_data
)	:
	cart2DDataConfig(i_cart2d_data.cart2DDataConfig),
	spectral_space_data(nullptr)
{
	if (i_cart2d_data.cart2DDataConfig == nullptr)
		return;

	alloc_data();

	operator=(i_cart2d_data);
}




DataSpectral::DataSpectral(
		DataSpectral &&i_cart2d_data
)	:
	cart2DDataConfig(i_cart2d_data.cart2DDataConfig),
	spectral_space_data(nullptr)
{
	if (i_cart2d_data.cart2DDataConfig == nullptr)
		return;

	std::swap(spectral_space_data, i_cart2d_data.spectral_space_data);
}


DataSpectral::DataSpectral(
		const DataGrid &i_cart2d_data_physical
)	:
	cart2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	setup(i_cart2d_data_physical.cart2DDataConfig);

	loadCart2DDataGrid(i_cart2d_data_physical);
}


DataSpectral::~DataSpectral()
{
	clear();
}



void DataSpectral::clear()
{
	if (spectral_space_data != nullptr)
	{
		sweet::Memory::MemBlockAlloc::free(spectral_space_data, cart2DDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));

		spectral_space_data = nullptr;
		cart2DDataConfig = nullptr;
	}
}



void DataSpectral::_validateRes(
		const Config *i_cart2DDataConfig
)	const
{
	SWEET_ASSERT(cart2DDataConfig->grid_res[0] == i_cart2DDataConfig->grid_res[0]);
	SWEET_ASSERT(cart2DDataConfig->grid_res[1] == i_cart2DDataConfig->grid_res[1]);

	SWEET_ASSERT(cart2DDataConfig->spectral_data_size[0] == i_cart2DDataConfig->spectral_data_size[0]);
	SWEET_ASSERT(cart2DDataConfig->spectral_data_size[1] == i_cart2DDataConfig->spectral_data_size[1]);
}


void DataSpectral::spectral_zeroAliasingModes()
{
//! ALWAYS run this to eliminate Nyquist Frequency even without dealiasing activated
//#if SWEET_USE_CART2D_SPECTRAL_DEALIASING || 1	//! ALWAYS run this to eliminate Nyquist Frequency even without dealiasing activated

	SWEET_THREADING_SPACE_PARALLEL_FOR
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
//#endif
}


void DataSpectral::spectral_debugCheckForZeroAliasingModes()	const
{
#if SWEET_DEBUG

	SWEET_THREADING_SPACE_PARALLEL_FOR
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
					std::complex<double> &data = spectral_space_data[cart2DDataConfig->getArrayIndexByModes(jj, ii)];

					double error = std::sqrt(data.real()*data.real() + data.imag()*data.imag());
					if (error >= 1e-9)
					{
						print_spectralData_zeroNumZero();
						std::cout << "Value at spectral coordinate " << jj << ", " << ii << " should be zero, but is " << data << std::endl;
						SWEETErrorFatal("EXIT");
					}
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
					std::complex<double> &data = spectral_space_data[cart2DDataConfig->getArrayIndexByModes(jj, ii)];

					double error = std::sqrt(data.real()*data.real() + data.imag()*data.imag());
					if (error >= 1e-9)
					{
						print_spectralData_zeroNumZero();
						std::cout << "Value at spectral coordinate " << jj << ", " << ii << " should be zero, but is " << data << std::endl;
						SWEETErrorFatal("EXIT");
					}
				}
		}
	}
#endif
}


DataSpectral& DataSpectral::load_nodealiasing(
		const DataSpectral &i_cart2d_data
)
{
	if (cart2DDataConfig == nullptr)
		SWEETErrorFatal("cart2DDataConfig not initialized");

	sweet::Memory::parmemcpy(spectral_space_data, i_cart2d_data.spectral_space_data, sizeof(Tcomplex)*cart2DDataConfig->spectral_array_data_number_of_elements);

	return *this;
}



DataSpectral& DataSpectral::operator=(
		DataSpectral &&i_cart2d_data
)
{
	if (cart2DDataConfig == nullptr)
		setup(i_cart2d_data.cart2DDataConfig);

	std::swap(spectral_space_data, i_cart2d_data.spectral_space_data);

	return *this;
}


DataSpectral DataSpectral::spectral_returnWithDifferentModes(
		const Config &i_cart2DDataConfig
)	const
{
	return spectral_returnWithDifferentModes(&i_cart2DDataConfig);
}


DataSpectral DataSpectral::spectral_returnWithDifferentModes(
		const Config *i_cart2DDataConfig
)	const
{
	DataSpectral out(i_cart2DDataConfig);

	/*
	 *  0 = invalid
	 * -1 = scale down
	 *  1 = scale up
	 */
	int scaling_mode = 0;

	if (cart2DDataConfig->spectral_modes[0] < out.cart2DDataConfig->spectral_modes[0])
	{
		scaling_mode = 1;
	}
	else if (cart2DDataConfig->spectral_modes[0] > out.cart2DDataConfig->spectral_modes[0])
	{
		scaling_mode = -1;
	}

	if (cart2DDataConfig->spectral_modes[1] < out.cart2DDataConfig->spectral_modes[1])
	{
		SWEET_ASSERT(scaling_mode != -1);
		scaling_mode = 1;
	}
	else if (cart2DDataConfig->spectral_modes[1] > out.cart2DDataConfig->spectral_modes[1])
	{
		SWEET_ASSERT(scaling_mode != 1);
		scaling_mode = -1;
	}

	if (scaling_mode == 0)
	{
		// Just copy the data
		out = *this;
		return out;
	}

	double rescale =
			(double)(out.cart2DDataConfig->grid_number_elements)
			/
			(double)(cart2DDataConfig->grid_number_elements);

	//rescale = 1.0;

	{
		if (scaling_mode == -1)
		{
			/*
			 * more modes -> less modes
			 */

			/*
			 * Region #1
			 *
			 * 00000000 7
			 * 00000000 6
			 * 00000000 5
			 * 00000000 4
			 * 00000000 3
			 * XXXXX000 2
			 * XXXXX000 1
			 * XXXXX000 0
			 */
			{
				const std::size_t* src_range_dim0 = &(cart2DDataConfig->spectral_data_iteration_ranges[0][0][0]);
				const std::size_t* src_range_dim1 = &(cart2DDataConfig->spectral_data_iteration_ranges[0][1][0]);

				const std::size_t* dst_range_dim0 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[0][0][0]);
				const std::size_t* dst_range_dim1 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[0][1][0]);

				SWEET_ASSERT(src_range_dim0[0] == 0);
				SWEET_ASSERT(dst_range_dim0[0] == 0);
				SWEET_ASSERT(src_range_dim1[0] == 0);
				SWEET_ASSERT(dst_range_dim1[0] == 0);

				std::size_t dst_size = dst_range_dim0[1];//-dst_range_dim0[0];

				SWEET_THREADING_SPACE_PARALLEL_FOR
				for (std::size_t j = 0; j < dst_range_dim1[1]; j++)
				{
					std::complex<double> *src = &spectral_space_data[cart2DDataConfig->spectral_data_size[0]*j];
					std::complex<double> *dst = &out.spectral_space_data[out.cart2DDataConfig->spectral_data_size[0]*j];

					for (std::size_t i = 0; i < dst_size; i++)
						dst[i] = src[i]*rescale;
				}
			}


			/*
			 * Region #2
			 *
			 * XXXXX000 7
			 * XXXXX000 6
			 * XXXXX000 5
			 * 00000000 4
			 * 00000000 3
			 * 00000000 2
			 * 00000000 1
			 * 00000000 0
			 */
			{
				const std::size_t* src_range_dim0 = &(cart2DDataConfig->spectral_data_iteration_ranges[1][0][0]);
				const std::size_t* src_range_dim1 = &(cart2DDataConfig->spectral_data_iteration_ranges[1][1][0]);

				const std::size_t* dst_range_dim0 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[1][0][0]);
				const std::size_t* dst_range_dim1 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[1][1][0]);

				SWEET_ASSERT(src_range_dim0[0] == 0);
				SWEET_ASSERT(dst_range_dim0[0] == 0);

				std::size_t dst_size = dst_range_dim0[1];//-dst_range_dim0[0];

				SWEET_THREADING_SPACE_PARALLEL_FOR
				for (std::size_t j = dst_range_dim1[0]; j < dst_range_dim1[1]; j++)
				{
					std::complex<double> *src = &spectral_space_data[cart2DDataConfig->spectral_data_size[0]*(src_range_dim1[1]-(dst_range_dim1[1]-j))];
					std::complex<double> *dst = &out.spectral_space_data[out.cart2DDataConfig->spectral_data_size[0]*j];

					for (std::size_t i = 0; i < dst_size; i++)
						dst[i] = src[i]*rescale;
				}
			}

			out.spectral_zeroAliasingModes();
		}
		else
		{
			/*
			 * less modes -> more modes
			 */

			/*
			 * Region #1
			 *
			 * 00000000 7
			 * 00000000 6
			 * 00000000 5
			 * 00000000 4
			 * 00000000 3
			 * XXXXX000 2
			 * XXXXX000 1
			 * XXXXX000 0
			 */
			out.spectral_setZero();

			{
				int r = 0;

				const std::size_t* src_range_dim0 = &(cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]);
				const std::size_t* src_range_dim1 = &(cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]);

				const std::size_t* dst_range_dim0 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]);
				const std::size_t* dst_range_dim1 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]);

				SWEET_ASSERT(src_range_dim0[0] == 0);
				SWEET_ASSERT(dst_range_dim0[0] == 0);
				SWEET_ASSERT(src_range_dim1[0] == 0);

				std::size_t src_size = src_range_dim0[1];//-dst_range_dim0[0];


				SWEET_THREADING_SPACE_PARALLEL_FOR
				for (std::size_t j = 0; j < src_range_dim1[1]; j++)
				{
					std::complex<double> *src = &spectral_space_data[cart2DDataConfig->spectral_data_size[0]*(j-src_range_dim1[0]+dst_range_dim1[0])];
					std::complex<double> *dst = &out.spectral_space_data[out.cart2DDataConfig->spectral_data_size[0]*j];

					for (std::size_t i = 0; i < src_size; i++)
						dst[i] = src[i]*rescale;
				}
			}


			/*
			 * Region #2
			 *
			 * XXXXX000 7
			 * XXXXX000 6
			 * XXXXX000 5
			 * 00000000 4
			 * 00000000 3
			 * 00000000 2
			 * 00000000 1
			 * 00000000 0
			 */
			{
				int r = 1;

				const std::size_t* src_range_dim0 = &(cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]);
				const std::size_t* src_range_dim1 = &(cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]);

				const std::size_t* dst_range_dim0 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]);
				const std::size_t* dst_range_dim1 = &(out.cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]);

				SWEET_ASSERT(src_range_dim0[0] == 0);
				SWEET_ASSERT(dst_range_dim0[0] == 0);

				std::size_t src_size0 = src_range_dim0[1];//-dst_range_dim0[0];
				std::size_t src_size1 = src_range_dim1[1]-src_range_dim1[0];


				SWEET_THREADING_SPACE_PARALLEL_FOR
				for (std::size_t j = src_range_dim1[0]; j < src_range_dim1[1]; j++)
				{
					std::complex<double> *src = &spectral_space_data[cart2DDataConfig->spectral_data_size[0]*j];
					std::complex<double> *dst = &out.spectral_space_data[out.cart2DDataConfig->spectral_data_size[0]*(dst_range_dim1[1]-src_size1+(j-src_range_dim1[0]))];

					for (std::size_t i = 0; i < src_size0; i++)
					{
#if SWEET_DEBUG
						SWEET_ASSERT((int)(out.spectral_space_data - &dst[i]) < (int)out.cart2DDataConfig->spectral_array_data_number_of_elements);
#endif
						dst[i] = src[i]*rescale;
					}
				}

			}
		}

	}

	return out;
}



DataSpectral DataSpectral::spectral_invert()	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx+=2)
	{
			// 1/(a+b*i) = (a-b*i)/(a*a + b*b)

			double re = t[idx];
			double im = t[idx+1];
			double den = re*re + im*im;

			out[idx] = re / den;
			out[idx+1] = im / den;
	}
#else
	CART2D_DATA_SPECTRAL_FOR_IDX(
			retval.spectral_space_data[idx] = 1.0/spectral_space_data[idx];
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}



void DataSpectral::loadCart2DDataGrid(
	const DataGrid &i_cart2DDataGrid
)
{
	cart2DDataConfig->fft_grid_2_spectral_OUTOFPLACE(i_cart2DDataGrid.grid_space_data, spectral_space_data);
	// ALWAYS zero aliasing modes after doing transformation to spectral space
	spectral_zeroAliasingModes();
}


DataGrid DataSpectral::getCart2DDataGrid()	const
{
	/**
	 * Warning: This is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 */
	DataSpectral tmp(*this);
	DataGrid retval(cart2DDataConfig);

	cart2DDataConfig->fft_spectral_2_grid_INPLACE(tmp.spectral_space_data, retval.grid_space_data);
	return retval;
}


DataGrid DataSpectral::toGrid()	const
{
	/*
	 * Warning: This is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 */
	DataSpectral tmp(*this);
	DataGrid retval(cart2DDataConfig);
	cart2DDataConfig->fft_spectral_2_grid_INPLACE(tmp.spectral_space_data, retval.grid_space_data);

	return retval;
}


Cart2DComplex::DataGrid DataSpectral::getCart2DDataGridComplex()	const
{
	Cart2DComplex::DataGrid out(cart2DDataConfig);

	/*
	 * WARNING:
	 * We have to use a temporary array here because of destructive FFTW transformations
	 */
	DataSpectral tmp_spectral(*this);
	DataGrid tmp_physical(cart2DDataConfig);
	cart2DDataConfig->fft_spectral_2_grid_INPLACE(tmp_spectral.spectral_space_data, tmp_physical.grid_space_data);

	sweet::Memory::parmemcpy(out.grid_space_data, tmp_physical.grid_space_data, sizeof(double)*cart2DDataConfig->grid_number_elements);

	return out;
}



void DataSpectral::_main_setup(
	const Config *i_cart2DDataConfig
)
{
	SWEET_ASSERT(cart2DDataConfig == nullptr);

	cart2DDataConfig = i_cart2DDataConfig;
	alloc_data();
}


void DataSpectral::setup(
	const Config *i_cart2DDataConfig
)
{
	_main_setup(i_cart2DDataConfig);
}


void DataSpectral::setup(
	const Config &i_cart2DDataConfig
)
{
	_main_setup(&i_cart2DDataConfig);
}


void DataSpectral::setup(
	const Config *i_cart2DDataConfig,
	double i_value
)
{
	_main_setup(i_cart2DDataConfig);

	spectral_setValue(i_value);
}


void DataSpectral::alloc_data()
{
	SWEET_ASSERT(spectral_space_data == nullptr);
	spectral_space_data = sweet::Memory::MemBlockAlloc::alloc<Tcomplex>(cart2DDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));
}


void DataSpectral::setup_if_required(
	const Config *i_cart2DDataConfig
)
{
	if (cart2DDataConfig != nullptr)
		return;

	_main_setup(i_cart2DDataConfig);
}


DataSpectral &DataSpectral::operator=(int i_value)
{
	spectral_setValue(std::complex<double>(i_value, 0));

	return *this;
}


DataSpectral &DataSpectral::operator=(double i_value)
{
	spectral_setValue(i_value);

	return *this;
}


DataSpectral& DataSpectral::operator=(
		const DataSpectral &i_cart2d_data
)
{
	if (i_cart2d_data.cart2DDataConfig == nullptr)
		return *this;

	if (cart2DDataConfig == nullptr)
		setup(i_cart2d_data.cart2DDataConfig);

	sweet::Memory::parmemcpy(spectral_space_data, i_cart2d_data.spectral_space_data, sizeof(Tcomplex)*cart2DDataConfig->spectral_array_data_number_of_elements);

	return *this;
}


DataSpectral DataSpectral::operator+(
		double i_value
)	const
{
	DataSpectral out(*this);

	out.spectral_space_data[0] += i_value * (double)cart2DDataConfig->grid_number_elements;
	out.spectral_zeroAliasingModes();

	return out;
}


DataSpectral DataSpectral::operator+(
		const DataSpectral &i_cart2d_data
)	const
{
	_validateRes(i_cart2d_data.cart2DDataConfig);

	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_cart2d_data.spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = t[idx] + in[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = spectral_space_data[idx] + i_cart2d_data.spectral_space_data[idx];
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


const DataSpectral& DataSpectral::operator+=(
		double i_value
)	const
{
	spectral_space_data[0] += i_value * (double)cart2DDataConfig->grid_number_elements;

	return *this;
}


DataSpectral& DataSpectral::operator+=(
		const DataSpectral &i_cart2d_data
)
{
	_validateRes(i_cart2d_data.cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_cart2d_data.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		t[idx] += in[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] += i_cart2d_data.spectral_space_data[idx];

#endif
	spectral_zeroAliasingModes();

	return *this;
}


DataSpectral DataSpectral::operator-(
		double i_value
)	const
{
	DataSpectral out(*this);
	out.spectral_space_data[0] -= i_value * (double)cart2DDataConfig->grid_number_elements;

	out.spectral_zeroAliasingModes();

	return out;
}


DataSpectral DataSpectral::operator-(
		const DataSpectral &i_cart2d_data
)	const
{
	_validateRes(i_cart2d_data.cart2DDataConfig);

	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_cart2d_data.spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = t[idx] - in[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = spectral_space_data[idx] - i_cart2d_data.spectral_space_data[idx];
#endif
	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral& DataSpectral::operator-=(
		const DataSpectral &i_cart2d_data
)
{
	_validateRes(i_cart2d_data.cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_cart2d_data.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		t[idx] -= in[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] -= i_cart2d_data.spectral_space_data[idx];
#endif

	spectral_zeroAliasingModes();

	return *this;
}


DataSpectral& DataSpectral::operator-=(
		double i_value
)
{
	spectral_space_data[0] -= i_value * (double)cart2DDataConfig->grid_number_elements;

	spectral_zeroAliasingModes();

	return *this;
}


DataSpectral DataSpectral::operator-(
		const DataGrid &i_cart2d_data_physical
)	const
{
	_validateRes(i_cart2d_data_physical.cart2DDataConfig);

	DataSpectral i_cart2d_data_spectral(cart2DDataConfig);
	DataSpectral retval(cart2DDataConfig);

	i_cart2d_data_spectral.loadCart2DDataGrid(i_cart2d_data_physical);

#if 1
	double *t = (double*)spectral_space_data;
	double *in = (double*)i_cart2d_data_spectral.spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = t[idx] - in[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = spectral_space_data[idx] - i_cart2d_data_spectral.spectral_space_data[idx];
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::operator-()	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = -t[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = -spectral_space_data[idx];
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::operator*(
		double i_value
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = t[idx]*i_value;
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = spectral_space_data[idx]*i_value;
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::operator*(
		const DataSpectral &i_cart2d_data
)	const
{
	_validateRes(i_cart2d_data.cart2DDataConfig);

	DataGrid a = toGrid();
	DataGrid b = i_cart2d_data.toGrid();

	DataGrid mul(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		mul.grid_space_data[i] = a.grid_space_data[i]*b.grid_space_data[i];
	}

	// Dealiasing is performed inside the following call
	DataSpectral out(mul);

	return out;
}


const DataSpectral& DataSpectral::operator*=(
		double i_value
)	const
{
#if 1
	double *t = (double*)spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		t[idx] *= i_value;
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] *= i_value;
#endif

	return *this;
}


const DataSpectral& DataSpectral::operator*=(
		const std::complex<double> &i_value
)	const
{
#if 1
	double *t = (double*)spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx+=2)
	{
		double re = t[idx];
		double im = t[idx+1];
		t[idx] = re*i_value.real() - im*i_value.imag();
		t[idx+1] = re*i_value.imag() + im*i_value.real();
	}
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] *= i_value;
#endif

	return *this;
}


const DataSpectral& DataSpectral::operator/=(
		double i_value
)	const
{
#if 1
	double *t = (double*)spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		t[idx] /= i_value;
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] /= i_value;
#endif

	return *this;
}


DataSpectral DataSpectral::operator/(
		const DataSpectral &i_cart2d_data
)	const
{
	_validateRes(i_cart2d_data.cart2DDataConfig);

	DataGrid a = toGrid();
	DataGrid b = i_cart2d_data.toGrid();

	DataGrid div(cart2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
		div.grid_space_data[i] = a.grid_space_data[i]/b.grid_space_data[i];

	DataSpectral out(div);

	return out;
}


DataGrid DataSpectral::multiplication_grid_space(
			const DataGrid &i_a,
			const DataGrid &i_b
) const
{
	_validateRes(i_a.cart2DDataConfig);
	_validateRes(i_b.cart2DDataConfig);

	DataGrid mul(cart2DDataConfig);
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->grid_number_elements; i++)
	{
		mul.grid_space_data[i] = i_a.grid_space_data[i]*i_b.grid_space_data[i];
	}

	DataSpectral out_spec(mul);
	return out_spec.toGrid();
}


DataSpectral DataSpectral::operator/(
		double i_value
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = t[idx]/i_value;
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = spectral_space_data[idx]/i_value;
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::operator_scalar_sub_this(
		double i_value
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out[idx] = -t[idx];
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
		retval.spectral_space_data[idx] = -spectral_space_data[idx];
#endif

	retval.spectral_space_data[0] = i_value*std::sqrt(4.0*M_PI) + retval.spectral_space_data[0];
	return retval;
}


DataSpectral DataSpectral::spectral_addScalarAll(
		const double &i_value
)	const
{
	DataSpectral retval(cart2DDataConfig);

#if 1
	double *t = (double*)spectral_space_data;
	double *out = (double*)retval.spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx+=2)
	{
			out[idx] = t[idx] + i_value;
			out[idx+1] = t[idx+1];
	}
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx++)
			retval.spectral_space_data[idx] = spectral_space_data[idx] + i_value;
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


	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx+=2)
	{
			double ar = t[idx];
			double ai = t[idx+1];
			double br = in[idx];
			double bi = in[idx+1];

			double den = br*br+bi*bi;

			if (den == 0)
			{
				out[idx] = 0;
				out[idx+1] = 0;
			}
			else
			{
				out[idx] = (ar*br + ai*bi)/den;
				out[idx+1] = (ai*br - ar*bi)/den;
			}
	}
#else
	CART2D_DATA_SPECTRAL_FOR_IDX(

			double ar = spectral_space_data[idx].real();
			double ai = spectral_space_data[idx].imag();
			double br = i_array_data.spectral_space_data[idx].real();
			double bi = i_array_data.spectral_space_data[idx].imag();

			double den = br*br+bi*bi;

			if (den == 0)
			{
				retval.spectral_space_data[idx].real(0);
				retval.spectral_space_data[idx].imag(0);
			}
			else
			{
				retval.spectral_space_data[idx].real((ar*br + ai*bi)/den);
				retval.spectral_space_data[idx].imag((ai*br - ar*bi)/den);
			}
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::spectral_solve_helmholtz(
		const double &i_a,
		const double &i_b,
		double r
)
{
	DataSpectral out(*this);

	const double a = i_a;
	const double b = i_b;

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			io_data /= (a + (-b* (m * m + n * n) )); //???
		}
	);

	return out;
}


DataSpectral DataSpectral::spectral_solve_laplace(
		double r
)
{
	DataSpectral out(*this);

	const double b = 1.0;

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			if (n == 0)
				io_data = 0;
			else
				io_data /= (-b*(m * m + n * n));  // ????
		}
	);

	return out;
}


const DataSpectral& DataSpectral::spectral_truncate()	const
{
	DataGrid tmp(cart2DDataConfig);

	DataSpectral tmps(*this);
	cart2DDataConfig->fft_spectral_2_grid_INPLACE(tmps.spectral_space_data, tmp.grid_space_data);
	cart2DDataConfig->fft_grid_2_spectral_OUTOFPLACE(tmp.grid_space_data, spectral_space_data);

	return *this;
}



void DataSpectral::spectral_update_lambda(
		std::function<void(int,int,Tcomplex&)> i_lambda
)
{
#if 1

	/*
	 * This version updates all elements in spectral space, also the ones who would be affected by aliasing
	 */
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int r = 0; r < 2; r++)								\
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
		for (std::size_t jj = 0; jj < cart2DDataConfig->spectral_data_size[1]; jj++)		\
		{
			for (std::size_t ii = 0; ii < cart2DDataConfig->spectral_data_size[0]; ii++)		\
			{
				std::size_t idx = jj*cart2DDataConfig->spectral_data_size[0]+ii;
				i_lambda(jj, ii, spectral_space_data[idx]);
			}			\
		}				\
	}

#elif 1

	/*
	 * This version only updates the elements in spectral space which are not related to aliased modes
	 */
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int r = 0; r < 2; r++)								\
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
		for (std::size_t jj = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
		{
			for (std::size_t ii = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
			{
				std::size_t idx = jj*cart2DDataConfig->spectral_data_size[0]+ii;
				i_lambda(jj, ii, spectral_space_data[idx]);
			}			\
		}				\
	}

#else

	CART2D_DATA_SPECTRAL_FOR_IDX(
				i_lambda(jj, ii, spectral_space_data[idx]);
			)
#endif
}


const std::complex<double>& DataSpectral::spectral_get(
		int i_n,
		int i_m
)	const
{
	SWEET_ASSERT(i_n >= 0 && i_m >= 0);
	SWEET_ASSERT(i_n <= (int)cart2DDataConfig->spectral_data_size[1]);
	SWEET_ASSERT(i_m <= (int)cart2DDataConfig->spectral_data_size[0]);
	///SWEET_ASSERT(i_m <= i_n);

	return spectral_space_data[cart2DDataConfig->getArrayIndexByModes(i_n, i_m)];
}


void DataSpectral::spectral_set(
		int i_n,
		int i_m,
		const std::complex<double> &i_data
)	const
{
#if SWEET_DEBUG
	if (i_n < 0 ||  i_m < 0)
		SWEETErrorFatal("Out of boundary a");

	if (i_m > (int)cart2DDataConfig->spectral_data_size[0])
		SWEETErrorFatal("Out of boundary b");

	if (i_n > (int)cart2DDataConfig->spectral_data_size[1])
		SWEETErrorFatal("Out of boundary c");
#endif

	spectral_space_data[cart2DDataConfig->getArrayIndexByModes(i_n, i_m)] = i_data;
}


void DataSpectral::spectral_setZero()
{

#if 1
	double *t = (double*)spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->spectral_array_data_number_of_elements; i++)
		t[i] = 0;
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->spectral_array_data_number_of_elements; i++)
		spectral_space_data[i] = {0,0};
#endif
}


void DataSpectral::spectral_setValue(
		const std::complex<double> &i_value
)
{
#if 1
	double *t = (double*)spectral_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*cart2DDataConfig->spectral_array_data_number_of_elements; i+=2)
	{
		t[i] = i_value.real();
		t[i+1] = i_value.imag();
	}
#else
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < cart2DDataConfig->spectral_array_data_number_of_elements; i++)
		spectral_space_data[i] = i_value;
#endif
}


void DataSpectral::spectral_add_grid_constant(
		double i_value
)	const
{
	spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
}


std::complex<double> DataSpectral::spectral_reduce_sum_quad_increasing()	const
{
	std::complex<double> sum = 0;
	std::complex<double> c = 0;
#if SWEET_THREADING_SPACE
//#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
#endif
	for (std::size_t i = 0; i < cart2DDataConfig->spectral_array_data_number_of_elements; i++)
	{
		std::complex<double> value = spectral_space_data[i]*(double)i;

		// Use Kahan summation
		std::complex<double> y = value - c;
		std::complex<double> t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	return sum;
}


std::complex<double> DataSpectral::spectral_reduce_sum_quad()	const
{
	std::complex<double> sum = 0;
	std::complex<double> c = 0;
#if SWEET_THREADING_SPACE
//#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
#endif
	for (std::size_t i = 0; i < cart2DDataConfig->spectral_array_data_number_of_elements; i++)
	{
		std::complex<double> value = spectral_space_data[i];

		// Use Kahan summation
		std::complex<double> y = value - c;
		std::complex<double> t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	sum -= c;

	return sum;
}


double DataSpectral::spectral_reduce_sum_sqr_quad()	const
{
	std::complex<double> sum = 0;

	for (std::size_t n = 0; n <= cart2DDataConfig->spectral_data_size[1]; n++)
		for (std::size_t m = 0; m <= cart2DDataConfig->spectral_data_size[0]; m++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes(n, m);
			sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
		}

	return sum.real();
}


double DataSpectral::spectral_reduce_sum_sq()
{
	double rms = 0.0;

	std::complex<double> sum = 0;

	for (int r = 0; r < 2; r++)
	{
		for (std::size_t jj = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
		{
			for (std::size_t ii = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
			{	
				std::size_t idx = jj*cart2DDataConfig->spectral_data_size[0]+ii;
				sum += (spectral_space_data[idx]*std::conj(spectral_space_data[idx]));
			}
		}
	}

	//sum = std::__complex_sqrt (sum/(double)(cart2DDataConfig->spectral_array_data_number_of_elements));

	if (sum.imag()>DBL_EPSILON)
		SWEETErrorFatal("Reduce operation of complex values (rms) error");

	rms = sum.real()/(double)(cart2DDataConfig->spectral_array_data_number_of_elements);
	//rms = (sum.real()*sum.real()+sum.imag()*sum.imag());
	return rms;

}


std::complex<double> DataSpectral::spectral_reduce_min()	const
{
	std::complex<double> error = std::numeric_limits<double>::infinity();

	for (std::size_t j = 0; j < cart2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		error.real(std::min(spectral_space_data[j].real(), error.real()));
		error.imag(std::min(spectral_space_data[j].imag(), error.imag()));
	}

	return error;
}


std::complex<double> DataSpectral::spectral_reduce_max()	const
{
	std::complex<double> error = -std::numeric_limits<double>::infinity();

	for (std::size_t j = 0; j < cart2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		error.real(std::max(spectral_space_data[j].real(), error.real()));
		error.imag(std::max(spectral_space_data[j].imag(), error.imag()));
	}

	return error;
}


double DataSpectral::spectral_reduce_max_abs()	const
{
	double error = -std::numeric_limits<double>::infinity();
	std::complex<double> w = {0,0};
	for (std::size_t j = 0; j < cart2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
		error = std::max(std::abs(w), error);
	}

	return error;
}


double DataSpectral::spectral_reduce_max_abs(std::size_t rnorm)	const
{

	assert (rnorm <= cart2DDataConfig->spectral_data_size[0]);
	assert (rnorm <= cart2DDataConfig->spectral_data_size[1]);

	double error = -std::numeric_limits<double>::infinity();
	std::complex<double> w = {0,0};

	for (std::size_t n = 0; n <= rnorm; n++)
		for (std::size_t m = 0; m <= rnorm; m++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes(n, m);
			w = spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
			error = std::max(std::abs(w), error);
		}

	return error;
}


double DataSpectral::spectral_reduce_min_abs()	const
{
	double error = std::numeric_limits<double>::infinity();
	std::complex<double> w = {0,0};
	for (std::size_t j = 0; j < cart2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
		error = std::min(std::abs(w), error);
	}

	return error;
}


double DataSpectral::spectral_reduce_rms()
{
	double rms = 0.0;

	std::complex<double> sum = 0;

	for (int r = 0; r < 2; r++)
	{
		for (std::size_t jj = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
		{
			for (std::size_t ii = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
			{	
				std::size_t idx = jj*cart2DDataConfig->spectral_data_size[0]+ii;
				sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
			}
		}
	}

	//sum = std::__complex_sqrt (sum/(double)(cart2DDataConfig->spectral_array_data_number_of_elements));

	if (sum.imag()>DBL_EPSILON)
		SWEETErrorFatal("Reduce operation of complex values (rms) error");

	rms = std::sqrt(sum.real()/(double)(cart2DDataConfig->spectral_array_data_number_of_elements));
	return rms;

}


bool DataSpectral::spectral_reduce_is_any_nan_or_inf()	const
{
	bool retval = false;

#if SWEET_THREADING_SPACE
	#pragma omp parallel for simd reduction(|:retval)
#endif
	for (std::size_t j = 0; j < cart2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		retval |= std::isnan(spectral_space_data[j].real());
		retval |= std::isinf(spectral_space_data[j].real());
		retval |= std::isnan(spectral_space_data[j].imag());
		retval |= std::isinf(spectral_space_data[j].imag());
	}

	return retval;
}


bool DataSpectral::spectral_is_first_nan_or_inf()	const
{
	bool retval = false;

	retval |= std::isnan(spectral_space_data[0].real());
	retval |= std::isinf(spectral_space_data[0].real());
	retval |= std::isnan(spectral_space_data[0].imag());
	retval |= std::isinf(spectral_space_data[0].imag());

	return retval;
}


void DataSpectral::spectral_print(
		int i_precision,
		double i_abs_threshold
)	const
{
	std::cout << std::setprecision(i_precision);

	std::cout << "m \\ n ----->" << std::endl;
	for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]; m++)
	{
		///std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]; n++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes(n, m);
			if (std::abs(spectral_space_data[idx]) < i_abs_threshold)
				std::cout << 0 << "\t";
			else
				std::cout << spectral_space_data[idx] << "\t";
			idx++;
		}
		std::cout << std::endl;
	}
}


void DataSpectral::print_spectralData()	const
{
	DataSpectral &rw_array_data = (DataSpectral&)*this;

	//for (std::size_t y = cart2DDataConfig->spectral_data_size[1]-1; y >= 0; y--) // https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
	for (std::size_t y = cart2DDataConfig->spectral_data_size[1]-1; y < cart2DDataConfig->spectral_data_size[1]; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
		{
			const std::complex<double> &value = rw_array_data.spectral_get(y, x);
			std::cout << "(" << value.real() << ", " << value.imag() << ")\t";
		}
		std::cout << std::endl;
	}
}



void DataSpectral::print_spectralData_zeroNumZero(double i_zero_threshold)	const
{
	DataSpectral &rw_array_data = (DataSpectral&)*this;

	//for (std::size_t y = cart2DDataConfig->spectral_data_size[1]-1; y >= 0; y--) //  https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
	for (std::size_t y = cart2DDataConfig->spectral_data_size[1]-1; y < cart2DDataConfig->spectral_data_size[1]; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
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


void DataSpectral::print_spectralIndex()	const
{
	DataSpectral &rw_array_data = (DataSpectral&)*this;

	//for (std::size_t y = cart2DDataConfig->spectral_data_size[1]-1; y >= 0; y--) //  https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
	for (int y = cart2DDataConfig->spectral_data_size[1]-1; y >= 0; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
		{
			const std::complex<double> &value = rw_array_data.spectral_get(y, x);
			std::cout << "(" << x << ", " << y << ", " << value.real() << ", " << value.imag() << ")\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

}


void DataSpectral::print_spectralNonZero()	const
{
	DataSpectral &rw_array_data = (DataSpectral&)*this;

	for (int y = cart2DDataConfig->spectral_data_size[1]-1; y >= 0; y--)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
		{
			const std::complex<double> &value = rw_array_data.spectral_get(y, x);
			if (value.real()*value.real()+value.imag()*value.imag() > 1.0e-13)
				std::cout << "(" << x << ", " << y << ", " << value.real() << ", " << value.imag() << ")" <<std::endl;;
		}
	}
}


void DataSpectral::spectral_structure_print(
		int i_precision,
		double i_abs_threshold
)	const
{
	std::cout << std::setprecision(i_precision);

	std::cout << "m \\ n ----->" << std::endl;
	for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]; m++)
	{
		//std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]; n++)
		{
			std::cout << "(" << m << "," << n << ")" << "\t";
		}
		std::cout << std::endl;
	}
}


void DataSpectral::spectrum_file_write(
		const std::string &i_filename,
		const char *i_title,
		int i_precision
)	const
{
	std::ofstream file(i_filename, std::ios_base::trunc);

	file << std::setprecision(i_precision);
	file << "#SWEET_CART2D_SPECTRAL_DATA_ASCII" << std::endl;

	file << "#TI " << i_title << std::endl;

	// Use 0 to make it processable by python
	file << "0\t";

	std::complex<double> w = {0,0};
	std::vector<double> sum(cart2DDataConfig->spectral_data_size[1]+1,0);
	std::vector<double> sum_squared(cart2DDataConfig->spectral_data_size[1]+1,0);
	std::vector<double> max(cart2DDataConfig->spectral_data_size[1]+1,0);

	for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]; m++)
	{
		///std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]; n++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes(n, m);
			w = spectral_space_data[idx];
			idx++;

			sum[n]         += std::abs(w);
			sum_squared[n] += std::abs(w * w);
			if (std::abs(w) >= max[n]) max[n] = std::abs(w);
		}
	}

	file << cart2DDataConfig->spectral_data_size[0] << " "
			<< cart2DDataConfig->spectral_data_size[1] << std::endl;
	for (std::size_t n = 0; n <= cart2DDataConfig->spectral_data_size[1]; n++)
		file << n << " " << sum[n] << " " << max[n] << " " << std::sqrt(sum_squared[n]) << std::endl;

	file.close();
}


void DataSpectral::spectrum_abs_file_write_line(
		const std::string &i_filename,
		const char *i_title,
		const double i_time,
		int i_precision,
		double i_abs_threshold,
		int i_reduce_mode_factor
)	const
{

	std::ofstream file;

	if (i_time == 0.0){
		file.open(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);
		file << "#SWEET_CART2D_SPECTRAL_ABS_EVOL_ASCII" << std::endl;
		file << "#TI " << i_title << std::endl;
		file << "0\t" << std::endl; // Use 0 to make it processable by python
		file << "(n_max=" <<cart2DDataConfig->spectral_data_size[1] << " m_max="
				<< cart2DDataConfig->spectral_data_size[1] << ")" << std::endl;
		file << "timestamp\t" ;
		for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
		{
			//std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
			for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
			{
				file << "(" << n << ";" << m << ")\t" ;
			}
		}
		file<< "SpectralSum" <<std::endl;
	}
	else{
		file.open(i_filename, std::ios_base::app);
		file << std::setprecision(i_precision);
	}

	std::complex<double> w = {0,0};
	double wabs = 0.0;
	double sum = 0.0;
	//std::cout << "n" << " " << "m" << " " << "norm" <<std::endl;
	file << i_time << "\t";
	for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
	{
		///std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes(n, m);
			w = spectral_space_data[idx];
			wabs = std::abs(w * std::conj(w));
			if ( m > 0 ) sum += 2*wabs; //term appears twice in the spectrum
			else sum += wabs;      // term appears only once

			if ( wabs < i_abs_threshold){
				//file << "(" << n << "," << m << ")\t" <<std::endl;
				file <<  0 << "\t"; //<<std::endl;
			}
			else{
				//file << "(" << n << "," << m << ")\t" <<std::endl;
				file <<  wabs << "\t"; //<<std::endl;;
				//std::cout << n << " " << m << " " << wabs <<std::endl;
			}
			idx++;
		}
	}
	file<< sum << std::endl;
	file.close();
}


void DataSpectral::spectrum_phase_file_write_line(
		const std::string &i_filename,
		const char *i_title,
		const double i_time,
		int i_precision,
		double i_abs_threshold,
		int i_reduce_mode_factor
)	const
{

	std::ofstream file;

	if (i_time == 0.0){
		file.open(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);
		file << "#SWEET_CART2D_SPECTRAL_PHASE_EVOL_ASCII" << std::endl;
		file << "#TI " << i_title << std::endl;
		file << "0\t" << std::endl; // Use 0 to make it processable by python
		file << "(n_max=" <<cart2DDataConfig->spectral_data_size[1] << " m_max="
				<< cart2DDataConfig->spectral_data_size[1] << ")" << std::endl;
		file << "timestamp\t" ;
		for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
		{
			//std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
			for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
			{
				file << "(" << n << ";" << m << ")\t" ;
			}
		}
		file << std::endl;
	}
	else{
		file.open(i_filename, std::ios_base::app);
		file << std::setprecision(i_precision);
	}

	std::complex<double> w = {0,0};
	double wphase = 0.0;

	file << i_time << "\t";
	for (std::size_t m = 0; m < cart2DDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
	{
		///std::size_t idx = cart2DDataConfig->getArrayIndexByModes(m, m);
		for (std::size_t n = 0; n < cart2DDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
		{
			std::size_t idx = cart2DDataConfig->getArrayIndexByModes(n, m);
			w = spectral_space_data[idx];
			wphase = std::arg(w); // std::abs(w * std::conj(w));

			file <<  wphase << "\t"; //<<std::endl;;

			idx++;
		}
	}
	file<< std::endl;
	file.close();
}


void DataSpectral::file_write_binary_spectral(
		const std::string &i_filename
)	const
{
	std::ofstream file(i_filename, std::ios_base::trunc | std::ios_base::binary);

	if (!file.is_open())
		SWEETErrorFatal("Error while opening file");

	file << "SWEET" << std::endl;
	file << "DATA_TYPE MODES_DATA" << std::endl;
	file << "MODES_M_MAX " << cart2DDataConfig->spectral_data_size[0] << std::endl;
	file << "MODES_N_MAX " << cart2DDataConfig->spectral_data_size[1] << std::endl;
	file << "GRID_TYPE GAUSSIAN" << std::endl;
	file << "NUM_ELEMENTS " << cart2DDataConfig->spectral_array_data_number_of_elements << std::endl;
	file << "FIN" << std::endl;

	file.write((const char*)spectral_space_data, sizeof(std::complex<double>)*cart2DDataConfig->spectral_array_data_number_of_elements);

	file.close();
}


void DataSpectral::loadDataFromFile(
		const std::string i_filename,
		const std::string i_filemode
)
{
	if (i_filemode == "csv")
		file_read_csv_grid(i_filename);
	else if (i_filemode == "bin")
		file_read_binary_spectral(i_filename);
	else
		SWEETErrorFatal("Unknown output file mode '"+i_filemode+"'");
}

void DataSpectral::file_read_csv_grid(
		const std::string &i_filename
)
{
	sweet::Data::Cart2D::DataGrid tmp(cart2DDataConfig);
	tmp.file_read_csv_grid(i_filename.c_str());
	loadCart2DDataGrid(tmp);
}

void DataSpectral::file_read_binary_spectral(
		const std::string &i_filename
)
{
	std::ifstream file(i_filename, std::ios_base::binary);

	if (!file.is_open())
		SWEETErrorFatal("Error while opening file");

	std::string magic;
	std::getline(file, magic);

	if (magic != "SWEET")
		SWEETErrorFatal("Magic code 'SWEET' not found");

	std::string data_type;
	int num_x = -1;
	int num_y = -1;
	int size = -1;

	while (true)
	{
		std::string buf;
		file >> buf;

		if (buf == "FIN")
			break;

		if (buf == "DATA_TYPE")
		{
			// load data type
			file >> data_type;
			std::cout << data_type << std::endl;
			continue;
		}

		if (buf == "NUM_X")
		{
			file >> buf;
			num_x = std::stoi(buf);
			std::cout << num_x << std::endl;
			continue;
		}

		if (buf == "NUM_Y")
		{
			file >> buf;
			num_y = std::stoi(buf);
			std::cout << num_y << std::endl;
			continue;
		}

		if (buf == "SIZE")
		{
			file >> buf;
			size = std::stoi(buf);
			std::cout << size << std::endl;
			continue;
		}

		SWEETErrorFatal("Unknown Tag '"+buf+"'");
	}

	// read last newline
	char nl;
	file.read(&nl, 1);
	std::cout << file.tellg() << std::endl;

	if (data_type != "MODES_DATA")
		SWEETErrorFatal("Unknown data type '"+data_type+"'");

	if (num_x != (int)cart2DDataConfig->spectral_data_size[0])
		SWEETErrorFatal("NUM_X "+std::to_string(num_x)+" doesn't match cart2DDataConfig");

	if (num_y != (int)cart2DDataConfig->spectral_data_size[1])
		SWEETErrorFatal("NUM_Y "+std::to_string(num_y)+" doesn't match cart2DDataConfig");

	file.read((char*)spectral_space_data, sizeof(std::complex<double>)*cart2DDataConfig->spectral_array_data_number_of_elements);

	file.close();
}


bool DataSpectral::file_spectral_saveData_ascii(
		const char *i_filename,
		char i_separator,
		int i_precision,
		int dimension
)	const
{

	std::ofstream file(i_filename, std::ios_base::trunc);
	file << std::setprecision(i_precision);

	file << "#SWEET_CART2D_SPECTRAL_CPLX_DATA_ASCII" << std::endl;

	size_t ymax = 0;
	if (dimension == 2)
		ymax = cart2DDataConfig->spectral_data_size[1];
	else
		ymax = 1;

	for (std::size_t y = 0; y < ymax; y++)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
		{
			const std::complex<double> &value = spectral_get(y, x);
			file << "(" << value.real() << ", " << value.imag() << ")";

			if (x < cart2DDataConfig->spectral_data_size[0]-1)
				file << i_separator;
			else
				file << std::endl;
		}
	}

	return true;
}


bool DataSpectral::file_spectral_abs_saveData_ascii(
		const char *i_filename,
		char i_separator,
		int i_precision,
		int dimension
)	const
{

	std::ofstream file(i_filename, std::ios_base::trunc);
	file << std::setprecision(i_precision);

	file << "#SWEET_CART2D_SPECTRAL_ABS_DATA_ASCII" << std::endl;

	size_t ymax = 0;
	if (dimension == 2)
		ymax = cart2DDataConfig->spectral_data_size[1];
	else
		ymax = 1;

	for (std::size_t y = 0; y < ymax; y++)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
		{
			file << spectral_return_amplitude(y, x);

			if (x < cart2DDataConfig->spectral_data_size[0]-1)
				file << i_separator;
			else
				file << std::endl;
		}
	}

	return true;
}


bool DataSpectral::file_spectral_arg_saveData_ascii(
		const char *i_filename,
		char i_separator,
		int i_precision,
		int dimension
)	const
{

	std::ofstream file(i_filename, std::ios_base::trunc);
	file << std::setprecision(i_precision);

	size_t ymax = 0;
	if (dimension == 2)
		ymax = cart2DDataConfig->spectral_data_size[1];
	else
		ymax = 1;

	for (std::size_t y = 0; y < ymax; y++)
	{
		for (std::size_t x = 0; x < cart2DDataConfig->spectral_data_size[0]; x++)
		{

			file << spectral_return_phase(y, x);

			if (x < cart2DDataConfig->spectral_data_size[0]-1)
				file << i_separator;
			else
				file << std::endl;
		}
	}

	return true;
}


double DataSpectral::get_average()	const
{
	return spectral_space_data[0].real()/(double)cart2DDataConfig->grid_number_elements;
}


double DataSpectral::spectral_return_amplitude(
		std::size_t j,
		std::size_t i
) const
{
	std::complex<double> val = spectral_get(j,i);
	val.real(val.real()*2/(cart2DDataConfig->grid_number_elements));
	val.imag(val.imag()*2/(cart2DDataConfig->grid_number_elements));

	return std::abs(val);
}


double DataSpectral::spectral_return_phase(
		std::size_t j,
		std::size_t i
) const
{
	std::complex<double> val = spectral_get(j,i);
	val.real(val.real()*2/(cart2DDataConfig->grid_number_elements));
	val.imag(val.imag()*2/(cart2DDataConfig->grid_number_elements));

	return std::arg(val);
}


void DataSpectral::normalize(
		const std::string &normalization
)
{
	if (normalization == "avg_zero")
	{
		// move average value to 0
		double phi_min = getCart2DDataGrid().grid_reduce_min();
		double phi_max = getCart2DDataGrid().grid_reduce_max();

		double avg = 0.5*(phi_max+phi_min);

		operator-=(avg);
	}
	else if (normalization == "min_zero")
	{
		// move minimum value to zero
		double phi_min = getCart2DDataGrid().grid_reduce_min();
		operator-=(phi_min);
	}
	else if (normalization == "max_zero")
	{
		// move maximum value to zero
		double phi_max = getCart2DDataGrid().grid_reduce_max();
		operator-=(phi_max);
	}
	else if (normalization == "")
	{
	}
	else
	{
		SWEETErrorFatal("Normalization not supported!");
	}
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

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (size_t idx = 0; idx < 2*cart2DDataConfig->spectral_array_data_number_of_elements; idx+=2)
	{
			double ar = t[idx];
			double ai = t[idx+1];
			double br = in[idx];
			double bi = in[idx+1];

			out[idx] = ar*br - ai*bi;
			out[idx+1] = ar*bi + ai*br;
	}
#else
	CART2D_DATA_SPECTRAL_FOR_IDX(
			retval.spectral_space_data[idx] = spectral_space_data[idx]*i_array_data.spectral_space_data[idx];
	);
#endif

	retval.spectral_zeroAliasingModes();

	return retval;
}


DataSpectral DataSpectral::restrict(
		const DataSpectral &i_array_data
)
{

	DataSpectral out = *this;
	out = i_array_data.spectral_returnWithDifferentModes(out.cart2DDataConfig);


	return out;
}


DataSpectral DataSpectral::pad_zeros(
		const DataSpectral &i_array_data
)
{

	DataSpectral out = *this;
	out = i_array_data.spectral_returnWithDifferentModes(out.cart2DDataConfig);

	return out;
}


}}}
