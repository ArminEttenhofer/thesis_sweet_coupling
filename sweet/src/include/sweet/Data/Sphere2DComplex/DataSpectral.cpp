
#include <complex>
#include <functional>
#include <array>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <functional>

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>

#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Memory/parmemcpy.hpp>


namespace sweet {
namespace Data {
namespace Sphere2DComplex {



DataSpectral::DataSpectral()	:
	sphere2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
}



DataSpectral::DataSpectral(
		const Sphere2D::Config *i_sphere2DDataConfig
)	:
	sphere2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_sphere2DDataConfig != 0);

	setup(i_sphere2DDataConfig);
}



DataSpectral::DataSpectral(
		const Sphere2DComplex::DataSpectral &i_sph_data
)	:
	sphere2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	setup(i_sph_data.sphere2DDataConfig);

	operator=(i_sph_data);
}





DataSpectral::DataSpectral(
		Sphere2DComplex::DataSpectral &&i_data
)	:
	sphere2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_data.sphere2DDataConfig);

	SWEET_ASSERT(i_data.spectral_space_data != nullptr);

	std::swap(spectral_space_data, i_data.spectral_space_data);
}



DataSpectral::DataSpectral(
		const Sphere2DComplex::DataGrid &i_sph_data
):
	sphere2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
	Sphere2DComplex::DataGrid tmp = i_sph_data;

	setup(i_sph_data.sphere2DDataConfig);

	/**
	 * Warning: This is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 */
	spat_cplx_to_SH(sphere2DDataConfig->shtns, tmp.grid_space_data, spectral_space_data);
}


 void DataSpectral::clear()
{
	if (spectral_space_data != nullptr)
	{
		sweet::Memory::MemBlockAlloc::free(spectral_space_data, sphere2DDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
		spectral_space_data = nullptr;

		sphere2DDataConfig = nullptr;
	}
}


DataSpectral::~DataSpectral()
{
	clear();
}


void DataSpectral::check_sphere2DDataConfig_identical_res(const Sphere2D::Config *i_sphere2DDataConfig)	const
{
	SWEET_ASSERT(sphere2DDataConfig->spectral_modes_m_max == i_sphere2DDataConfig->spectral_modes_m_max);
	SWEET_ASSERT(sphere2DDataConfig->spectral_modes_n_max == i_sphere2DDataConfig->spectral_modes_n_max);
}



DataSpectral DataSpectral::spectral_returnWithTruncatedModes(
		const Sphere2D::Config *i_sphere2DDataConfigTargetTruncation
)	const
{
	return spectral_returnWithDifferentModes(i_sphere2DDataConfigTargetTruncation).spectral_returnWithDifferentModes(sphere2DDataConfig);
}


DataSpectral DataSpectral::spectral_returnWithDifferentModes(
		const Sphere2D::Config *i_sphere2DDataConfigNew
)	const
{
	Sphere2DComplex::DataSpectral out(i_sphere2DDataConfigNew);

	/*
	 *  0 = invalid
	 * -1 = scale down
	 *  1 = scale up
	 */
	int scaling_mode = 0;

	if (sphere2DDataConfig->spectral_modes_m_max < out.sphere2DDataConfig->spectral_modes_m_max)
	{
		scaling_mode = 1;
	}
	else if (sphere2DDataConfig->spectral_modes_m_max > out.sphere2DDataConfig->spectral_modes_m_max)
	{
		scaling_mode = -1;
	}


	if (sphere2DDataConfig->spectral_modes_n_max < out.sphere2DDataConfig->spectral_modes_n_max)
	{
		SWEET_ASSERT(scaling_mode != -1);
		scaling_mode = 1;
	}
	else if (sphere2DDataConfig->spectral_modes_n_max > out.sphere2DDataConfig->spectral_modes_n_max)
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

	if (scaling_mode == -1)
	{
		/*
		 * more modes -> less modes
		 */


SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= out.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			int src_idx = sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
			int dst_idx = out.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);

			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[dst_idx] = spectral_space_data[src_idx];
				src_idx++;
				dst_idx++;
			}
		}
	}
	else
	{
		/*
		 * less modes -> more modes
		 */

		// zero all values
		out.spectral_setZero();


SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			int src_idx = sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
			int dst_idx = out.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);

			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[dst_idx] = spectral_space_data[src_idx];
				src_idx++;
				dst_idx++;
			}
		}
	}

	return out;
}


void DataSpectral::loadSphere2DDataGrid(
		const Sphere2DComplex::DataGrid &i_sphere2DDataGrid
)
{
	/**
	 * Warning: The sphat_2_SH function is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 * Hence, we create a copy
	 */
	Sphere2DComplex::DataGrid tmp(i_sphere2DDataGrid);
	spat_cplx_to_SH(sphere2DDataConfig->shtns, tmp.grid_space_data, spectral_space_data);
}



void DataSpectral::swap(
	Sphere2DComplex::DataSpectral &i_sphere2DData
)
{
	SWEET_ASSERT(sphere2DDataConfig == i_sphere2DData.sphere2DDataConfig);

	std::swap(spectral_space_data, i_sphere2DData.spectral_space_data);
}




Sphere2DComplex::DataGrid DataSpectral::toGrid()	const
{
	return getSphere2DDataGridComplex();
}

__attribute__((__deprecated__))
const std::complex<double>& DataSpectral::operator[](std::size_t i)	const
{
	return spectral_space_data[i];
}

__attribute__((__deprecated__))
std::complex<double>& DataSpectral::operator[](std::size_t i)
{
	return spectral_space_data[i];
}

Sphere2DComplex::DataGrid DataSpectral::getSphere2DDataGridComplex()	const
{
	Sphere2DComplex::DataGrid out(sphere2DDataConfig);

	/*
	 * WARNING:
	 * We have to use a temporary array here because of destructive SH transformations
	 */
	Sphere2DComplex::DataSpectral tmp = *this;
	SH_to_spat_cplx(sphere2DDataConfig->shtns, tmp.spectral_space_data, out.grid_space_data);

	return out;
}


DataSpectral& DataSpectral::operator=(
		const Sphere2DComplex::DataSpectral &i_sph_data
)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_sph_data.sphere2DDataConfig);

	SWEET_ASSERT(i_sph_data.spectral_space_data);
	sweet::Memory::parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(Tcomplex)*sphere2DDataConfig->spectral_complex_array_data_number_of_elements);

	return *this;
}


DataSpectral& DataSpectral::operator=(
		Sphere2DComplex::DataSpectral &&i_sph_data
)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_sph_data.sphere2DDataConfig);

	SWEET_ASSERT(i_sph_data.spectral_space_data);
	std::swap(spectral_space_data, i_sph_data.spectral_space_data);

	return *this;
}



DataSpectral DataSpectral::operator+(
		double i_value
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(*this);

	out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_sph_data;
}


DataSpectral& DataSpectral::operator+=(
		double i_value
)
{
	spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
	return *this;
}



DataSpectral DataSpectral::operator+(
		const std::complex<double> &i_value
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(*this);

	out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_sph_data;
}


Sphere2DComplex::DataSpectral DataSpectral::operator+(
		const Sphere2DComplex::DataSpectral &i_sph_data
)	const
{
	check_sphere2DDataConfig_identical_res(i_sph_data.sphere2DDataConfig);

	Sphere2DComplex::DataSpectral out_sph_data(sphere2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

	return out_sph_data;
}


DataSpectral& DataSpectral::operator+=(
		const Sphere2DComplex::DataSpectral &i_sph_data
)
{
	check_sphere2DDataConfig_identical_res(i_sph_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

	return *this;
}


DataSpectral DataSpectral::operator-(
		const std::complex<double> &i_value
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(*this);

	out_sph_data.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);

	return out_sph_data;
}


DataSpectral DataSpectral::operator-(
		const Sphere2DComplex::DataSpectral &i_sph_data
)	const
{
	check_sphere2DDataConfig_identical_res(i_sph_data.sphere2DDataConfig);

	Sphere2DComplex::DataSpectral out_sph_data(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

	return out_sph_data;
}


DataSpectral& DataSpectral::operator-=(
		const Sphere2DComplex::DataSpectral &i_sph_data
)
{
	check_sphere2DDataConfig_identical_res(i_sph_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

	return *this;
}


DataSpectral DataSpectral::operator-()	const
{
	Sphere2DComplex::DataSpectral out_sph_data(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = -spectral_space_data[idx];

	return out_sph_data;
}


DataSpectral DataSpectral::operator*(
		const double i_value
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR

	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

	return out_sph_data;
}


DataSpectral DataSpectral::operator*(
		const std::complex<double> &i_value
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(sphere2DDataConfig);


SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

	return out_sph_data;
}


const Sphere2DComplex::DataSpectral& DataSpectral::operator*=(
		const std::complex<double> &i_value
)	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] *= i_value;

	return *this;
}


DataSpectral DataSpectral::operator/(
		double i_value
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(sphere2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

	return out_sph_data;
}


const Sphere2DComplex::DataSpectral& DataSpectral::operator/=(
		const std::complex<double> &i_value
)	const
{

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		spectral_space_data[idx] /= i_value;

	return *this;
}


void DataSpectral::setup_if_required(
	const Sphere2D::Config *i_sphere2DDataConfig
)
{
	if (sphere2DDataConfig != nullptr)
		return;

	setup(i_sphere2DDataConfig);
}


void DataSpectral::setup(
		const Sphere2D::Config *i_sphere2DConfig
)
{
	// assure that the initialization is not done twice!
	SWEET_ASSERT(sphere2DDataConfig == nullptr);

	sphere2DDataConfig = i_sphere2DConfig;

	spectral_space_data = sweet::Memory::MemBlockAlloc::alloc<Tcomplex>(sphere2DDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
}

DataSpectral DataSpectral::spectral_solve_helmholtz(
		const std::complex<double> &i_a,
		const std::complex<double> &i_b,
		double r
)	const
{
	Sphere2DComplex::DataSpectral out(*this);

	const std::complex<double> a = i_a;
	const std::complex<double> b = i_b/(r*r);

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			io_data /= (a + (-b*(double)n*((double)n+1.0)));
		}
	);

	return out;
}


void DataSpectral::spectral_update_lambda(
		std::function<void(int,int,Tcomplex&)> i_lambda
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			i_lambda(n, m, spectral_space_data[idx]);
			idx++;
		}
	}
}


const std::complex<double>& DataSpectral::spectral_get_(
		int in,
		int im
)	const
{
	SWEET_ASSERT(in >= 0);
	SWEET_ASSERT(in <= sphere2DDataConfig->spectral_modes_n_max);
	SWEET_ASSERT(std::abs(im) <= sphere2DDataConfig->spectral_modes_m_max);
	SWEET_ASSERT(std::abs(im) <= in);

	return spectral_space_data[sphere2DDataConfig->getArrayIndexByModes_Complex(in, im)];
}


void DataSpectral::spectral_setZero()
{

SWEET_THREADING_SPACE_PARALLEL_FOR

	for (int n = 0; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			spectral_space_data[idx] = 0;
			idx++;
		}
	}
}


void DataSpectral::spectral_print(
		int i_precision
)	const
{
	std::cout << std::setprecision(i_precision);

	/**
	 * WARNING: This follows a different order contrast to how it is stored
	 */
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			std::cout << spectral_space_data[idx] << "\t";
		}
		std::cout << std::endl;
	}
}

}}}



sweet::Data::Sphere2DComplex::DataSpectral operator+(
		const double i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
)
{
	return ((sweet::Data::Sphere2DComplex::DataSpectral&)i_array_data)+i_value;
}


sweet::Data::Sphere2DComplex::DataSpectral operator+(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
)
{
	return i_array_data+i_value;
}


sweet::Data::Sphere2DComplex::DataSpectral operator-(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
)
{
	sweet::Data::Sphere2DComplex::DataSpectral out_sph_data(i_array_data.sphere2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < i_array_data.sphere2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = -i_array_data.spectral_space_data[idx];

	out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_sph_data;

}

sweet::Data::Sphere2DComplex::DataSpectral operator*(
		const double i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
)
{
	return i_array_data*i_value;
}


sweet::Data::Sphere2DComplex::DataSpectral operator*(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
)
{
	return i_array_data*i_value;
}


