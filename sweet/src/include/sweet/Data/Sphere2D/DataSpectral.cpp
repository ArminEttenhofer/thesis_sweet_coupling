

#include <cmath>
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


#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>

#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>

#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Memory/parmemcpy.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Fatal.hpp>



namespace sweet {
namespace Data {
namespace Sphere2D {


DataSpectral::DataSpectral(
		const Sphere2D::Config *i_sphere2DDataConfig
)	:
	sphere2DDataConfig(i_sphere2DDataConfig),
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_sphere2DDataConfig != 0);

	setup(i_sphere2DDataConfig);
}


DataSpectral::DataSpectral(
		const Sphere2D::Config &i_sphere2DDataConfig
)	:
	Sphere2D::DataSpectral(&i_sphere2DDataConfig)
{
}


DataSpectral::DataSpectral(
		const Sphere2D::Config *i_sphere2DDataConfig,
		const double &i_value
)	:
	sphere2DDataConfig(i_sphere2DDataConfig),
	spectral_space_data(nullptr)
{
	SWEET_ASSERT(i_sphere2DDataConfig != 0);

	setup(i_sphere2DDataConfig);
	spectral_setValue(i_value);
}



DataSpectral::DataSpectral()	:
	sphere2DDataConfig(nullptr),
	spectral_space_data(nullptr)
{
}




DataSpectral::DataSpectral(
		const Sphere2D::DataSpectral &i_sph_data
)	:
	sphere2DDataConfig(i_sph_data.sphere2DDataConfig),
	spectral_space_data(nullptr)
{
	if (i_sph_data.sphere2DDataConfig == nullptr)
		return;

	alloc_data();

	operator=(i_sph_data);
}




DataSpectral::DataSpectral(
		Sphere2D::DataSpectral &&i_sph_data
)	:
	sphere2DDataConfig(i_sph_data.sphere2DDataConfig),
	spectral_space_data(nullptr)
{
	if (i_sph_data.sphere2DDataConfig == nullptr)
		return;

	std::swap(spectral_space_data, i_sph_data.spectral_space_data);
}




/**
 * Run validation checks to make sure that the physical and spectral spaces match in size
 */

 DataSpectral::DataSpectral(
		const Sphere2D::DataGrid &i_sphere2d_data_physical
)
{
	setup(i_sphere2d_data_physical.sphere2DDataConfig);

	loadSphere2DDataGrid(i_sphere2d_data_physical);
}



DataSpectral::~DataSpectral()
{
	clear();
}


void DataSpectral::clear()
{
	if (spectral_space_data != nullptr)
	{
		sweet::Memory::MemBlockAlloc::free(spectral_space_data, sphere2DDataConfig->spectral_array_data_number_of_elements * sizeof(TComplex));
		spectral_space_data = nullptr;

		sphere2DDataConfig = nullptr;
	}
}




void DataSpectral::swapWithConfig(
		Sphere2D::DataSpectral &i_sphere2DData
)
{
	std::swap(sphere2DDataConfig, i_sphere2DData.sphere2DDataConfig);

	std::swap(spectral_space_data, i_sphere2DData.spectral_space_data);
}



void DataSpectral::swap(
		Sphere2D::DataSpectral &i_sphere2DData
)
{
	SWEET_ASSERT(sphere2DDataConfig == i_sphere2DData.sphere2DDataConfig);

	std::swap(spectral_space_data, i_sphere2DData.spectral_space_data);
}

void DataSpectral::check(
		const Sphere2D::Config *i_sphere2DDataConfig
)	const
{
	SWEET_ASSERT(sphere2DDataConfig->grid_num_lat == i_sphere2DDataConfig->grid_num_lat);
	SWEET_ASSERT(sphere2DDataConfig->grid_num_lon == i_sphere2DDataConfig->grid_num_lon);

	SWEET_ASSERT(sphere2DDataConfig->spectral_modes_m_max == i_sphere2DDataConfig->spectral_modes_m_max);
	SWEET_ASSERT(sphere2DDataConfig->spectral_modes_n_max == i_sphere2DDataConfig->spectral_modes_n_max);
}



DataSpectral& DataSpectral::operator=(
		const Sphere2D::DataSpectral &i_sph_data
)
{
	assert (i_sph_data.sphere2DDataConfig != nullptr);

	if (sphere2DDataConfig == nullptr)
		setup(i_sph_data.sphere2DDataConfig);

	sweet::Memory::parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(TComplex)*sphere2DDataConfig->spectral_array_data_number_of_elements);

	return *this;
}




/**
 * This function implements copying the spectral data only.
 *
 * This becomes handy if coping with data which should be only transformed without dealiasing.
 */

DataSpectral& DataSpectral::load_nodealiasing(
		const Sphere2D::DataSpectral &i_sph_data
)
{
	if (sphere2DDataConfig == nullptr)
		SWEETErrorFatal("sphere2DDataConfig not initialized");

	sweet::Memory::parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(TComplex)*sphere2DDataConfig->spectral_array_data_number_of_elements);

	return *this;
}



DataSpectral& DataSpectral::operator=(
		Sphere2D::DataSpectral &&i_sph_data
)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_sph_data.sphere2DDataConfig);

	std::swap(spectral_space_data, i_sph_data.spectral_space_data);

	return *this;
}



DataSpectral DataSpectral::spectral_returnWithDifferentModes(
		const Sphere2D::Config *i_sphere2DDataConfig
)	const
{
	Sphere2D::DataSpectral out(i_sphere2DDataConfig);

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
		for (int m = 0; m <= out.sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			TComplex *dst = &out.spectral_space_data[out.sphere2DDataConfig->getArrayIndexByModes(m, m)];
			TComplex *src = &spectral_space_data[sphere2DDataConfig->getArrayIndexByModes(m, m)];

			std::size_t size = sizeof(TComplex)*(out.sphere2DDataConfig->spectral_modes_n_max-m+1);
			sweet::Memory::parmemcpy(dst, src, size);
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
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			TComplex *dst = &out.spectral_space_data[out.sphere2DDataConfig->getArrayIndexByModes(m, m)];
			TComplex *src = &spectral_space_data[sphere2DDataConfig->getArrayIndexByModes(m, m)];

			std::size_t size = sizeof(TComplex)*(sphere2DDataConfig->spectral_modes_n_max-m+1);
			sweet::Memory::parmemcpy(dst, src, size);
		}
	}

	return out;
}



/*
 * Setup spectral sphere2D data based on data in physical space
 */
void DataSpectral::loadSphere2DDataGrid(
		const Sphere2D::DataGrid &i_sphere2DDataGrid
)
{
	/**
	 * Warning: The sphat_2_SH function is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 * Hence, we create a copy
	 */
	Sphere2D::DataGrid tmp(i_sphere2DDataGrid);
	spat_to_SH(sphere2DDataConfig->shtns, tmp.grid_space_data, spectral_space_data);
}



/*
 * Return the data converted to physical space
 */
Sphere2D::DataGrid DataSpectral::getSphere2DDataGrid()	const
{
	/**
	 * Warning: This is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 */
	Sphere2D::DataSpectral tmp(*this);
	Sphere2D::DataGrid retval(sphere2DDataConfig);
	SH_to_spat(sphere2DDataConfig->shtns, tmp.spectral_space_data, retval.grid_space_data);

	return retval;
}


/*
 * Return the data converted to physical space
 *
 * Alias for "getSphere2DDataGrid"
 */
Sphere2D::DataGrid DataSpectral::toGrid()	const
{
	/**
	 * Warning: This is an in-situ operation.
	 * Therefore, the data in the source array will be destroyed.
	 */
	Sphere2D::DataSpectral tmp(*this);
	Sphere2D::DataGrid retval(sphere2DDataConfig);

	SH_to_spat(sphere2DDataConfig->shtns, tmp.spectral_space_data, retval.grid_space_data);

	return retval;
}


Sphere2DComplex::DataGrid DataSpectral::getSphere2DDataGridComplex()	const
{
	Sphere2DComplex::DataGrid out(sphere2DDataConfig);

	/*
	 * WARNING:
	 * We have to use a temporary array here because of destructive SH transformations
	 */
	Sphere2D::DataSpectral tmp_spectral(*this);
	Sphere2D::DataGrid tmp_physical(sphere2DDataConfig);
	SH_to_spat(sphere2DDataConfig->shtns, tmp_spectral.spectral_space_data, tmp_physical.grid_space_data);

	sweet::Memory::parmemcpy(out.grid_space_data, tmp_physical.grid_space_data, sizeof(double)*sphere2DDataConfig->grid_number_elements);

	return out;
}



std::complex<double>& DataSpectral::operator[](std::size_t i)
{
	return spectral_space_data[i];
}

const std::complex<double>& DataSpectral::operator[](std::size_t i)	const
{
	return spectral_space_data[i];
}


DataSpectral DataSpectral::operator+(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	check(i_sph_data.sphere2DDataConfig);

	Sphere2D::DataSpectral out(sphere2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

	return out;
}



DataSpectral& DataSpectral::operator+=(
		const Sphere2D::DataSpectral &i_sph_data
)
{
	check(i_sph_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

	return *this;
}


DataSpectral DataSpectral::operator+(
		double i_value
)	const
{
	Sphere2D::DataSpectral out(*this);
	out.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
	return out;
}


const Sphere2D::DataSpectral& DataSpectral::operator+=(
		double i_value
)	const
{
	spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
	return *this;
}




DataSpectral DataSpectral::operator-(
		double i_value
)	const
{
	Sphere2D::DataSpectral out(*this);
	out.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);
	return out;
}



DataSpectral DataSpectral::operator-()	const
{
	Sphere2D::DataSpectral out(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out.spectral_space_data[idx] = -spectral_space_data[idx];

	return out;
}



DataSpectral DataSpectral::operator-(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	check(i_sph_data.sphere2DDataConfig);

	Sphere2D::DataSpectral out(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

	return out;
}



DataSpectral& DataSpectral::operator-=(
		double i_value
)
{
	spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);
	return *this;
}



DataSpectral& DataSpectral::operator-=(
		const Sphere2D::DataSpectral &i_sph_data
)
{
	check(i_sph_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

	return *this;
}


DataSpectral DataSpectral::operator*(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	check(i_sph_data.sphere2DDataConfig);

	Sphere2D::DataGrid a = getSphere2DDataGrid();
	Sphere2D::DataGrid b = i_sph_data.toGrid();

	Sphere2D::DataGrid mul(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		mul.grid_space_data[i] = a.grid_space_data[i]*b.grid_space_data[i];

	Sphere2D::DataSpectral out(mul);

	return out;
}



DataSpectral DataSpectral::operator/(
		double i_value
)	const
{
	Sphere2D::DataSpectral out(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

	return out;
}


const Sphere2D::DataSpectral& DataSpectral::operator/=(
		double i_value
)	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] /= i_value;

	return *this;
}





DataSpectral DataSpectral::operator/(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	check(i_sph_data.sphere2DDataConfig);

	Sphere2D::DataGrid a = getSphere2DDataGrid();
	Sphere2D::DataGrid b = i_sph_data.toGrid();

	Sphere2D::DataGrid div(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		div.grid_space_data[i] = a.grid_space_data[i]/b.grid_space_data[i];

	Sphere2D::DataSpectral out(div);

	return out;
}



DataSpectral DataSpectral::operator*(
		double i_value
)	const
{
	Sphere2D::DataSpectral out(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

	return out;
}




const Sphere2D::DataSpectral& DataSpectral::operator*=(
		double i_value
)	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] *= i_value;

	return *this;
}




const Sphere2D::DataSpectral& DataSpectral::operator*=(
		const std::complex<double> &i_value
)	const
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		spectral_space_data[idx] *= i_value;

	return *this;
}


DataSpectral DataSpectral::spectralElementwiseMultiplication(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	check(i_sph_data.sphere2DDataConfig);

	Sphere2D::DataSpectral out(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
		out.spectral_space_data[i] = this->spectral_space_data[i] * i_sph_data.spectral_space_data[i];

	return out;
}


DataSpectral DataSpectral::operator_scalar_sub_this(
		double i_value
)	const
{
	Sphere2D::DataSpectral out(sphere2DDataConfig);
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int idx = 0; idx < sphere2DDataConfig->spectral_array_data_number_of_elements; idx++)
		out.spectral_space_data[idx] = -spectral_space_data[idx];

	out.spectral_space_data[0] = i_value*std::sqrt(4.0*M_PI) + out.spectral_space_data[0];
	return out;
}



bool DataSpectral::setup(
	const Sphere2D::Config *i_sphere2DDataConfig
)
{
	sphere2DDataConfig = i_sphere2DDataConfig;
	return alloc_data();
}



bool DataSpectral::isSetup()
{
	return sphere2DDataConfig != nullptr;
}


bool DataSpectral::setup(
	const Sphere2D::Config &i_sphere2DDataConfig
)
{
	return setup(&i_sphere2DDataConfig);
}



bool DataSpectral::setup(
	const Sphere2D::Config *i_sphere2DDataConfig,
	double i_value
)
{
	sphere2DDataConfig = i_sphere2DDataConfig;
	bool retval = alloc_data();
	spectral_setValue(i_value);

	return retval;
}



bool DataSpectral::alloc_data()
{
	SWEET_ASSERT(spectral_space_data == nullptr);
	spectral_space_data = sweet::Memory::MemBlockAlloc::alloc<TComplex>(sphere2DDataConfig->spectral_array_data_number_of_elements * sizeof(TComplex));
	return true;
}



void DataSpectral::setup_if_required(
	const Sphere2D::Config *i_sphere2DDataConfig
)
{
	if (sphere2DDataConfig != nullptr)
		return;

	setup(i_sphere2DDataConfig);
}




DataSpectral DataSpectral::spectral_solve_helmholtz(
		const double &i_a,
		const double &i_b,
		double r
)
{
	Sphere2D::DataSpectral out(*this);

	const double a = i_a;
	const double b = i_b/(r*r);

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


DataSpectral DataSpectral::spectral_solve_helmholtz_higher_order(
		const double & a,
		const std::array<double, 4> & b,
		double r
)
{
	Sphere2D::DataSpectral out(*this);

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			double laplace_op_2 = - (double)n*((double)n+1.0)/(r*r);
			double laplace_op_4 = laplace_op_2 * laplace_op_2;
			double laplace_op_6 = laplace_op_4 * laplace_op_2;
			double laplace_op_8 = laplace_op_4 * laplace_op_4;
			io_data /= (a + b[0] * laplace_op_2 + b[1] * laplace_op_4 + b[2] * laplace_op_6 + b[3] * laplace_op_8);
		}
	);

	return out;
}


DataSpectral DataSpectral::spectral_solve_laplace(
		double r
)
{
	Sphere2D::DataSpectral out(*this);

	const double b = 1.0/(r*r);

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			if (n == 0)
				io_data = 0;
			else
				io_data /= (-b*(double)n*((double)n+1.0));
		}
	);

	return out;
}


const Sphere2D::DataSpectral& DataSpectral::spectral_truncate()	const
{
	Sphere2D::DataGrid tmp(sphere2DDataConfig);

	SH_to_spat(sphere2DDataConfig->shtns, spectral_space_data, tmp.grid_space_data);
	spat_to_SH(sphere2DDataConfig->shtns, tmp.grid_space_data, spectral_space_data);

	return *this;
}


void DataSpectral::spectral_update_lambda(
		std::function<void(int,int,TComplex&)> i_lambda
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			i_lambda(n, m, spectral_space_data[idx]);
			idx++;
		}
	}
}

const std::complex<double>& DataSpectral::spectral_get_DEPRECATED(
		int i_n,
		int i_m
)	const
{
	static const std::complex<double> zero = {0,0};

	if (i_n < 0 ||  i_m < 0)
		return zero;

	if (i_n > sphere2DDataConfig->spectral_modes_n_max)
		return zero;

	if (i_m > sphere2DDataConfig->spectral_modes_m_max)
		return zero;

	if (i_m > i_n)
		return zero;

	assert (i_m <= sphere2DDataConfig->spectral_modes_m_max);
	return spectral_space_data[sphere2DDataConfig->getArrayIndexByModes(i_n, i_m)];
}


const std::complex<double>& DataSpectral::spectral_get_(
		int i_n,
		int i_m
)	const
{
	SWEET_ASSERT(i_n >= 0 && i_m >= 0);
	SWEET_ASSERT(i_n <= sphere2DDataConfig->spectral_modes_n_max);
	SWEET_ASSERT(i_m <= sphere2DDataConfig->spectral_modes_m_max);
	SWEET_ASSERT(i_m <= i_n);

	return spectral_space_data[sphere2DDataConfig->getArrayIndexByModes(i_n, i_m)];
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

	if (i_n > sphere2DDataConfig->spectral_modes_n_max)
		SWEETErrorFatal("Out of boundary b");

	if (i_m > sphere2DDataConfig->spectral_modes_m_max)
		SWEETErrorFatal("Out of boundary c");

	if (i_m > i_n)
		SWEETErrorFatal("Out of boundary d");

	assert (i_m <= sphere2DDataConfig->spectral_modes_m_max);
#endif

	spectral_space_data[sphere2DDataConfig->getArrayIndexByModes(i_n, i_m)] = i_data;
}


void DataSpectral::spectral_setZero()
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int i = 0; i < sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
		spectral_space_data[i] = {0,0};
}


void DataSpectral::spectral_setValue(
		const std::complex<double> &i_value
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int i = 0; i < sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
		spectral_space_data[i] = i_value;
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
	for (int i = 0; i < sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
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

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (int i = 0; i < sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
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
	//std::complex<double> c = 0;

	//m=0 case - weight 1
	std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(0, 0);
	for (int n = 0; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
		idx++;
	}

	//m>0 case - weight 2, as they appear twice
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			sum += 2.0*spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
			idx++;
		}
	}

	return sum.real();
}


std::complex<double> DataSpectral::spectral_reduce_min()	const
{
	std::complex<double> error = std::numeric_limits<double>::infinity();

	for (int j = 0; j < sphere2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		error.real(std::min(spectral_space_data[j].real(), error.real()));
		error.imag(std::min(spectral_space_data[j].imag(), error.imag()));
	}

	return error;
}


std::complex<double> DataSpectral::spectral_reduce_max()	const
{
	std::complex<double> error = -std::numeric_limits<double>::infinity();

	for (int j = 0; j < sphere2DDataConfig->spectral_array_data_number_of_elements; j++)
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
	for (int j = 0; j < sphere2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
		error = std::max(std::abs(w), error);
	}

	return error;
}


double DataSpectral::spectral_reduce_max_abs(std::size_t rnorm)	const
{

	assert ((int)rnorm <= sphere2DDataConfig->spectral_modes_m_max);
	assert ((int)rnorm <= sphere2DDataConfig->spectral_modes_n_max);

	double error = -std::numeric_limits<double>::infinity();
	std::complex<double> w = {0,0};

	for (std::size_t m = 0; m <= rnorm; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (std::size_t n = m; n <= rnorm; n++)
		{
			w = spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
			error = std::max(std::abs(w), error);
			idx++;
		}
	}

	return error;
}


double DataSpectral::spectral_reduce_min_abs()	const
{
	double error = std::numeric_limits<double>::infinity();
	std::complex<double> w = {0,0};
	for (int j = 0; j < sphere2DDataConfig->spectral_array_data_number_of_elements; j++)
	{
		w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
		error = std::min(std::abs(w), error);
	}

	return error;
}

bool DataSpectral::spectral_reduce_is_any_nan_or_inf()	const
{
	bool retval = false;

#if SWEET_THREADING_SPACE
	#pragma omp parallel for simd reduction(|:retval)
#endif
	for (int j = 0; j < sphere2DDataConfig->spectral_array_data_number_of_elements; j++)
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
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			if (std::abs(spectral_space_data[idx]) < i_abs_threshold)
				std::cout << 0 << "\t";
			else
				std::cout << spectral_space_data[idx] << "\t";
			idx++;
		}
		std::cout << std::endl;
	}
}

void DataSpectral::spectral_structure_print(
		int i_precision,
		double i_abs_threshold
)	const
{
	std::cout << std::setprecision(i_precision);

	std::cout << "m \\ n ----->" << std::endl;
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		//std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			std::cout << "(" << m << "," << n << ")" << "\t";
		}
		std::cout << std::endl;
	}
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

	file << "#SWEET_SPHERE2D_SPECTRAL_CPLX_DATA_ASCII" << std::endl;


	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			const std::complex<double> &value = spectral_space_data[idx];
			file << "(" << value.real() << ", " << value.imag() << ")";
			idx++;

			if (n < sphere2DDataConfig->spectral_modes_n_max)
				file << i_separator;
			else
				file << std::endl;

		}
	}

	return true;
}


void DataSpectral::spectrum_file_write(
		const std::string &i_filename,
		const char *i_title,
		int i_precision
)	const
{
	std::ofstream file(i_filename, std::ios_base::trunc);

	file << std::setprecision(i_precision);
	file << "#SWEET_SPHERE2D_SPECTRAL_DATA_ASCII" << std::endl;

	file << "#TI " << i_title << std::endl;

	// Use 0 to make it processable by python
	file << "0\t";

	std::complex<double> w = {0,0};
	std::vector<double> sum(sphere2DDataConfig->spectral_modes_n_max+1,0);
	std::vector<double> sum_squared(sphere2DDataConfig->spectral_modes_n_max+1,0);
	std::vector<double> max(sphere2DDataConfig->spectral_modes_n_max+1,0);

	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			w = spectral_space_data[idx];
			idx++;

			sum[n]         += std::abs(w);
			sum_squared[n] += std::abs(w * w);
			if (std::abs(w) >= max[n]) max[n] = std::abs(w);
		}
	}

	file << sphere2DDataConfig->spectral_modes_m_max << " "
			<< sphere2DDataConfig->spectral_modes_n_max << std::endl;
	for (int n = 0; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
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
		file << "#SWEET_SPHERE2D_SPECTRAL_ABS_EVOL_ASCII" << std::endl;
		file << "#TI " << i_title << std::endl;
		file << "0\t" << std::endl; // Use 0 to make it processable by python
		file << "(n_max=" <<sphere2DDataConfig->spectral_modes_n_max << " m_max="
				<< sphere2DDataConfig->spectral_modes_n_max << ")" << std::endl;
		file << "timestamp\t" ;
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
		{
			//std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
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
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
		{
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
		file << "#SWEET_SPHERE2D_SPECTRAL_PHASE_EVOL_ASCII" << std::endl;
		file << "#TI " << i_title << std::endl;
		file << "0\t" << std::endl; // Use 0 to make it processable by python
		file << "(n_max=" <<sphere2DDataConfig->spectral_modes_n_max << " m_max="
				<< sphere2DDataConfig->spectral_modes_n_max << ")" << std::endl;
		file << "timestamp\t" ;
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
		{
			//std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
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

	//std::cout << "n" << " " << "m" << " " << "norm" <<std::endl;
	file << i_time << "\t";
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
		{
			w = spectral_space_data[idx];
			wphase = std::arg(w); // std::abs(w * std::conj(w));

			//file << "(" << n << "," << m << ")\t" <<std::endl;
			file <<  wphase << "\t"; //<<std::endl;;
			//std::cout << n << " " << m << " " << wabs <<std::endl;

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
	file << "DATA_TYPE SH_DATA" << std::endl;
	file << "MODES_M_MAX " << sphere2DDataConfig->spectral_modes_m_max << std::endl;
	file << "MODES_N_MAX " << sphere2DDataConfig->spectral_modes_n_max << std::endl;
	file << "GRID_TYPE GAUSSIAN" << std::endl;
	file << "NUM_ELEMENTS " << sphere2DDataConfig->spectral_array_data_number_of_elements << std::endl;
	file << "FIN" << std::endl;

	file.write((const char*)spectral_space_data, sizeof(std::complex<double>)*sphere2DDataConfig->spectral_array_data_number_of_elements);

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
	sweet::Data::Sphere2D::DataGrid tmp(sphere2DDataConfig);
	tmp.file_read_csv_grid(i_filename.c_str());
	loadSphere2DDataGrid(tmp);
}


void DataSpectral::file_read_binary_spectral(
		const std::string &i_filename
)
{
	std::ifstream file(i_filename, std::ios_base::binary);

	if (!file.is_open())
		SWEETErrorFatal("Error while opening file " + i_filename);

	std::string magic;
	std::getline(file, magic);

	if (magic != "SWEET")
		SWEETErrorFatal("Magic code 'SWEET' not found");

	std::string data_type;
	int num_lon = -1;
	int num_lat = -1;

	while (true)
	{
		std::string grid_type = "";
		std::string buf;
		file >> buf;

		if (buf == "FIN")
			break;

		if (buf == "DATA_TYPE")
		{
			// load data type
			file >> data_type;
			continue;
		}

		if (buf == "MODES_M_MAX")
		{
			file >> buf;
			num_lon = std::stoi(buf);
			continue;
		}

		if (buf == "MODES_N_MAX")
		{
			file >> buf;
			num_lat = std::stoi(buf);
			continue;
		}

		if (buf == "GRID_TYPE")
		{
			file >> buf;
			grid_type = buf;
			continue;
		}

		if (buf == "NUM_ELEMENTS")
		{
			file >> buf;
			continue;
		}

		SWEETErrorFatal("Unknown Tag '"+buf+"'");
	}

	// read last newline
	char nl;
	file.read(&nl, 1);

	if (data_type != "SH_DATA")
		SWEETErrorFatal("Unknown data type '"+data_type+"'");

	if (num_lon != sphere2DDataConfig->spectral_modes_m_max)
		SWEETErrorFatal("NUM_LON "+std::to_string(num_lon)+" doesn't match Sphere2DDataConfig");

	if (num_lat != sphere2DDataConfig->spectral_modes_n_max)
		SWEETErrorFatal("NUM_LAT "+std::to_string(num_lat)+" doesn't match Sphere2DDataConfig");

	file.read((char*)spectral_space_data, sizeof(std::complex<double>)*sphere2DDataConfig->spectral_array_data_number_of_elements);

	file.close();
}


void DataSpectral::normalize(
		const std::string &i_normalization
)
{
	if (i_normalization == "avg_zero")
	{
		// move average value to 0
		double phi_min = getSphere2DDataGrid().grid_reduce_min();
		double phi_max = getSphere2DDataGrid().grid_reduce_max();

		double avg = 0.5*(phi_max+phi_min);

		operator-=(avg);
	}
	else if (i_normalization == "min_zero")
	{
		// move minimum value to zero
		double phi_min = getSphere2DDataGrid().grid_reduce_min();
		operator-=(phi_min);
	}
	else if (i_normalization == "max_zero")
	{
		// move maximum value to zero
		double phi_max = getSphere2DDataGrid().grid_reduce_max();
		operator-=(phi_max);
	}
	else if (i_normalization == "")
	{
	}
	else
	{
		SWEETErrorFatal("Normalization not supported!");
	}
}


DataSpectral DataSpectral::restrict(
		const Sphere2D::DataSpectral &i_array_data
)
{
	Sphere2D::DataSpectral out = *this;
	out = i_array_data.spectral_returnWithDifferentModes(out.sphere2DDataConfig);
	return out;

}


DataSpectral DataSpectral::pad_zeros(
		const Sphere2D::DataSpectral &i_array_data
)
{

	Sphere2D::DataSpectral out = *this;
	out = i_array_data.spectral_returnWithDifferentModes(out.sphere2DDataConfig);
	return out;

}



}}}



sweet::Data::Sphere2D::DataSpectral operator*(
		double i_value,
		const sweet::Data::Sphere2D::DataSpectral &i_array_data
)
{
	return ((sweet::Data::Sphere2D::DataSpectral&)i_array_data)*i_value;
}


sweet::Data::Sphere2D::DataSpectral operator+(
		double i_value,
		const sweet::Data::Sphere2D::DataSpectral &i_array_data
)
{
	return ((sweet::Data::Sphere2D::DataSpectral&)i_array_data)+i_value;
}


sweet::Data::Sphere2D::DataSpectral operator-(
		double i_value,
		const sweet::Data::Sphere2D::DataSpectral &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}
