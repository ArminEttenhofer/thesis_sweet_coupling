
#include <complex>
#include <functional>
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
#include <sweet/Data/Sphere2D/DataGrid.hpp>

#include "DataGrid.hpp"

namespace sweet {
namespace Data {
namespace Sphere2DComplex {



void DataGrid::swap(
		Sphere2DComplex::DataGrid &i_sphere2DData
)
{
	SWEET_ASSERT(sphere2DDataConfig == i_sphere2DData.sphere2DDataConfig);

	std::swap(grid_space_data, i_sphere2DData.grid_space_data);
}


DataGrid::DataGrid(
		const Sphere2D::Config *i_sphere2DDataConfig
)	:
	sphere2DDataConfig(i_sphere2DDataConfig),
	grid_space_data(nullptr)
{
	setup(i_sphere2DDataConfig);
}



DataGrid::DataGrid()	:
	sphere2DDataConfig(nullptr),
	grid_space_data(nullptr)
{
}



DataGrid::DataGrid(
		const Sphere2DComplex::DataGrid &i_data
)	:
	sphere2DDataConfig(i_data.sphere2DDataConfig),
	grid_space_data(nullptr)
{
	setup(i_data.sphere2DDataConfig);

	operator=(i_data);
}



DataGrid::DataGrid(
		const Sphere2D::DataGrid &i_data
)	:
	sphere2DDataConfig(i_data.sphere2DDataConfig),
	grid_space_data(nullptr)
{
	setup(i_data.sphere2DDataConfig);

	operator=(i_data);
}


DataGrid::~DataGrid()
{
	if (grid_space_data != nullptr)
		sweet::Memory::MemBlockAlloc::free(grid_space_data, sphere2DDataConfig->grid_number_elements * sizeof(std::complex<double>));
}





void DataGrid::setup(
		const Sphere2D::Config *i_sphere2DDataConfig
)
{
	sphere2DDataConfig = i_sphere2DDataConfig;

	grid_space_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(sphere2DDataConfig->grid_number_elements * sizeof(std::complex<double>));
}




void DataGrid::setup_if_required(
		const Sphere2D::Config *i_sphere2DDataConfig
)
{
	if (sphere2DDataConfig != nullptr)
	{
		SWEET_ASSERT(grid_space_data != nullptr);
		return;
	}

	sphere2DDataConfig = i_sphere2DDataConfig;
	grid_space_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(sphere2DDataConfig->grid_number_elements * sizeof(std::complex<double>));
}




void DataGrid::loadRealImag(
		const Sphere2D::DataGrid &i_re,
		const Sphere2D::DataGrid &i_im
)
{

#if 1

	double *t = (double*)grid_space_data;
	double *in_re = (double*)i_re.grid_space_data;
	double *in_im = (double*)i_im.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
	{
		t[2*i] = in_re[i];
		t[2*i+1] = in_im[i];
	}

#else

	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
	{
		grid_space_data[i].real(i_re.grid_space_data[i]);
		grid_space_data[i].imag(i_im.grid_space_data[i]);
	}
#endif
}


 void DataGrid::check(
		const Sphere2D::Config *i_sphere2DDataConfig
)	const
{
	SWEET_ASSERT(sphere2DDataConfig->grid_num_lat == i_sphere2DDataConfig->grid_num_lat);
	SWEET_ASSERT(sphere2DDataConfig->grid_num_lon == i_sphere2DDataConfig->grid_num_lon);
}




Sphere2DComplex::DataGrid& DataGrid::operator=(
		const Sphere2DComplex::DataGrid &i_data
)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_data.sphere2DDataConfig);

	sweet::Memory::parmemcpy(grid_space_data, i_data.grid_space_data, sizeof(std::complex<double>)*sphere2DDataConfig->grid_number_elements);

	return *this;
}



Sphere2DComplex::DataGrid& DataGrid::operator=(
		Sphere2DComplex::DataGrid &&i_data
)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_data.sphere2DDataConfig);

	std::swap(grid_space_data, i_data.grid_space_data);

	return *this;
}




Sphere2DComplex::DataGrid& DataGrid::operator=(
		const Sphere2D::DataGrid &i_data
)
{
	if (sphere2DDataConfig == nullptr)
		setup(i_data.sphere2DDataConfig);

#if 1

	/*
	 * We implemented this double precision version
	 * since SIMDization with std::complex doesn't seem to work on LLVM
	 */
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
	{
		t[2*i] = in[i];
		t[2*i+1] = 0;
	}

#else

	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] = i_data.grid_space_data[i];

#endif

	return *this;
}


Sphere2DComplex::DataGrid DataGrid::operator+(
		const Sphere2DComplex::DataGrid &i_data
)	const
{
	check(i_data.sphere2DDataConfig);

	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	/*
	 * We implemented this double precision version
	 * since SIMDization with std::complex doesn't seem to work on LLVM
	 */
	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		out[i] = t[i] + in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i] + i_data.grid_space_data[i];

#endif

	return retval;
}



Sphere2DComplex::DataGrid& DataGrid::operator+=(
		const Sphere2DComplex::DataGrid &i_data
)
{
	check(i_data.sphere2DDataConfig);

#if 1

	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		t[i] += in[i];

#else


	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] += i_data.grid_space_data[i];

#endif

	return *this;
}


Sphere2DComplex::DataGrid DataGrid::operator+(
		const std::complex<double> &i_value
)	const
{
	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i] + i_value.real();
		out[i+1] = t[i+1] + i_value.imag();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]+i_value;

	#endif

	return retval;
}



Sphere2DComplex::DataGrid& DataGrid::operator-=(
		const Sphere2DComplex::DataGrid &i_data
)
{
	check(i_data.sphere2DDataConfig);

#if 1

	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		t[i] -= in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] -= i_data.grid_space_data[i];

#endif

	return *this;
}



Sphere2DComplex::DataGrid DataGrid::operator-(
		const Sphere2DComplex::DataGrid &i_data
)	const
{
	check(i_data.sphere2DDataConfig);

	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		out[i] = t[i] - in[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i] - i_data.grid_space_data[i];
#endif

	return retval;
}



Sphere2DComplex::DataGrid DataGrid::operator-()
{
	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		out[i] = -t[i];

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = -grid_space_data[i];

#endif
	return retval;
}



Sphere2DComplex::DataGrid DataGrid::operator-(
		const std::complex<double> &i_value
)	const
{
	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i] - i_value.real();
		out[i+1] = t[i+1] - i_value.imag();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]-i_value;

#endif

	return retval;
}





Sphere2DComplex::DataGrid DataGrid::operator*(
		const Sphere2DComplex::DataGrid &i_data
)	const
{
	check(i_data.sphere2DDataConfig);

	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;
	double *in = (double*)i_data.grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i]*in[i] - t[i+1]*in[i+1];
		out[i+1] = t[i]*in[i+1] + t[i+1]*in[i];
	}


#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t in = 0; in < sphere2DDataConfig->grid_number_elements; in++)
		retval.grid_space_data[in] = grid_space_data[in]*i_data.grid_space_data[in];

#endif

	return retval;
}



Sphere2DComplex::DataGrid DataGrid::operator*(
		const double i_value
)	const
{
	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		out[i] = t[i]*i_value;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]*i_value;

#endif

	return retval;
}



Sphere2DComplex::DataGrid DataGrid::operator*(
		const std::complex<double> &i_value
)	const
{
	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *out = (double*)retval.grid_space_data;
	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i+=2)
	{
		out[i] = t[i]*i_value.real() - t[i+1]*i_value.imag();
		out[i+1] = t[i]*i_value.imag() + t[i+1]*i_value.real();
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]*i_value;

#endif

	return retval;
}




const Sphere2DComplex::DataGrid& DataGrid::operator*=(
		const double i_value
)	const
{

#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		t[i] *= i_value;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] *= i_value;

#endif

	return *this;
}



Sphere2DComplex::DataGrid DataGrid::operator/(
		const Sphere2DComplex::DataGrid &i_data
)	const
{
	check(i_data.sphere2DDataConfig);

	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

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
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i+=2)
	{
		double &a = in[i];
		double &b = in[i+1];
		double denom = a*a + b*b;

		out[i] = (t[i]*a + t[i+1]*b)/denom;
		out[i+1] = (-t[i]*b + t[i+1]*a)/denom;
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t in = 0; in < sphere2DDataConfig->grid_number_elements; in++)
		retval.grid_space_data[in] = grid_space_data[in]/i_data.grid_space_data[in];

#endif

	return retval;
}



Sphere2DComplex::DataGrid DataGrid::operator/(
		double i_value
)	const
{
	Sphere2DComplex::DataGrid retval(sphere2DDataConfig);

#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		t[i] /= i_value;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		retval.grid_space_data[i] = grid_space_data[i]/i_value;

#endif

	return retval;
}



/*
 * Set values for all latitude and longitude degrees
 *
 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
 */
void DataGrid::grid_update_lambda(
		std::function<void(double,double,std::complex<double>&)> i_lambda
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR

#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS

	for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
	{
		double lon_degree = ((double)i/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

		for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
		{
			//double colatitude = acos(shtns->ct[j]);

			/*
			 * Colatitude is 0 at the north pole and 180 at the south pole
			 *
			 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
			 */
			//double lat_degree = M_PI*0.5 - colatitude;
			double lat_degree = sphere2DDataConfig->lat[j];

			i_lambda(lon_degree, lat_degree, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
		}
	}

#else

	for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
	{
		double lat_degree = sphere2DDataConfig->lat[jlat];

		for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++)
		{
			double lon_degree = ((double)ilon/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

			//double colatitude = acos(shtns->ct[j]);

			/*
			 * Colatitude is 0 at the north pole and 180 at the south pole
			 *
			 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
			 */
			//double lat_degree = M_PI*0.5 - colatitude;

			i_lambda(lon_degree, lat_degree, grid_space_data[jlat*sphere2DDataConfig->grid_num_lon + ilon]);
		}
	}

#endif
}


void DataGrid::grid_update_lambda_array(
		std::function<void(int,int,std::complex<double>&)> i_lambda
)
{
#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
	{
		for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
		{
			i_lambda(i, j, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
		}
	}

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
	{
		for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++)
		{
			i_lambda(ilon, jlat, grid_space_data[jlat*sphere2DDataConfig->grid_num_lon + ilon]);
		}
	}

#endif
}


void DataGrid::grid_update_lambda_gaussian_grid(
		std::function<void(double,double,std::complex<double>&)> i_lambda
)
{

#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
	{
		double lon_degree = ((double)i/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

		for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
		{
			double sin_phi = sphere2DDataConfig->lat_gaussian[j];

			i_lambda(lon_degree, sin_phi, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
		}
	}
#else

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
	{
		double sin_phi = sphere2DDataConfig->lat_gaussian[jlat];

		for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++)
		{
			double lon_degree = ((double)ilon/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

			i_lambda(lon_degree, sin_phi, grid_space_data[jlat*sphere2DDataConfig->grid_num_lon + ilon]);
		}
	}
#endif
}



void DataGrid::grid_update_lambda_cogaussian_grid(
		std::function<void(double,double,std::complex<double>&)> i_lambda
)
{


#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int i = 0; i < sphere2DDataConfig->grid_num_lon; i++)
	{
		double lon_degree = (((double)i)/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

		for (int j = 0; j < sphere2DDataConfig->grid_num_lat; j++)
		{
			double cos_phi = sphere2DDataConfig->lat_cogaussian[j];

			/*
			 * IDENTITAL FORMULATION
			double mu = shtns->ct[j];
			double comu = sqrt(1.0-mu*mu);
			*/

			i_lambda(lon_degree, cos_phi, grid_space_data[i*sphere2DDataConfig->grid_num_lat + j]);
		}
	}
#else

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
	{
		double cos_phi = sphere2DDataConfig->lat_cogaussian[jlat];

		for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++)
		{
			double lon_degree = (((double)ilon)/(double)sphere2DDataConfig->grid_num_lon)*2.0*M_PI;

			/*
			 * IDENTITAL FORMULATION
			double mu = shtns->ct[j];
			double comu = sqrt(1.0-mu*mu);
			*/

			i_lambda(lon_degree, cos_phi, grid_space_data[jlat*sphere2DDataConfig->grid_num_lon + ilon]);
		}
	}
#endif
}


void DataGrid::grid_update_lambda_sinphi_grid(
		std::function<void(double,double,std::complex<double>&)> i_lambda
)
{
	grid_update_lambda_gaussian_grid(i_lambda);
}

void DataGrid::grid_update_lambda_cosphi_grid(
		std::function<void(double,double,std::complex<double>&)> i_lambda
)
{
	grid_update_lambda_cogaussian_grid(i_lambda);
}



void DataGrid::grid_setZero()
{
#if 1

	double *t = (double*)grid_space_data;

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < 2*sphere2DDataConfig->grid_number_elements; i++)
		t[i] = 0;

#else

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] = 0;

#endif
}


void DataGrid::grid_setValue(
		std::complex<double> &i_value
)
{
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	for (std::size_t i = 0; i < sphere2DDataConfig->grid_number_elements; i++)
		grid_space_data[i] = i_value;
}



void DataGrid::grid_setValue(
		int i_lon_idx,
		int i_lat_idx,
		std::complex<double> &i_value
)
{
#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS
	grid_space_data[i_lon_idx*sphere2DDataConfig->grid_num_lat + i_lat_idx] = i_value;
#else
	grid_space_data[i_lat_idx*sphere2DDataConfig->grid_num_lon + i_lon_idx] = i_value;
#endif
}


double DataGrid::grid_reduce_max(
		const Sphere2DComplex::DataGrid &i_data
)
{
	check(i_data.sphere2DDataConfig);

	double error = -1;

	for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++)
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

	for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++)
	{
		std::complex<double> &d = grid_space_data[j];
		error += d.real()*d.real() + d.imag()*d.imag();
	}

	return std::sqrt(error / (double)sphere2DDataConfig->grid_number_elements);
}


double DataGrid::grid_reduce_max_abs()	const
{
	double error = -1;

	for (std::size_t j = 0; j < sphere2DDataConfig->grid_number_elements; j++)
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


}}}
