/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_DATAGRID_HPP

#include "DataGrid.hpp"
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <complex>


namespace sweet {
namespace Data {
namespace Sphere2DComplex {


class DataGrid
{
	friend class DataSpectral;

public:
	const Sphere2D::Config *sphere2DDataConfig;

public:
	std::complex<double> *grid_space_data;


	void swap(
			Sphere2DComplex::DataGrid &i_sphere2DData
	);

public:
	DataGrid(
			const Sphere2D::Config *i_sphere2DDataConfig
	);


public:
	DataGrid();


public:
	DataGrid(
			const Sphere2DComplex::DataGrid &i_sph_data
	);


public:
	DataGrid(
			const Sphere2D::DataGrid &i_sph_data
	);


	/*
	 * load real and imaginary data from physical arrays
	 */
	void loadRealImag(
			const Sphere2D::DataGrid &i_re,
			const Sphere2D::DataGrid &i_im
	);


	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const Sphere2D::Config *i_sphere2DDataConfig
	)	const;



public:
	Sphere2DComplex::DataGrid& operator=(
			const Sphere2DComplex::DataGrid &i_sph_data
	);


public:
	Sphere2DComplex::DataGrid& operator=(
			Sphere2DComplex::DataGrid &&i_sph_data
	);



public:
	Sphere2DComplex::DataGrid& operator=(
			const Sphere2D::DataGrid &i_sph_data
	);

	Sphere2DComplex::DataGrid operator+(
			const Sphere2DComplex::DataGrid &i_sph_data
	)	const;


	Sphere2DComplex::DataGrid& operator+=(
			const Sphere2DComplex::DataGrid &i_sph_data
	);


	Sphere2DComplex::DataGrid& operator-=(
			const Sphere2DComplex::DataGrid &i_sph_data
	);



	Sphere2DComplex::DataGrid operator-(
			const Sphere2DComplex::DataGrid &i_sph_data
	)	const;


	Sphere2DComplex::DataGrid operator-();


	Sphere2DComplex::DataGrid operator*(
			const Sphere2DComplex::DataGrid &i_sph_data
	)	const;


	Sphere2DComplex::DataGrid operator/(
			const Sphere2DComplex::DataGrid &i_sph_data
	)	const;


	Sphere2DComplex::DataGrid operator*(
			const double i_value
	)	const;


	Sphere2DComplex::DataGrid operator*(
			const std::complex<double> &i_value
	)	const;

	const Sphere2DComplex::DataGrid& operator*=(
			const double i_value
	)	const;


	Sphere2DComplex::DataGrid operator/(
			double i_value
	)	const;


	Sphere2DComplex::DataGrid operator+(
			const std::complex<double> &i_value
	)	const;

	Sphere2DComplex::DataGrid operator-(
			const std::complex<double> &i_value
	)	const;


public:
	void setup(
			const Sphere2D::Config *i_sphere2DDataConfig
	);


public:
	void setup_if_required(
			const Sphere2D::Config *i_sphere2DDataConfig
	);


public:
	~DataGrid();



	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
	void grid_update_lambda(
			std::function<void(double,double,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);

	void grid_update_lambda_array(
			std::function<void(int,int,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);

	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude sin(phi) \in [-1;1])
	 */
	void grid_update_lambda_gaussian_grid(
			std::function<void(double,double,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);


	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters:
	 *   (longitude \in [0;2*pi], Cogaussian latitude cos(phi) \in [0;1])
	 */
	void grid_update_lambda_cogaussian_grid(
			std::function<void(double,double,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);

	void grid_update_lambda_sinphi_grid(
			std::function<void(double,double,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);

	void grid_update_lambda_cosphi_grid(
			std::function<void(double,double,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
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
			int i_lon_idx,
			int i_lat_idx,
			std::complex<double> &i_value
	);

	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double grid_reduce_max(
			const Sphere2DComplex::DataGrid &i_sph_data
	);

	double grid_reduce_rms();



	/**
	 * Return the maximum error norm
	 */
	double grid_reduce_max_abs()	const;
};

}}}


/**
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
sweet::Data::Sphere2DComplex::DataGrid operator*(
		const double i_value,
		const sweet::Data::Sphere2DComplex::DataGrid &i_array_data
)
{
	return ((sweet::Data::Sphere2DComplex::DataGrid&)i_array_data)*i_value;
}


inline
static
sweet::Data::Sphere2DComplex::DataGrid operator*(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataGrid &i_array_data
)
{
	return ((sweet::Data::Sphere2DComplex::DataGrid&)i_array_data)*i_value;
}


#endif
