/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_DATASPECTRAL_HPP

#include <complex>
#include <functional>
#include <array>
#include <string.h>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <functional>

#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Memory/parmemcpy.hpp>


namespace sweet {
namespace Data {
namespace Sphere2DComplex {



class DataSpectral
{
public:
	const Sphere2D::Config *sphere2DDataConfig;

	typedef std::complex<double> Tcomplex;

public:
	std::complex<double> *spectral_space_data;

	std::complex<double>& operator[](std::size_t i);

	const std::complex<double>& operator[](std::size_t i)	const;

public:
	void swap(
		Sphere2DComplex::DataSpectral &i_sphere2DData
	);



public:
	DataSpectral();


public:
	DataSpectral(
			const Sphere2D::Config *i_sphere2DDataConfig
	);


public:
	DataSpectral(
			const Sphere2DComplex::DataSpectral &i_sph_data
	);



public:
	DataSpectral(
			Sphere2DComplex::DataSpectral &&i_data
	);



	DataSpectral(
			const Sphere2DComplex::DataGrid &i_sph_data
	);


	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	 void check_sphere2DDataConfig_identical_res(const Sphere2D::Config *i_sphere2DDataConfig)	const;



public:
	DataSpectral& operator=(
			const Sphere2DComplex::DataSpectral &i_sph_data
	);


public:
	DataSpectral& operator=(
			Sphere2DComplex::DataSpectral &&i_sph_data
	);

	DataSpectral spectral_returnWithTruncatedModes(
			const Sphere2D::Config *i_sphere2DDataConfigTargetTruncation
	)	const;

public:
	DataSpectral spectral_returnWithDifferentModes(
			const Sphere2D::Config *i_sphere2DDataConfigNew
	)	const;

	/*
	 * Setup spectral sphere2D data based on data in physical space
	 */
	void loadSphere2DDataGrid(
			const Sphere2DComplex::DataGrid &i_sphere2DDataGrid
	);



	Sphere2DComplex::DataGrid toGrid()	const;

	Sphere2DComplex::DataGrid getSphere2DDataGridComplex()	const;

	DataSpectral operator+(
			const Sphere2DComplex::DataSpectral &i_sph_data
	) const;


	DataSpectral& operator+=(
			const Sphere2DComplex::DataSpectral &i_sph_data
	);

	DataSpectral& operator-=(
			const Sphere2DComplex::DataSpectral &i_sph_data
	);


	DataSpectral operator-(
			const Sphere2DComplex::DataSpectral &i_sph_data
	)	const;


	DataSpectral operator-()	const;



	DataSpectral operator*(
			const double i_value
	)	const;

	const Sphere2DComplex::DataSpectral& operator*=(
			const std::complex<double> &i_value
	)	const;


	DataSpectral operator/(
			double i_value
	)	const;


	const Sphere2DComplex::DataSpectral& operator/=(
			const std::complex<double> &i_value
	)	const;


	DataSpectral operator+(
			double i_value
	)	const;


	DataSpectral& operator+=(
			double i_value
	);


	DataSpectral operator+(
			const std::complex<double> &i_value
	)	const;

	DataSpectral operator-(
			const std::complex<double> &i_value
	)	const;


	DataSpectral operator*(
			const std::complex<double> &i_value
	)	const;



public:
	void setup_if_required(
		const Sphere2D::Config *i_sphere2DDataConfig
	);

public:
	void setup(
			const Sphere2D::Config *i_sphere2DConfig
	);

public:
	void clear();

public:
	~DataSpectral();


	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */

	DataSpectral spectral_solve_helmholtz(
			const std::complex<double> &i_a,
			const std::complex<double> &i_b,
			double r
	)	const;




	void spectral_update_lambda(
			std::function<void(int,int,Tcomplex&)> i_lambda
	);



	const std::complex<double>& spectral_get_(
			int in,
			int im
	)	const;


	/*
	 * Set all values to zero in spectral space
	 */
	void spectral_setZero();


	void spectral_print(
			int i_precision = 8
	)	const;
};





}}}


sweet::Data::Sphere2DComplex::DataSpectral operator*(
		const double i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
);



sweet::Data::Sphere2DComplex::DataSpectral operator*(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
);


/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */
sweet::Data::Sphere2DComplex::DataSpectral operator+(
		const double i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
);


sweet::Data::Sphere2DComplex::DataSpectral operator+(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
);


sweet::Data::Sphere2DComplex::DataSpectral operator-(
		const std::complex<double> &i_value,
		const sweet::Data::Sphere2DComplex::DataSpectral &i_array_data
);


#endif
