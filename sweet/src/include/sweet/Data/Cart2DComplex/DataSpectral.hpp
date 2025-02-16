/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATASPECTRAL_HPP

#include <complex>

#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>

namespace sweet {
namespace Data {
namespace Cart2DComplex {


/*!
 * \brief Data storage class for complex value in spectral Fourier space for a 2D Cartesian grid with complex values
 */
class DataSpectral
{
public:
	const Cart2D::Config *cart2DDataConfig;

	typedef std::complex<double> Tcomplex;

public:
	std::complex<double> *spectral_space_data;

	std::complex<double>& operator[](std::size_t i);

	const std::complex<double>& operator[](std::size_t i)	const;


public:
	DataSpectral();


public:
	DataSpectral(
			const Cart2D::Config *i_cart2DDataConfig
	);

public:
	DataSpectral(
			const Cart2D::Config &i_cart2DDataConfig
	);


public:
	DataSpectral(
			const DataSpectral &i_cart2d_data
	);

public:
	DataSpectral(
			DataSpectral &&i_data
	);


	DataSpectral(
			const DataGrid &i_cart2d_data
	);


	/*!
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	 void check_cart2DDataConfig_identical_res(
			 const Cart2D::Config *i_cart2DDataConfig
	) const;



public:
	DataSpectral& operator=(
			const DataSpectral &i_cart2d_data
	);


public:
	DataSpectral& operator=(
			DataSpectral &&i_cart2d_data
	);




	DataGrid toGrid()	const;

	DataGrid getCart2DDataGridComplex()	const;


	DataSpectral operator+(
			const DataSpectral &i_cart2d_data
	)	const;



	DataSpectral& operator+=(
			const DataSpectral &i_cart2d_data
	);


	DataSpectral& operator-=(
			const DataSpectral &i_cart2d_data
	);



	DataSpectral operator-(
			const DataSpectral &i_cart2d_data
	)	const;


	DataSpectral operator-()	const;


	DataSpectral operator*(
			const double i_value
	)	const;

	const DataSpectral& operator*=(
			const std::complex<double> &i_value
	)	const;

	DataSpectral operator*(
			const DataSpectral &i_cart2d_data
	)	const;

	DataSpectral operator/(
			double i_value
	)	const;

	DataSpectral operator/(
			const DataSpectral &i_cart2d_data
	)	const;

	const DataSpectral& operator/=(
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


	/*!
	 * Apply a linear operator given by this class to the input data array.
	 */

	DataSpectral operator()(
			const DataSpectral &i_array_data
	)	const;


public:
	void setup_if_required(
		const Cart2D::Config *i_cart2DDataConfig
	);

public:
	void setup(
			const Cart2D::Config *i_cart2dConfig
	);

public:
	~DataSpectral();


	/*!
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

	void spectral_update_lambda_modes(
			std::function<void(int,int,std::complex<double>&)> i_lambda	//!< lambda function to return value for lat/mu
	);



	const std::complex<double>& spectral_get(
			int in,
			int im
	)	const;


	/*
	 * Set all values to zero in spectral space
	 */
	void spectral_setZero();

	void spectral_set(
			int i_n,
			int i_m,
			double i_real,
			double i_imag
	)	const;

	void spectral_set(
			int i_n,
			int i_m,
			std::complex<double> i_data
	)	const;

	/*!
	 * Add scalar to all spectral modes
	 */
	DataSpectral spectral_addScalarAll(
			const double &i_value
	)	const;


	/*!
	 * Add scalar to all spectral modes
	 */
	DataSpectral spectral_addScalarAll(
			const std::complex<double> &i_value
	)	const;

	/*!
	 * Invert the application of a linear operator in spectral space.
	 * The operator is given in i_array_data
	 */
	DataSpectral spectral_div_element_wise(
			const DataSpectral &i_array_data	//!< operator
	)	const;


	/*!
	 * Return Cart2D Array with all spectral coefficients a+bi --> 1/(a+bi)
	 */
	DataSpectral spectral_invert()	const;


	void spectral_print(
			int i_precision = 8
	)	const;


	void spectral_zeroAliasingModes()	const;

	/*!
	 * Test for real values only in physical space
	 * by checking for certain symmetries
	 */
	void test_realphysical()	const;

	/*
	 * Setup spectral sphere2D data based on data in physical space
	 */
	void loadCart2DDataGrid(
			const DataGrid &i_cart2DDataGrid
	);


	void print_spectralData()	const;


	/*!
	 * print spectral data and zero out values which are numerically close to zero
	 */
	void print_spectralData_zeroNumZero(double i_zero_threshold = 1e-13)	const;

};

#include <sweet/Data/Cart2DComplex/DataSpectral_Operators_inc.hpp>

}}}

#include <sweet/Data/Cart2DComplex/DataSpectral_Operators_inc.hpp>

#endif
