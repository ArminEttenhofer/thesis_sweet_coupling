/*
 * Sphere2DData.hpp
 *
 *  Created on: 9 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_DATASPECTRAL_HPP

#include <sweet/Data/Sphere2D/DataGrid.hpp>
//#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
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
#include <sweet/Memory/parmemcpy.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Fatal.hpp>


namespace sweet {
namespace Data {
namespace Sphere2DComplex {
	class DataGrid;
	class DataSpectral;
}}}


namespace sweet {
namespace Data {
namespace Sphere2D {

class DataGrid;

class DataSpectral
{
//	friend class Sphere2DComplex::DataSpectral;

	typedef std::complex<double> TComplex;

public:
	const Sphere2D::Config *sphere2DDataConfig = nullptr;

public:
	std::complex<double> *spectral_space_data = nullptr;

public:
	std::complex<double>& operator[](std::size_t i);

	const std::complex<double>& operator[](std::size_t i)	const;

public:
	void swap(
			Sphere2D::DataSpectral &i_sphere2DData
	);

public:
	void swapWithConfig(
			Sphere2D::DataSpectral &i_sphere2DData
	);


public:
	DataSpectral(
			const Sphere2D::Config *i_sphere2DDataConfig
	);

public:
	DataSpectral(
			const Sphere2D::Config &i_sphere2DDataConfig
	);


public:
	DataSpectral(
			const Sphere2D::Config *i_sphere2DDataConfig,
			const double &i_value
	);

public:
	DataSpectral();



public:
	DataSpectral(
			const Sphere2D::DataSpectral &i_sph_data
	);



public:
	DataSpectral(
			Sphere2D::DataSpectral &&i_sph_data
	);




	/*!
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	 void check(
			const Sphere2D::Config *i_sphere2DDataConfig
	)	const;

public:
	DataSpectral& operator=(
			const Sphere2D::DataSpectral &i_sph_data
	);




	/*1
	 * This function implements copying the spectral data only.
	 *
	 * This becomes handy if coping with data which should be only transformed without dealiasing.
	 */
public:
	DataSpectral& load_nodealiasing(
			const Sphere2D::DataSpectral &i_sph_data		//!< data to be converted to sphere2DDataConfig_nodealiasing
	);

public:
	DataSpectral& operator=(
			Sphere2D::DataSpectral &&i_sph_data
	);

public:
	DataSpectral spectral_returnWithDifferentModes(
			const Sphere2D::Config *i_sphere2DDataConfig
	)	const;

	/*!
	 * Setup spectral sphere2D data based on data in physical space
	 */
public:
	void loadSphere2DDataGrid(
			const Sphere2D::DataGrid &i_sphere2DDataGrid
	);

	/*!
	 * Return the data converted to physical space
	 */
public:
	Sphere2D::DataGrid getSphere2DDataGrid()	const;

	/*!
	 * Return the data converted to physical space
	 *
	 * Alias for "getSphere2DDataGrid"
	 */
public:
	Sphere2D::DataGrid toGrid()	const;


public:
	Sphere2DComplex::DataGrid getSphere2DDataGridComplex()	const;


public:
	DataSpectral(
			const Sphere2D::DataGrid &i_sphere2d_data_physical
	);


public:
	DataSpectral operator+(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

	DataSpectral& operator+=(
			const Sphere2D::DataSpectral &i_sph_data
	);

	DataSpectral& operator-=(
			const Sphere2D::DataSpectral &i_sph_data
	);

	DataSpectral operator-(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

	DataSpectral operator-()	const;

	DataSpectral operator*(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

	DataSpectral operator/(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

	DataSpectral operator*(
			double i_value
	)	const;

	const Sphere2D::DataSpectral& operator*=(
			double i_value
	)	const;

	const Sphere2D::DataSpectral& operator/=(
			double i_value
	)	const;

	const Sphere2D::DataSpectral& operator*=(
			const std::complex<double> &i_value
	)	const;

	DataSpectral operator/(
			double i_value
	)	const;

	DataSpectral operator+(
			double i_value
	)	const;

	DataSpectral spectralElementwiseMultiplication(
			const Sphere2D::DataSpectral &i_sph_data
	)	const;

public:
	DataSpectral operator_scalar_sub_this(
			double i_value
	)	const;

	const Sphere2D::DataSpectral& operator+=(
			double i_value
	)	const;

	DataSpectral operator-(
			double i_value
	)	const;


	DataSpectral& operator-=(
			double i_value
	);

public:
	bool setup(
		const Sphere2D::Config *i_sphere2DDataConfig
	);

public:
	bool isSetup();

public:
	bool setup(
		const Sphere2D::Config &i_sphere2DDataConfig
	);

public:
	bool setup(
		const Sphere2D::Config *i_sphere2DDataConfig,
		double i_value
	);

private:
	bool alloc_data();

public:
	void setup_if_required(
		const Sphere2D::Config *i_sphere2DDataConfig
	);

public:
	void clear();

public:
	~DataSpectral();


	/*!
	 * Solve a Helmholtz problem given by
	 *
	 * \f$$
	 * 	(a + b D^2) x = rhs
	 * \f$$
	 */
	DataSpectral spectral_solve_helmholtz(
			const double &i_a,
			const double &i_b,
			double r
	);


	/*!
	 * Solve a Helmholtz problem given by
	 *
	 * \f$$
	 * 	(a + b0 D^2 + b1 D^4 + b2 D^6 + b3 D^8) x = rhs
	 * \f$$
	 */
	DataSpectral spectral_solve_helmholtz_higher_order(
			const double & a,
			const std::array<double, 4> & b,
			double r
	);


	/*!
	 * Solve a Laplace problem given by
	 *
	 * \f$$
	 * 	(D^2) x = rhs
	 * \f$$
	 */
public:
	DataSpectral spectral_solve_laplace(
			double r
	);

	/*!
	 * Truncate modes which are not representable in spectral space
	 */
public:
	const Sphere2D::DataSpectral& spectral_truncate()	const;

public:
	void spectral_update_lambda(
			std::function<void(int,int,TComplex&)> i_lambda
	);

public:
	const std::complex<double>& spectral_get_DEPRECATED(
			int i_n,
			int i_m
	)	const;

public:
	const std::complex<double>& spectral_get_(
			int i_n,
			int i_m
	)	const;

public:
	void spectral_set(
			int i_n,
			int i_m,
			const std::complex<double> &i_data
	)	const;


	/*!
	 * Set all values to zero
	 */
public:
	void spectral_setZero();


	/*!
	 * Set all values to the given value
	 */
public:
	void spectral_setValue(
			const std::complex<double> &i_value
	);

	/*!
	 * Add a constant in physical space by adding the corresponding value in spectral space
	 */
public:
	void spectral_add_grid_constant(
			double i_value
	)	const;

	/*!
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
public:
	std::complex<double> spectral_reduce_sum_quad_increasing()	const;

	/*!
	 * return the sum of all values, use quad precision for reduction
	 */
public:
	std::complex<double> spectral_reduce_sum_quad()	const;

	/*!
	 * return the sum of squares of all values, use quad precision for reduction
	 * Important: Since m=0  modes appear only once and m>0 appear twice in the full spectrum
	 */
public:
	double spectral_reduce_sum_sqr_quad()	const;

	/*!
	 * Return the minimum value
	 */
public:
	std::complex<double> spectral_reduce_min()	const;


	/*!
	 * Return the max value
	 */
public:
	std::complex<double> spectral_reduce_max()	const;


	/*!
	 * Return the max abs value
	 */
public:
	double spectral_reduce_max_abs()	const;


	/*!
	 * Return the max abs value for the first rnorm spectral coefficients
	 */
public:
	double spectral_reduce_max_abs(std::size_t rnorm)	const;


	/*!
	 * Return the minimum abs value
	 */
public:
	double spectral_reduce_min_abs()	const;

public:
	bool spectral_reduce_is_any_nan_or_inf()	const;

public:
	bool spectral_is_first_nan_or_inf()	const;


public:
	void spectral_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const;

public:
	void spectral_structure_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const;


	/*!
	 * Write spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 *
	 * Per default, a tab separator is used in each line to separate the values.
	 */
public:
	bool file_spectral_saveData_ascii(
			const char *i_filename,		//!< Name of file to store data to
			char i_separator = '\t',	//!< separator to use for each line
			int i_precision = 12,		//!< number of floating point digits
			int dimension = 2			//!< store 1D or 2D
	)	const;


public:
 	void spectrum_file_write(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const;


public:
	void spectrum_abs_file_write_line(
			const std::string &i_filename,
			const char *i_title = "",
			const double i_time = 0.0,
			int i_precision = 20,
			double i_abs_threshold = -1,
			int i_reduce_mode_factor = 1
	)	const;

public:
	void spectrum_phase_file_write_line(
			const std::string &i_filename,
			const char *i_title = "",
			const double i_time = 0.0,
			int i_precision = 20,
			double i_abs_threshold = -1,
			int i_reduce_mode_factor = 1
	)	const;

  	/*!
  	 * Write the spectral data to a file in binary format
  	 */
public:
  	void file_write_binary_spectral(
			const std::string &i_filename
	)	const;

public:
	void loadDataFromFile(
		const std::string i_filename,
		const std::string i_filemode
	);

public:
	void file_read_csv_grid(
			const std::string &i_filename
	);

public:
  	void file_read_binary_spectral(
			const std::string &i_filename
	);

public:
	void normalize(
			const std::string &normalization = ""
	);


	/*!
	 * Interpolate from a finer mesh. Remove highest frequency modes.
	 */
public:
	DataSpectral restrict(
			const Sphere2D::DataSpectral &i_array_data
	);

	/*!
	 * Interpolate from a coarser mesh. Pad zeros corresponding to highest frequency modes.
	 */
public:
	DataSpectral pad_zeros(
			const Sphere2D::DataSpectral &i_array_data
	);
};


}}}



/*!
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
sweet::Data::Sphere2D::DataSpectral operator*(
		double i_value,
		const sweet::Data::Sphere2D::DataSpectral &i_array_data
);



/*!
 * operator to support operations such as:
 *
 * 1.5 + arrayData
 *
 */
sweet::Data::Sphere2D::DataSpectral operator+(
		double i_value,
		const sweet::Data::Sphere2D::DataSpectral &i_array_data
);


/*!
 * operator to support operations such as:
 *
 * 1.5 - arrayData
 *
 */
sweet::Data::Sphere2D::DataSpectral operator-(
		double i_value,
		const sweet::Data::Sphere2D::DataSpectral &i_array_data
);

#endif
