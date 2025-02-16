/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_CART2D_DATASPECTRAL_HPP

#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <complex>
#include <functional>

#if !SWEET_USE_LIBFFT
#error "LIBFFT not activated, but spectral cart2d data compiled in"
#endif


namespace sweet {
namespace Data {
namespace Cart2D {



/*!
 * \brief Data storage class for complex value in spectral Fourier space for a 2D Cartesian grid
 */
class DataSpectral
{
	typedef std::complex<double> Tcomplex;

public:
	const Config *cart2DDataConfig = nullptr;

public:
	std::complex<double> *spectral_space_data = nullptr;

	std::complex<double>& operator[](std::size_t i);

	const std::complex<double>& operator[](std::size_t i)	const;

	void swap(
			DataSpectral &i_cart2DData
	);


public:
	DataSpectral(
			const Config *i_cart2DDataConfig
	);


public:
	DataSpectral(
			const Config *i_cart2DDataConfig,
			const std::complex<double> &i_value
	);



public:
	DataSpectral(
			const Config *i_cart2DDataConfig,
			double &i_value
	);



	/*!
	 * Without setup where we need to call setup(...) later on
	 */
public:
	DataSpectral();


	/*!
	 * dummy initialization by handing over an unused integer
	 */
public:
	DataSpectral(int i);


public:
	DataSpectral(
			const DataSpectral &i_cart2d_data
	);



public:
	DataSpectral(
			DataSpectral &&i_cart2d_data
	);


	/*!
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	void _validateRes(
			const Config *i_cart2DDataConfig
	) const;


public:
	DataSpectral& operator=(
			const DataSpectral &i_cart2d_data
	);

public:
	/*!
	 * assignment operator
	 */
	DataSpectral &operator=(double i_value);

public:
	/*!
	 * assignment operator
	 */
	DataSpectral &operator=(int i_value);


public:
	void spectral_zeroAliasingModes();

public:
	void spectral_debugCheckForZeroAliasingModes()	const;


	/*!
	 * This function implements copying the spectral data only.
	 *
	 * This becomes handy if coping with data which should be only transformed without dealiasing.
	 */
public:
	DataSpectral& load_nodealiasing(
			const DataSpectral &i_cart2d_data		//!< data to be converted to cart2DDataConfig_nodealiasing
	);


public:
	DataSpectral& operator=(
			DataSpectral &&i_cart2d_data
	);


	DataSpectral spectral_returnWithDifferentModes(
			const Config &i_cart2DDataConfig
	)	const;

public:
	DataSpectral spectral_returnWithDifferentModes(
			const Config *i_cart2DDataConfig
	)	const;

	/*!
	 * Return Cart2D Array with all spectral coefficients a+bi --> 1/(a+bi)
	 */
	DataSpectral spectral_invert()	const;



	/*!
	 * Setup spectral sphere2D data based on data in physical space
	 */
	void loadCart2DDataGrid(
		const DataGrid &i_cart2DDataGrid
	);


	/*
	 * Return the data converted to physical space
	 */
	DataGrid getCart2DDataGrid()	const;


	/*
	 * Return the data converted to physical space
	 *
	 * alias for "getCart2DDataGrid"
	 */
	DataGrid toGrid()	const;


	Cart2DComplex::DataGrid getCart2DDataGridComplex()	const;


	DataSpectral(
			const DataGrid &i_cart2d_data_physical
	);



private:
	void _main_setup(
		const Config *i_cart2DDataConfig
	);


public:
	void setup(
		const Config *i_cart2DDataConfig
	);


	/*
	 * Wrapper for main setup
	 */
public:
	void setup(
		const Config &i_cart2DDataConfig
	);

public:
	void setup(
		const Config *i_cart2DDataConfig,
		double i_value
	);

private:
	void alloc_data();


public:
	void setup_if_required(
		const Config *i_cart2DDataConfig
	);


public:
	void clear();


public:
	~DataSpectral();


	DataSpectral operator+(
			const DataSpectral &i_cart2d_data
	) const;


	DataSpectral& operator+=(
			const DataSpectral &i_cart2d_data
	);


	DataSpectral& operator-=(
			const DataSpectral &i_cart2d_data
	);



	DataSpectral operator-(
			const DataSpectral &i_cart2d_data
	)	const;


	DataSpectral operator-(
			const DataGrid &i_cart2d_data_physical
	)	const;



	DataSpectral operator-()	const;



	DataSpectral operator*(
			const DataSpectral &i_cart2d_data
	)	const;


	DataGrid multiplication_grid_space(
				const DataGrid &i_a,
				const DataGrid &i_b
	) const;

	DataSpectral operator/(
			const DataSpectral &i_cart2d_data
	)	const;


	DataSpectral operator*(
			double i_value
	)	const;


	const DataSpectral& operator*=(
			double i_value
	)	const;


	const DataSpectral& operator/=(
			double i_value
	)	const;


	const DataSpectral& operator*=(
			const std::complex<double> &i_value
	)	const;


	DataSpectral operator/(
			double i_value
	)	const;


	DataSpectral operator+(
			double i_value
	)	const;


	DataSpectral operator_scalar_sub_this(
			double i_value
	)	const;


	const DataSpectral& operator+=(
			double i_value
	)	const;

	DataSpectral operator-(
			double i_value
	)	const;


	DataSpectral& operator-=(
			double i_value
	);

public:
	/*!
	 * Add scalar to all spectral modes
	 */
	DataSpectral spectral_addScalarAll(
			const double &i_value
	)	const;


	/*!
	 * Invert the application of a linear operator in spectral space.
	 * The operator is given in i_array_data
	 *
	 * Solving for phi in
	 *
	 * 		\Delta^2 phi = delta_phi
	 *
	 * can be accomplished by
	 *
	 * 		phi = delta_phi.spectral_div_element_wise(\Delta^2)
	 */
	DataSpectral spectral_div_element_wise(
			const DataSpectral &i_array_data	//!< operator
	)	const;


	/*!
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	DataSpectral spectral_solve_helmholtz(
			const double &i_a,
			const double &i_b,
			double r
	);


	/*!
	 * Solve a Laplace problem given by
	 *
	 * (D^2) x = rhs
	 */
	DataSpectral spectral_solve_laplace(
			double r
	);

	/*!
	 * Truncate modes which are not representable in spectral space
	 */
	const DataSpectral& spectral_truncate()	const;


	void spectral_update_lambda(
			std::function<void(int,int,Tcomplex&)> i_lambda
	);


	const std::complex<double>& spectral_get(
			int i_n,
			int i_m
	)	const;

	void spectral_set(
			int i_n,
			int i_m,
			const std::complex<double> &i_data
	)	const;


	/*
	 * Set all values to zero
	 */
	void spectral_setZero();


	/*
	 * Set all values to value
	 */
	void spectral_setValue(
			const std::complex<double> &i_value
	);


	/*
	 * Add a constant in physical space by adding the corresponding value in spectral space
	 */
	void spectral_add_grid_constant(
			double i_value
	)	const;


	/*!
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	std::complex<double> spectral_reduce_sum_quad_increasing()	const;


	/*!
	 * return the sum of all values, use quad precision for reduction
	 */
	std::complex<double> spectral_reduce_sum_quad()	const;

	/*!
	 * return the sum of squares of all values, use quad precision for reduction
	 * Important: Since m=0  modes appear only once and m>0 appear twice in the full spectrum
	 */
	double spectral_reduce_sum_sqr_quad()	const;

	/*!
	 * reduce to sum square in spectrum
	 */
	double spectral_reduce_sum_sq();


	/*!
	 * Return the minimum value
	 */
	std::complex<double> spectral_reduce_min()	const;

	/*!
	 * Return the max value
	 */
	std::complex<double> spectral_reduce_max()	const;


	/*!
	 * Return the max abs value
	 */
	double spectral_reduce_max_abs()	const;


	/*!
	 * Return the max abs value for the first rnorm spectral coefficients
	 */
	double spectral_reduce_max_abs(std::size_t rnorm)	const;



	/*!
	 * Return the minimum abs value
	 */
	double spectral_reduce_min_abs()	const;

	/*!
	 * reduce to root mean square in spectrum
	 */
	double spectral_reduce_rms();


	bool spectral_reduce_is_any_nan_or_inf()	const;

	bool spectral_is_first_nan_or_inf()	const;


	void spectral_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const;


	void print_spectralData()	const;

	/*!
	 * print spectral data and zero out values which are numerically close to zero
	 */
	void print_spectralData_zeroNumZero(double i_zero_threshold = 1e-13)	const;

	void print_spectralIndex()	const;

	void print_spectralNonZero()	const;


	void spectral_structure_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const;

  
  	void spectrum_file_write(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const;


	void spectrum_abs_file_write_line(
			const std::string &i_filename,
			const char *i_title = "",
			const double i_time = 0.0,
			int i_precision = 20,
			double i_abs_threshold = -1,
			int i_reduce_mode_factor = 1
	)	const;


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
  	void file_write_binary_spectral(
			const std::string &i_filename
	)	const;

	void loadDataFromFile(
		const std::string i_filename,
		const std::string i_filemode
	);

	void file_read_csv_grid(
			const std::string &i_filename
	);

	void file_read_binary_spectral(
			const std::string &i_filename
	);


	/*!
	 * Write spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_spectral_saveData_ascii(
			const char *i_filename,		//!< Name of file to store data to
			char i_separator = '\t',	//!< separator to use for each line
			int i_precision = 12,		//!< number of floating point digits
			int dimension = 2			//!< store 1D or 2D
	)	const;

	/*!
	 * Write amplitude of spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_spectral_abs_saveData_ascii(
			const char *i_filename,		//!< Name of file to store data to
			char i_separator = '\t',	//!< separator to use for each line
			int i_precision = 12,		//!< number of floating point digits
			int dimension = 2			//!< store 1D or 2D
	)	const;

	/*!
	 * Write spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_spectral_arg_saveData_ascii(
			const char *i_filename,		//!< Name of file to store data to
			char i_separator = '\t',	//!< separator to use for each line
			int i_precision = 12,		//!< number of floating point digits
			int dimension = 2			//!< store 1D or 2D
	)	const;


	/*!
	 * Return average which is given by the first mode
	 */
	double get_average()	const;

	double spectral_return_amplitude(
			std::size_t j,
			std::size_t i
	) const;

	double spectral_return_phase(
			std::size_t j,
			std::size_t i
	) const;


	void normalize(
			const std::string &normalization = ""
	);


	/*!
	 * Apply a linear operator given by this class to the input data array.
	 */
	DataSpectral operator()(
			const DataSpectral &i_array_data
	)	const;


	/*!
	 * Interpolate from a finer mesh. Remove highest frequency modes.
	 */
	DataSpectral restrict(
			const DataSpectral &i_array_data
	);


	/*!
	 * Interpolate from a coarser mesh. Pad zeros corresponding to highest frequency modes.
	 */
	DataSpectral pad_zeros(
			const DataSpectral &i_array_data
	);

};

#include <sweet/Data/Cart2D/DataSpectral_Operators_inc.hpp>

}}}

#include <sweet/Data/Cart2D/DataSpectral_Operators_inc.hpp>

#endif
