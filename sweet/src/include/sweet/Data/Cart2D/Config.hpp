/*
 * Cart2DDataConfig.hpp
 *
 *  Created on: 17 Oct 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_CONFIG_HPP
#define INCLUDE_SWEET_DATA_CART2D_CONFIG_HPP



#include <iostream>
#include <complex>

#include <fftw3.h>

#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>

#include <sweet/Tools/TransformationPlans.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>


namespace sweet {
namespace Data {
namespace Cart2D {


/*!
 * \brief Configuration class for Cart2D Data
 *
 * This is used to setup the transformation plans and keep all configuration in one place.
 */
class Config
{
public:
	//! Error handling
	Error::Base error;

	int ndim = 2;

	//! Grid resolution: Number of cells in each dimension
	std::size_t grid_res[2];

	//! Size of physical storage, is identical to physical resolution
	std::size_t grid_data_size[2];

	//! number of real valued data
	std::size_t grid_number_elements;

	//! Different strategies to cope with transformation plans
	sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE reuse_spectral_transformation_plans;


public:

#if SWEET_USE_LIBFFT

	/*!
	 * Number of spectral modes
	 *
	 * This is not related to the storage size!
	 *
	 * Also, this includes the modes which cannot be represented according
	 * to the Shannon-Nyquist theorem!
	 */
	std::size_t spectral_modes[2];

	/*!
	 * Number of spectral modes which are realted to real, representable, modes.
	 */
	std::size_t spectral_representable_modes[2];

	//! allocated size for spectral data for each modes
	//! This storage size is for the real-to-complex transformation
	std::size_t spectral_data_size[2];

	//! iteration ranges for updating data in spectrum
	//! 1st index (left): the id of the range,
	//! 2nd index (middle): dimension the range,
	//! 3rd index (last one): start and end (exclusive) index
	std::size_t spectral_data_iteration_ranges[2][2][2];

	//! total number of complex-valued data elements in spectral space
	std::size_t spectral_array_data_number_of_elements;

private:
	/*
	 * FFTW related stuff
	 */
	fftw_plan	_fftw_plan_forward;
	fftw_plan	_fftw_plan_backward;

	//! FFTW scaling related stuff for backward transformation
	//! WARNING: FFTW doesn't implement a symmetric FFTW
	//! We only to the rescaling for the backward transformation
	double fftw_backward_scale_factor;


public:
	//! allocated size for spectral data in case of complex data in physical space
	std::size_t spectral_complex_data_size[2];

	//! total number of elements in spectrum
	std::size_t spectral_complex_array_data_number_of_elements;

	//! iteration range for complex valued space
	//! 1st index (left): the id of the range,
	//! 2nd index (middle): dimension the range,
	//! 3rd index (last one): start and end (exclusive) index
	std::size_t spectral_complex_ranges[4][2][2];

	//! Forward FFTW plan
	fftw_plan	_fftw_plan_complex_forward;

	//! Backward FFTW plan
	fftw_plan	_fftw_plan_complex_backward;

#endif

	/*!
	 * True if it's already initialized
	 */
private:
	bool fftw_initialized;


public:
	Config()
	{
		grid_res[0] = 0;
		grid_res[1] = 0;

#if SWEET_USE_LIBFFT
		spectral_modes[0] = 0;
		spectral_modes[1] = 0;
#endif

		fftw_initialized = false;
	}


	/*!
	 * Return a unique string for the current configuration
	 */
public:
	std::string getUniqueIDString()	const
	{
		return getConfigInformationString();
	}

	/*!
	 * Get a human readable information string about the current config
	 */
	std::string getConfigInformationString()	const
	{
		std::ostringstream buf;
		buf <<

#if SWEET_USE_LIBFFT
				"M" << spectral_modes[0] << "," << spectral_modes[1] << "_" <<
#endif
				"N" << grid_res[0] << "," << grid_res[1];

		return buf.str();
	}


	/*!
	 * Print various information about this configuration
	 */
public:
	void printInformation()	const
	{
		std::cout << std::endl;
		std::cout << "grid_res: " << grid_res[0] << ", " << grid_res[1] << std::endl;
		std::cout << "grid_data_size: " << grid_data_size[0] << ", " << grid_data_size[1] << std::endl;
		std::cout << "grid_number_elements: " << grid_number_elements << std::endl;

#if SWEET_USE_LIBFFT
		std::cout << std::endl;
		std::cout << "spectral_modes: " << spectral_modes[0] << ", " << spectral_modes[1] << std::endl;
		std::cout << "spectral_representable_modes: " << spectral_representable_modes[0] << ", " << spectral_representable_modes[1] << std::endl;
		std::cout << "spectral_data_size: " << spectral_data_size[0] << ", " << spectral_data_size[1] << std::endl;
		std::cout << "spectral_array_data_number_of_elements: " << spectral_array_data_number_of_elements << std::endl;
		std::cout << std::endl;
		std::cout << "spectral_data_iteration_ranges [0][0]: " << spectral_data_iteration_ranges[0][0][0] << ", " << spectral_data_iteration_ranges[0][0][1] << std::endl;
		std::cout << "spectral_data_iteration_ranges [0][1]: " << spectral_data_iteration_ranges[0][1][0] << ", " << spectral_data_iteration_ranges[0][1][1] << std::endl;
		std::cout << "spectral_data_iteration_ranges [1][0]: " << spectral_data_iteration_ranges[1][0][0] << ", " << spectral_data_iteration_ranges[1][0][1] << std::endl;
		std::cout << "spectral_data_iteration_ranges [1][1]: " << spectral_data_iteration_ranges[1][1][0] << ", " << spectral_data_iteration_ranges[1][1][1] << std::endl;
		std::cout << std::endl;
		std::cout << "spectral_complex_data_size: " << spectral_complex_data_size[0] << ", " << spectral_complex_data_size[1] << std::endl;
		std::cout << "spectral_complex_array_data_number_of_elements: " << spectral_complex_array_data_number_of_elements << std::endl;
		std::cout << std::endl;
		std::cout << "spectral_complex_ranges [0][0]: " << spectral_complex_ranges[0][0][0] << ", " << spectral_complex_ranges[0][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [0][1]: " << spectral_complex_ranges[0][1][0] << ", " << spectral_complex_ranges[0][1][1] << std::endl;
		std::cout << "spectral_complex_ranges [1][0]: " << spectral_complex_ranges[1][0][0] << ", " << spectral_complex_ranges[1][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [1][1]: " << spectral_complex_ranges[1][1][0] << ", " << spectral_complex_ranges[1][1][1] << std::endl;
		std::cout << "spectral_complex_ranges [2][0]: " << spectral_complex_ranges[2][0][0] << ", " << spectral_complex_ranges[2][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [2][1]: " << spectral_complex_ranges[2][1][0] << ", " << spectral_complex_ranges[2][1][1] << std::endl;
		std::cout << "spectral_complex_ranges [3][0]: " << spectral_complex_ranges[3][0][0] << ", " << spectral_complex_ranges[3][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [3][1]: " << spectral_complex_ranges[3][1][0] << ", " << spectral_complex_ranges[3][1][1] << std::endl;
		std::cout << std::endl;
#endif
	}


	/*!
	 * Load wisdom from predefined file 'sweet_fftw'
	 */
public:
	static
	bool loadWisdom(
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_planCacheMode	//!< Mode of plan cache
	)
	{
		static const char *wisdom_file = "sweet_fftw";

		int wisdom_plan_loaded = fftw_import_wisdom_from_filename(wisdom_file);
		if (wisdom_plan_loaded == 0)
		{
			std::cerr << "Failed to load FFTW wisdom from file '" << wisdom_file << "'" << std::endl;
			if (i_planCacheMode & sweet::Tools::TransformationPlans::REQUIRE_LOAD)
			{
				std::cerr << "IMPORTANT: Usage of FFTW wisdoms file enforced, hence exiting here" << std::endl;
				exit(1);
			}
		}

		return true;
	}


	/*!
	 * Store wisdom to the file 'sweet_fftw'
	 */
public:
	static
	bool storeWisdom()
	{
		static const char *wisdom_file = "sweet_fftw";

		int wisdom_plan_loaded = fftw_export_wisdom_to_filename(wisdom_file);
		if (wisdom_plan_loaded == 0)
		{
			std::cerr << "Failed to store FFTW wisdom to file " << wisdom_file << std::endl;
			exit(1);
		}

		return true;
	}


	/*!
	 * Reference counter of all generated FFTW plans.
	 * We use this to entirely also shutdown FFTW if there are no active plan anymore.
	 */
private:
	int& _refCounterFftwPlans()
	{
#if SWEET_THREADING_SPACE && SWEET_DEBUG
		if (omp_get_level() != 0)
			SWEETErrorFatal("Cart2DDataConfig is not threadsafe, but called inside parallel region with more than one thread!!!");
#endif

		static int ref_counter = 0;
		return ref_counter;
	}



	/*!
	 * Setup internal data of this config based on other initialized fields.
	 */
private:
	bool _setupInternalData(
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		reuse_spectral_transformation_plans = i_reuse_spectral_transformation_plans;

		if (fftw_initialized)
		{
			//
			// refCounter()--;
			// is inside cleanup!
			_clear_data();
		}
		else
		{
			// use REF counter in case of multiple plans
			// this allows a clean cleanup of fftw library
			fftw_initialized = true;
		}

		// FFTW PLANS are allocated below
		_refCounterFftwPlans()++;

		grid_data_size[0] = grid_res[0];
		grid_data_size[1] = grid_res[1];

		grid_number_elements = grid_res[0]*grid_res[1];

#if SWEET_USE_LIBFFT
		if (
			grid_res[0] < spectral_modes[0]	||
			grid_res[1] < spectral_modes[1]
		)
			return error.set("Lower physical resolution than spectral resolution not supported!");

		SWEET_ASSERT(grid_res[0] > 0);
		SWEET_ASSERT(grid_res[1] > 0);

		SWEET_ASSERT(spectral_modes[0] > 0);
		SWEET_ASSERT(spectral_modes[1] > 0);

#if SWEET_THREADING_SPACE

		// Is this the first instance?
		if (_refCounterFftwPlans() == 1)
		{
			// initialise FFTW with spatial parallelization
			int retval = fftw_init_threads();
			if (retval == 0)
				return error.set("fftw_init_threads() failed");

			int nthreads = omp_get_max_threads();

			fftw_plan_with_nthreads(nthreads);

		}
#endif

		if (_refCounterFftwPlans() == 1)
		{
			// load wisdom the first time
			// this must be done after initializing the threading!
			if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::LOAD)
				loadWisdom(reuse_spectral_transformation_plans);
		}


		unsigned int flags = 0;

		// allow destroying input for faster transformations
		//flags |= FFTW_DESTROY_INPUT;

		if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::QUICK)
		{
			flags |= FFTW_ESTIMATE;
		}
		else
		{
			// estimation base don workload
			unsigned int cells = grid_data_size[0]*grid_data_size[1];

			if (cells < 32*32)
				//flags |= FFTW_EXHAUSTIVE;
				flags |= FFTW_MEASURE;
			else if (cells < 256*256)
				flags |= FFTW_MEASURE;
			else
				flags |= FFTW_MEASURE;
				//flags |= FFTW_PATIENT;

			if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::REQUIRE_LOAD)
			{
				flags |= FFTW_WISDOM_ONLY;	// only allow plans from wisdom
			}
		}


		/*
		 * REAL PHYSICAL SPACE DATA (REAL to COMPLEX FFT)
		 */
		{
			// real-to-complex storage representation
			spectral_data_size[0] = grid_data_size[0]/2+1;
			spectral_data_size[1] = grid_data_size[1];

			if ((grid_data_size[0] & 1) == 1)
				return error.set("Unsupported odd resolution in x-direction");

			if ((grid_data_size[1] & 1) == 1)
				return error.set("Unsupported odd resolution in y-direction");

#if SWEET_USE_CART2D_SPECTRAL_DEALIASING

			/*
			 * For more information, have a look at
			 * doc/software_development_discussions/antialiasing/implementation_strategy.pdf
			 */

			if (	spectral_modes[0] == grid_res[0] ||
					spectral_modes[1] == grid_res[1]
			)
				return error.set("Aliasing doesn't make sense since physical resolution is identical to spectral");

			spectral_data_iteration_ranges[0][0][0] = 0;
			spectral_data_iteration_ranges[0][0][1] = (grid_data_size[0]-1)/3;
			spectral_data_iteration_ranges[0][1][0] = 0;
			spectral_data_iteration_ranges[0][1][1] = (grid_data_size[1]-1)/3;


#else

			/*
			 * The central mode is the Shannon-Nyquist one
			 * => Remove this one, since this is not tracked in transformations
			 */
			spectral_data_iteration_ranges[0][0][0] = 0;
			spectral_data_iteration_ranges[0][0][1] = spectral_data_size[0]-1;		// Shannon-Nyquist
			spectral_data_iteration_ranges[0][1][0] = 0;
			spectral_data_iteration_ranges[0][1][1] = spectral_data_size[1]/2;

#endif


			spectral_data_iteration_ranges[1][0][0] = spectral_data_iteration_ranges[0][0][0];
			spectral_data_iteration_ranges[1][0][1] = spectral_data_iteration_ranges[0][0][1];
			spectral_data_iteration_ranges[1][1][0] = spectral_data_size[1] - spectral_data_iteration_ranges[0][1][1] + 1;
			spectral_data_iteration_ranges[1][1][1] = spectral_data_size[1];

			spectral_representable_modes[0] = spectral_data_iteration_ranges[0][0][1];
			spectral_representable_modes[1] = spectral_data_iteration_ranges[0][1][1];

			spectral_array_data_number_of_elements = spectral_data_size[0]*spectral_data_size[1];



			/*
			 * Grid space data
			 */
			double *data_physical = sweet::Memory::MemBlockAlloc::alloc<double>(grid_number_elements*sizeof(double));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < grid_number_elements; i++)
				data_physical[i] = 1;	// dummy data

			/*
			 * Spectral space data
			 */
			std::complex<double> *data_spectral = sweet::Memory::MemBlockAlloc::alloc< std::complex<double> >(spectral_array_data_number_of_elements*sizeof(std::complex<double>));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < spectral_array_data_number_of_elements; i++)
				data_spectral[i] = 1;	// dummy data

			_fftw_plan_forward =
				fftw_plan_dft_r2c_2d(
					grid_data_size[1],	// n0 = ny
					grid_data_size[0],	// n1 = nx
					data_physical,
					(fftw_complex*)data_spectral,
					flags | FFTW_PRESERVE_INPUT
				);

			if (_fftw_plan_forward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to get forward plan dft_r2c fftw" << std::endl;
				std::cerr << "r2c preverse_input forward " << grid_res[0] << " x " << grid_res[1] << std::endl;
				std::cerr << "FFTW-wisdom plan: rf" << grid_res[0] << "x" << grid_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::REQUIRE_LOAD)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			_fftw_plan_backward =
					fftw_plan_dft_c2r_2d(
						grid_res[1],	// n0 = ny
						grid_res[0],	// n1 = nx
						(fftw_complex*)data_spectral,
						data_physical,
						flags
					);

			if (_fftw_plan_backward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to get backward plan dft_c2r fftw" << std::endl;
				std::cerr << "r2c backward " << grid_res[0] << " x " << grid_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: rb" << grid_res[0] << "x" << grid_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::REQUIRE_LOAD)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			sweet::Memory::MemBlockAlloc::free(data_physical, grid_number_elements*sizeof(double));
			sweet::Memory::MemBlockAlloc::free(data_spectral, spectral_array_data_number_of_elements*sizeof(std::complex<double>));

			// Backward scaling factor
			fftw_backward_scale_factor = 1.0/((double)(grid_data_size[0]*grid_data_size[1]));
		}


		/*
		 * COMPLEX PHYSICAL SPACE DATA
		 */
		{
			spectral_complex_data_size[0] = grid_data_size[0];
			spectral_complex_data_size[1] = grid_data_size[1];

			if ((spectral_complex_data_size[0] & 1) == 1)
				SWEETErrorFatal("Not supported c");

			if ((spectral_complex_data_size[1] & 1) == 1)
				SWEETErrorFatal("Not supported d");

			spectral_complex_ranges[0][0][0] = spectral_data_iteration_ranges[0][0][0];
			spectral_complex_ranges[0][0][1] = spectral_data_iteration_ranges[0][0][1];
			spectral_complex_ranges[0][1][0] = spectral_data_iteration_ranges[0][1][0];
			spectral_complex_ranges[0][1][1] = spectral_data_iteration_ranges[0][1][1];

			spectral_complex_ranges[1][0][0] = spectral_complex_ranges[0][0][0];
			spectral_complex_ranges[1][0][1] = spectral_complex_ranges[0][0][1];
			spectral_complex_ranges[1][1][0] = spectral_data_size[1] - spectral_complex_ranges[0][1][1] + 1;
			spectral_complex_ranges[1][1][1] = spectral_data_size[1];

			spectral_complex_ranges[2][0][0] = spectral_complex_data_size[0] - spectral_complex_ranges[0][0][1] + 1;
			spectral_complex_ranges[2][0][1] = spectral_complex_data_size[0];
			spectral_complex_ranges[2][1][0] = 0;
			spectral_complex_ranges[2][1][1] = spectral_complex_ranges[0][1][1];

			spectral_complex_ranges[3][0][0] = spectral_complex_data_size[0] - spectral_complex_ranges[0][0][1] + 1;
			spectral_complex_ranges[3][0][1] = spectral_complex_data_size[0];
			spectral_complex_ranges[3][1][0] = spectral_complex_data_size[1] - spectral_complex_ranges[0][1][1] + 1;
			spectral_complex_ranges[3][1][1] = spectral_complex_data_size[1];

			spectral_complex_array_data_number_of_elements = spectral_complex_data_size[0]*spectral_complex_data_size[1];

			/*
			 * Grid space data
			 */
			std::complex<double> *data_physical = sweet::Memory::MemBlockAlloc::alloc< std::complex<double> >(grid_number_elements*sizeof(std::complex<double>));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < grid_number_elements; i++)
				data_physical[i] = 1;	// dummy data

			/*
			 * Spectral space data
			 */
			std::complex<double> *data_spectral = sweet::Memory::MemBlockAlloc::alloc< std::complex<double> >(spectral_complex_array_data_number_of_elements*sizeof(std::complex<double>));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < spectral_complex_array_data_number_of_elements; i++)
				data_spectral[i] = 1;	// dummy data


			_fftw_plan_complex_forward =
					fftw_plan_dft_2d(
						grid_res[1],
						grid_res[0],
						(fftw_complex*)data_physical,
						(fftw_complex*)data_spectral,
						FFTW_FORWARD,
						flags
					);

			if (_fftw_plan_complex_forward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to create complex forward plan for fftw" << std::endl;
				std::cerr << "complex forward preverse_input forward " << grid_res[0] << " x " << grid_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << grid_res[0] << "x" << grid_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans == 2)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			_fftw_plan_complex_backward =
					fftw_plan_dft_2d(
						grid_res[1],
						grid_res[0],
						(fftw_complex*)data_spectral,
						(fftw_complex*)data_physical,
						FFTW_BACKWARD,
						flags
					);

			if (_fftw_plan_complex_backward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to create complex backward plan for fftw" << std::endl;
				std::cerr << "complex backward preverse_input forward " << grid_res[0] << " x " << grid_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << grid_res[0] << "x" << grid_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans == 2)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			sweet::Memory::MemBlockAlloc::free(data_physical, grid_number_elements*sizeof(std::complex<double>));
			sweet::Memory::MemBlockAlloc::free(data_spectral, spectral_complex_array_data_number_of_elements*sizeof(std::complex<double>));
		}
#endif
		return true;
	}



#if SWEET_USE_CART2D_SPECTRAL_SPACE

	/*
	 * Return the spectral range
	 */
public:
	std::size_t getSpectralIterationRangeArea(int i)	const
	{
		SWEET_ASSERT(i >= 0);
		SWEET_ASSERT(i <= 2);
		return	(spectral_data_iteration_ranges[i][0][1] - spectral_data_iteration_ranges[i][0][0])*
				(spectral_data_iteration_ranges[i][1][1] - spectral_data_iteration_ranges[i][1][0]);
	}

#endif



#if SWEET_USE_LIBFFT
	/*!
	 * For real-valued physical data, convert from physical to spectral space
	 */
public:
	void fft_grid_2_spectral_OUTOFPLACE(
			double *i_grid_data,
			std::complex<double> *o_spectral_data
	)	const
	{
		fftw_execute_dft_r2c(
				_fftw_plan_forward,
				i_grid_data,
				(fftw_complex*)o_spectral_data
			);
	}



	/*!
	 * For real-valued physical data, convert spectral to physical space
	 *
	 * WARNING: This is an in-place transformation
	 * where the spectral data will be overwritten.
	 */
	void fft_spectral_2_grid_INPLACE(
			std::complex<double> *i_spectral_data,
			double *o_grid_data
	)	const
	{
		SWEET_ASSERT(_fftw_plan_backward != nullptr);
		SWEET_ASSERT(i_spectral_data != nullptr);
		SWEET_ASSERT(o_grid_data != nullptr);

		fftw_execute_dft_c2r(
				_fftw_plan_backward,
				(fftw_complex*)i_spectral_data,
				o_grid_data
			);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < grid_number_elements; i++)
			o_grid_data[i] *= fftw_backward_scale_factor;
	}



	/*!
	 * For complex-valued physical data, convert from physical to spectral space
	 */
	void fft_complex_grid_2_spectral_OUTOFPLACE(
			std::complex<double> *i_grid_data,
			std::complex<double> *o_spectral_data
	)	const
	{
		fftw_execute_dft(
				_fftw_plan_complex_forward,
				(fftw_complex*)i_grid_data,
				(fftw_complex*)o_spectral_data
			);
	}




	/*!
	 * For complex-valued physical data, convert from spectral to physical space
	 */
	void fft_complex_spectral_2_grid_OUTOFPLACE(
			std::complex<double> *i_spectral_data,
			std::complex<double> *o_grid_data
	)	const
	{
		fftw_execute_dft(
				_fftw_plan_complex_backward,
				(fftw_complex*)i_spectral_data,
				(fftw_complex*)o_grid_data
			);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < grid_number_elements; i++)
			o_grid_data[i] *= fftw_backward_scale_factor;
	}
#endif

	/*
	 * Setup Cart2DDataConfig with concretely given data
	 */
public:
	bool setup(
			int i_grid_res_x,
			int i_grid_res_y,

			int i_spectral_modes_x,
			int i_spectral_modes_y,

			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		// we always clear the data - just in case
		clear();

		grid_res[0] = i_grid_res_x;
		grid_res[1] = i_grid_res_y;

#if SWEET_USE_LIBFFT
		spectral_modes[0] = i_spectral_modes_x;
		spectral_modes[1] = i_spectral_modes_y;
#endif

		return _setupInternalData(i_reuse_spectral_transformation_plans);
	}



	/*
	 * Just a convenience function
	 */
public:
	bool setupAuto(
			int io_grid_res[2],
			int io_spectral_modes[2],
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		if (io_grid_res[0] > 0 && io_spectral_modes[0] > 0)
		{
			return setup(	io_grid_res[0],
					io_grid_res[1],
					io_spectral_modes[0],
					io_spectral_modes[1],
					i_reuse_spectral_transformation_plans
				);
		}

		if (io_grid_res[0] > 0)
		{
			if (!setupAutoSpectralSpaceFromGrid(
					io_grid_res[0],
					io_grid_res[1],
					i_reuse_spectral_transformation_plans
				))
				return false;

#if SWEET_USE_LIBFFT
			io_spectral_modes[0] = spectral_modes[0];
			io_spectral_modes[1] = spectral_modes[1];
#endif
			return true;
		}

		if (io_spectral_modes[0] > 0)
		{
#if SWEET_USE_LIBFFT
			setupAutoGridSpaceFromSpectral(
					io_spectral_modes[0],
					io_spectral_modes[1],
					i_reuse_spectral_transformation_plans
				);

			io_grid_res[0] = grid_res[0];
			io_grid_res[1] = grid_res[1];
			return true;
#else
			return Error.set("Setup with spectral modes not enabled");
#endif
		}

		return error.set("No resolution/modes selected");
	}


	bool setupAuto(Shack &i_shackCart2DDataOps)
	{
		return setupAuto(
				i_shackCart2DDataOps.space_res_physical,
				i_shackCart2DDataOps.space_res_spectral,
				i_shackCart2DDataOps.reuse_spectral_transformation_plans
			);
	}
	bool setupAuto(Shack *i_shackCart2DDataOps)
	{
		return setupAuto(
				i_shackCart2DDataOps->space_res_physical,
				i_shackCart2DDataOps->space_res_spectral,
				i_shackCart2DDataOps->reuse_spectral_transformation_plans
			);
	}

	/*
	 * Just a convenience function
	 */
public:
	bool setupAutoSpectralSpaceFromGrid(
			int i_grid_res[2],
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		return setupAutoSpectralSpaceFromGrid(
				i_grid_res[0],
				i_grid_res[1],
				i_reuse_spectral_transformation_plans
		);
	}



public:
	bool setupAutoSpectralSpaceFromGrid(
			int i_grid_res_x,
			int i_grid_res_y,
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		return setup(
				i_grid_res_x,
				i_grid_res_y,

#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
				// REDUCTION IN EFFECTIVE SPECTRAL MODE RESOLUTION TO CUT OFF ANTI-ALIASED MODES
				(i_grid_res_x*2)/3,
				(i_grid_res_y*2)/3,
#else
	#if SWEET_USE_LIBFFT
				i_grid_res_x,
				i_grid_res_y,
	#else
				0,
				0,
	#endif
#endif
				i_reuse_spectral_transformation_plans
		);
	}


	bool setupAutoSpectralSpaceFromGrid(
			int i_grid_res_x,
			int i_grid_res_y,
			int *o_spectral_res_x,
			int *o_spectral_res_y,
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		if (!setupAutoSpectralSpaceFromGrid(
				i_grid_res_x,
				i_grid_res_y,
				i_reuse_spectral_transformation_plans
			))
			return false;

#if SWEET_USE_LIBFFT
		*o_spectral_res_x = spectral_modes[0];
		*o_spectral_res_y = spectral_modes[1];
#else
		*o_spectral_res_x = 0;
		*o_spectral_res_y = 0;
#endif
		return true;
	}



#if SWEET_USE_LIBFFT

public:
	bool setupAutoGridSpaceFromSpectral(
			int i_spectral_res_x,
			int i_spectral_res_y,
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		return setup(
#if SWEET_USE_CART2D_SPECTRAL_DEALIASING

				// REDUCTION IN EFFECTIVE SPECTRAL MODE RESOLUTION TO CUT OFF ANTI-ALIASED MODES
				(i_spectral_res_x*3+1)/2,
				(i_spectral_res_y*3+1)/2,
#else
				i_spectral_res_x,
				i_spectral_res_y,
#endif

				i_spectral_res_x,
				i_spectral_res_y,

				i_reuse_spectral_transformation_plans
		);
	}


public:
	bool setupAutoGridSpaceFromSpectral(
			int i_spectral_res_x,
			int i_spectral_res_y,
			int *o_grid_res_x,
			int *o_grid_res_y,
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		if (!setupAutoGridSpaceFromSpectral(
				i_spectral_res_x,
				i_spectral_res_y,
				i_reuse_spectral_transformation_plans
		))
			return false;

		*o_grid_res_x = grid_res[0];
		*o_grid_res_y = grid_res[1];

		return true;
	}

	bool setupAdditionalModes(
			Config *i_cart2dConfig,
			int i_additional_modes_x,
			int i_additional_modes_y,
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		return setupAutoGridSpaceFromSpectral(
				i_cart2dConfig->spectral_modes[0] + i_additional_modes_x,
				i_cart2dConfig->spectral_modes[1] + i_additional_modes_y,
				i_reuse_spectral_transformation_plans
		);
	}



	// TODO: CHECK THIS
	inline
	std::size_t getArrayIndexByModes(
			int n,
			int m
	)	const
	{
		SWEET_ASSERT(n >= 0);

		return n * spectral_data_size[0] + m;
	}


	inline
	std::size_t getArrayIndexByModes_Complex(
			int n,
			int m
	)	const
	{
		SWEET_ASSERT(n >= 0);

		int idx =  n * spectral_complex_data_size[0] + m;
		return idx;
	}

#endif



	void _clear_data()
	{
	}


	void clear()
	{
		if (!fftw_initialized)
			return;

#if SWEET_USE_LIBFFT
		fftw_destroy_plan(_fftw_plan_forward);
		fftw_destroy_plan(_fftw_plan_backward);

		fftw_destroy_plan(_fftw_plan_complex_forward);
		fftw_destroy_plan(_fftw_plan_complex_backward);

		_refCounterFftwPlans()--;
		SWEET_ASSERT(_refCounterFftwPlans() >= 0);

		if (_refCounterFftwPlans() == 0)
		{
			// backup wisdom
			if (reuse_spectral_transformation_plans == 1)
				storeWisdom();

#if SWEET_THREADING_SPACE
			fftw_cleanup_threads();
#endif
			fftw_cleanup();
		}
#endif

		grid_res[0] = 0;
		grid_res[1] = 0;

#if SWEET_USE_LIBFFT
		spectral_modes[0] = 0;
		spectral_modes[1] = 0;
#endif

		fftw_initialized = false;
	}



	~Config()
	{
		clear();
	}
};

}}}

#endif
