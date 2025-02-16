/*
 * SPHSetup.hpp
 *
 *  Created on: 12 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_CONFIG_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_CONFIG_HPP

#include <sweet/Data/Sphere2D/Shack.hpp>

#include <fftw3.h>
#include <iostream>
#include <cmath>
#include <sweet/IO/FileOperations.hpp>
#include <sweet/Tools/TransformationPlans.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>
#include <sweet/LibMath/shtns_inc.hpp>
#include <stdexcept>

#if SWEET_MPI
#	include <mpi.h>
#endif


#define SPHERE2D_DATA_LON_CONTIGUOUS	SHT_PHI_CONTIGUOUS
#define SPHERE2D_DATA_LAT_CONTIGUOUS	SHT_THETA_CONTIGUOUS

// SWEET
//#define SPHERE2D_DATA_GRID_LAYOUT	SPHERE2D_DATA_LAT_CONTIGUOUS

// SHTNS shallow water example
#define SPHERE2D_DATA_GRID_LAYOUT		SPHERE2D_DATA_LON_CONTIGUOUS

namespace sweet {
namespace Data {
namespace Sphere2D {


/*!
 * \brief Configuration class for Sphere2D including transformations
 */
class Config
{
	friend class Operators;
	friend class Sphere2DData;

public:
	Error::Base error;
	shtns_cfg shtns;

	int ndim = 2;

	/*
	 * Verbosity of outputs
	 */
	int verbosity;

	/*
	 * Number of threads
	 */
	int numThreads;

	/*!
	 * Number of longitudes
	 */
public:
	int grid_num_lon;

	/*!
	 * Number of latitudes
	 */
public:
	int grid_num_lat;

	/*!
	 * Number of total longitudes and latitudes
	 */
public:
	std::size_t grid_number_elements;


	/*!
	 * Number of modes
	 */
public:
	int spectral_modes_n_max;
	int spectral_modes_m_max;

	double shtns_error = 0;
//	double shtns_error = 1.e-10;
//	double shtns_error = 1.e-6;

	/*!
	 * Total number of modes (complex valued)
	 */
	int spectral_array_data_number_of_elements;

	/*!
	 * Number of elements for SPH which is based on
	 * the *** complex-valued physical data ***
	 */
	std::size_t spectral_complex_array_data_number_of_elements;


	/*!
	 * Array with latitude phi angle values
	 *
	 * WARNING: Phi is not the phi from SHTNS!
	 */
public:
	double *lat;


	/*!
	 * Array with mu = sin(phi) values
	 */
public:
	double *lat_gaussian;


	/*!
	 * Array with comu = cos(phi) values
	 */
public:
	double *lat_cogaussian;



public:
	Config()	:
		shtns(nullptr),
		grid_num_lon(-1),
		grid_num_lat(-1),
		grid_number_elements(-1),

		spectral_modes_n_max(-1),
		spectral_modes_m_max(-1),
		spectral_array_data_number_of_elements(-1),
		spectral_complex_array_data_number_of_elements(-1),

		lat(nullptr),
		lat_gaussian(nullptr),
		lat_cogaussian(nullptr)
	{
	}

	const
	std::string getUniqueIDString()	const
	{
		return getConfigInformationString();
	}


	const
	std::string getConfigInformationString()	const
	{
		std::ostringstream buf;
		buf << "M" << spectral_modes_m_max << "," << spectral_modes_n_max << "_N" << grid_num_lon << "," << grid_num_lat;
		buf << " total_spec_modes: " << spectral_array_data_number_of_elements;
		return buf.str();
	}


	inline
	std::size_t getArrayIndexByModes(
			int n,
			int m
	)	const
	{
		SWEET_ASSERT(n >= 0);
		SWEET_ASSERT(n >= m);

		return (m*(2*spectral_modes_n_max-m+1)>>1)+n;
	}



	inline
	std::size_t getArrayIndexByModes_Complex(
			int n,
			int m
	)	const
	{
		SWEET_ASSERT(n >= 0);
		SWEET_ASSERT(n >= std::abs(m));

		int idx = n*n+(m+n);
		return idx;
	}


	/*!
	 * Return indices with N-variables compactly stored
	 */
	inline
	std::size_t getArrayIndexByModes_Complex_NCompact(
			int n,
			int m
	)	const
	{
		SWEET_ASSERT(n >= 0);
		SWEET_ASSERT(n >= std::abs(m));

		if (m < 0)
		{
			int minc = spectral_modes_m_max+m;
			int row_idx = (((minc-1)*minc)>>1) + minc;
			int idx =  row_idx + (n+m);
			return idx;
		}
		else
		{
			int rel_idx = (m*(2*spectral_modes_n_max-m+1)>>1)+n;
			return ((spectral_modes_n_max*(spectral_modes_n_max+1))>>1)+rel_idx;
		}
	}


private:
	bool setup_data()
	{
		if (verbosity > 0)
			shtns_print_cfg(shtns);

		grid_num_lat = shtns->nlat;
		grid_num_lon = shtns->nphi;
		grid_number_elements = shtns->nspat;

		spectral_modes_n_max = shtns->lmax;
		spectral_modes_m_max = shtns->mmax;
		spectral_array_data_number_of_elements = shtns->nlm;
		spectral_complex_array_data_number_of_elements = (spectral_modes_n_max+1)*(spectral_modes_m_max+1);

		if (spectral_modes_n_max != spectral_modes_m_max)
		{
			std::cerr << "only spec_n_max == spec_m_max currently supported!" << std::endl;
			SWEET_ASSERT(false);
			exit(1);
		}

#if SWEET_DEBUG
		/*
		 * Some safety checks to make sure that we really get what we've asked for
		 */

		/*
		 * TEST: iteration over the modes n,m for real-valued physical space
		 */
		{
			int idx = 0;
			for (int m = 0; m <= spectral_modes_m_max; m++)
			{
				int test_idx = getArrayIndexByModes(m,m);

				if (test_idx != idx)
				{
					std::cerr << "IDX TEST NOT SUCCESSFUL (real-valued physical transformation)" << std::endl;
					std::cout << "n=" << m << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
					exit(1);
				}

				for (int n = m; n <= spectral_modes_n_max; n++)
				{

					int test_idx2 = getArrayIndexByModes(n,m);
					if (test_idx2 != idx)
					{
						std::cerr << "IDX TEST2 NOT SUCCESSFUL (real-valued physical transformation)" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}
					idx++;
				}
			}

			if (idx != spectral_array_data_number_of_elements)
			{
				std::cerr << "INTERNAL SPH ERROR (real-valued physical transformation)" << std::endl;
				SWEET_ASSERT(false);
				exit(1);
			}
		}



		/*
		 * TEST: iteration over the modes n,m for complex-valued physical space
		 *
		 * Note: In SHTNS, the m-coefficients are compactly stored for individual n'l
		 */
		{
			std::size_t idx = 0;
			for (int n = 0; n <= spectral_modes_n_max; n++)
			{
				std::size_t test_idx = getArrayIndexByModes_Complex(n,-n);

				if (test_idx != idx)
				{
					std::cerr << "IDX TEST NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
					std::cout << "n=" << n << ", m=" << -n << "     " << idx << ", " << test_idx << std::endl;
					exit(1);
				}

				for (int m = -n; m <= n; m++)
				{
					std::size_t test_idx2 = getArrayIndexByModes_Complex(n,m);

					if (test_idx2 != idx)
					{
						std::cerr << "IDX TEST2 NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}

					idx++;
				}
			}


			if (idx != spectral_complex_array_data_number_of_elements)
			{
				std::cerr << "INTERNAL SPH ERROR" << std::endl;
				SWEET_ASSERT(false);
				exit(1);
			}
		}


		/*
		 * TEST: iteration over the modes n,m for complex-valued physical space
		 *
		 * This version tests for the n-coefficients compactly stored for individual m's
		 */
		{
			int idx = 0;
			for (int m = -spectral_modes_m_max; m <= spectral_modes_m_max; m++)
			{
				int test_idx = getArrayIndexByModes_Complex_NCompact(std::abs(m),m);

				if (test_idx != idx)
				{
					std::cerr << "IDX TEST NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
					std::cout << "n=" << m << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
					exit(1);
				}


				for (int n = std::abs(m); n <= spectral_modes_n_max; n++)
				{
					int test_idx2 = getArrayIndexByModes_Complex_NCompact(n,m);

					if (test_idx2 != idx)
					{
						std::cerr << "IDX TEST2 NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}

					idx++;
				}
			}
		}
#endif

		lat = (double*)fftw_malloc(sizeof(double)*grid_num_lat);

		/*
		 * Colatitude is 0 at the north pole and 180 at the south pole
		 *
		 * WARNING: The latitude degrees are not spaced equidistantly in the angles!!!!
		 * Those points are computed for optimal Gauss quadrature.
		 * They are close to equidistant spacing, but not fully equidistant.
		 *
		 * We have to use the shtns->ct lookup table
		 */
		for (int i = 0; i < grid_num_lat; i++)
			lat[i] = M_PI_2 - std::acos(shtns->ct[i]);

		lat_gaussian = (double*)fftw_malloc(sizeof(double)*shtns->nlat);
		for (int i = 0; i < grid_num_lat; i++)
			lat_gaussian[i] = shtns->ct[i];		//! sin(phi) (SHTNS stores cos(phi))

		lat_cogaussian = (double*)fftw_malloc(sizeof(double)*shtns->nlat);
		for (int i = 0; i < grid_num_lat; i++)
			lat_cogaussian[i] = shtns->st[i];	//! cos(phi) (SHTNS stores sin(phi))


		return true;
	}


	int getFlags(
			int i_reuse_spectral_transformation_plans
	)
	{
		std::string cache_filename = "shtns_cfg";
		std::string cache_filename_fftw = "shtns_cfg_fftw";

		int flags = 0;

		// lat or lon continue contiguously stored in grid space
		flags |= SPHERE2D_DATA_GRID_LAYOUT;

		if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::QUICK)
		{
			/*
			 * No special initialization, no caching
			 */
			flags |= sht_quick_init;
		}

		if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::LOAD)
		{
			if (verbosity > 0)
				std::cout << " + Trying to load plans" << std::endl;

			bool plans_exist = true;
			if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::REQUIRE_LOAD)
			{
				if (!sweet::IO::FileOperations::file_exists(cache_filename))
				{
					if ((i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::REQUIRE_LOAD) == 0)
						throw std::runtime_error(std::string("File '"+cache_filename+"' does not exist"));

					plans_exist = false;
				}

				if (!sweet::IO::FileOperations::file_exists(cache_filename_fftw))
				{
					if ((i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::REQUIRE_LOAD) == 0)
						throw std::runtime_error(std::string("File '"+cache_filename_fftw+"' does not exist"));

					plans_exist = false;
				}
			}

			if (!plans_exist)
			{
				if (verbosity > 0)
					std::cout << " + WARNING: No existing plan found" << std::endl;
			}
			else
			{
				flags |= SHT_LOAD_SAVE_CFG;
			}
		}
		else if (i_reuse_spectral_transformation_plans & sweet::Tools::TransformationPlans::SAVE)
		{
			if (verbosity > 0)
				std::cout << " + Generating and storing SH transformation plans" << std::endl;

			flags |= SHT_LOAD_SAVE_CFG;
		}


		return flags;
	}


	void setupSHTNSBoilerplateCode()
	{
		shtns_verbose(verbosity);

#if SWEET_THREADING_SPACE
		// enable multi-threaded transforms (if supported).
		if (numThreads <= 0)
		{
			shtns_use_threads(0);
		}
		else
		{
			shtns_use_threads(numThreads);
		}
#else
		shtns_use_threads(1);	// value of 1 disables threading
#endif
	}


public:
	bool setup(
		int nphi,	// physical
		int nlat,

		int mmax,	// spectral
		int nmax,

		sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_transformation_plans,
		int i_verbosity,
		int i_numThreads
	)
	{
		verbosity = i_verbosity;
		numThreads = i_numThreads;

		mmax--;
		nmax--;

		setupSHTNSBoilerplateCode();

		shtns = shtns_create(
				nmax,
				mmax,
				1,
				(shtns_norm)((int)sht_orthonormal + SHT_NO_CS_PHASE)
			);

#if SWEET_MPI
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		if (mpi_rank == 0 && i_reuse_transformation_plans)
			MPI_Barrier(MPI_COMM_WORLD);
#endif

		shtns_set_grid(
				shtns,
				(shtns_type)getFlags(i_reuse_transformation_plans),
				shtns_error,
				nlat,		// number of latitude grid points
				nphi		// number of longitude grid points
			);

#if SWEET_MPI
		if (mpi_rank > 0 && i_reuse_transformation_plans)
			MPI_Barrier(MPI_COMM_WORLD);
#endif

		return setup_data();
	}


	/*!
	 * Setup with given modes.
	 * Spatial resolution will be determined automatically
	 */
	bool setupAutoGridSpace(
			int i_mmax,		//! longitude modes
			int i_nmax,		//! latitude modes
			int *o_nphi,	//! physical resolution along longitude
			int *o_nlat,	//! physical resolution along latitude
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_transformation_plans,
			int i_verbosity,
			int i_numThreads
	)
	{
		verbosity = i_verbosity;
		numThreads = i_numThreads;

		i_mmax--;
		i_nmax--;

		setupSHTNSBoilerplateCode();

		shtns = shtns_create(
				i_nmax,
				i_mmax,
				1,
				(shtns_norm)((int)sht_orthonormal + SHT_NO_CS_PHASE)
			);

		*o_nphi = 0;
		*o_nlat = 0;

#if SWEET_MPI
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		if (mpi_rank == 0 && i_reuse_transformation_plans)
		{
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif

		shtns_set_grid_auto(
				shtns,
				(shtns_type)getFlags(i_reuse_transformation_plans),
				shtns_error,
				2,		// use order 2
				o_nlat,
				o_nphi
			);

#if SWEET_MPI
		if (mpi_rank > 0 && i_reuse_transformation_plans)
		{
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif

		return setup_data();
	}


	/*!
	 * Setup with given modes.
	 * Spatial resolution will be determined automatically
	 */
	bool setupAutoGridSpace(
			int i_mmax,		//!< longitude modes
			int i_nmax,		//!< latitude modes
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_transformation_plans,
			int i_verbosity,
			int i_numThreads
	)
	{
		verbosity = i_verbosity;
		numThreads = i_numThreads;

		i_mmax--;
		i_nmax--;

		setupSHTNSBoilerplateCode();

		shtns = shtns_create(
				i_nmax,
				i_mmax,
				1,
				(shtns_norm)((int)sht_orthonormal + SHT_NO_CS_PHASE)
			);

		grid_num_lat = 0;
		grid_num_lon = 0;

#if SWEET_MPI
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		if (mpi_rank == 0 && i_reuse_transformation_plans)
		{
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif

		shtns_set_grid_auto(
				shtns,
				(shtns_type)getFlags(i_reuse_transformation_plans),
				shtns_error,
				2,		// use order 2
				&grid_num_lat,
				&grid_num_lon
			);

#if SWEET_MPI
		if (mpi_rank > 0 && i_reuse_transformation_plans)
		{
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif

		return setup_data();
	}


public:
	bool setupAuto(
			int io_grid_res[2],
			int io_spectral_modes[2],
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE &i_reuse_transformation_plans,
			int i_verbosity,
			int i_numThreads
	)
	{
		if (io_grid_res[0] > 0 && io_spectral_modes[0] > 0)
		{
			return setup(	io_grid_res[0],
					io_grid_res[1],
					io_spectral_modes[0],
					io_spectral_modes[1],
					i_reuse_transformation_plans,
					i_verbosity,
					i_numThreads
				);
		}

		if (io_spectral_modes[0] > 0)
		{
			bool retval = setupAutoGridSpace(
					io_spectral_modes[0],
					io_spectral_modes[1],
					i_reuse_transformation_plans,
					i_verbosity,
					i_numThreads
				);

			io_grid_res[0] = grid_num_lon;
			io_grid_res[1] = grid_num_lat;

			return retval;
		}

		error.set("No resolution/modes selected");
		return false;
	}


	bool setupAdditionalModes(
			const Sphere2D::Config *i_sphere2DDataConfig,
			int i_additional_modes_longitude,
			int i_additional_modes_latitude,
			Shack *shackSphere2DDataOps
	)
	{
		return setupAdditionalModes(
				i_sphere2DDataConfig,
				i_additional_modes_longitude,
				i_additional_modes_latitude,
				shackSphere2DDataOps->reuse_spectral_transformation_plans,
				shackSphere2DDataOps->sh_setup_verbosity,
				shackSphere2DDataOps->sh_setup_num_threads
			);
	}

public:
	bool setupAdditionalModes(
			const Sphere2D::Config *i_sphere2DDataConfig,
			int i_additional_modes_longitude,
			int i_additional_modes_latitude,
			sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE i_plan_load_save,
			int i_verbosity,
			int i_numThreads
	)
	{
		SWEET_ASSERT(shtns == nullptr);

		return setupAutoGridSpace(
				i_sphere2DDataConfig->spectral_modes_m_max + i_additional_modes_longitude,
				i_sphere2DDataConfig->spectral_modes_n_max + i_additional_modes_latitude,
				&grid_num_lon,
				&grid_num_lat,
				i_plan_load_save,
				i_verbosity,
				i_numThreads
		);
	}


public:
	bool setupAuto(Sphere2D::Shack *i_shackSphere2DDataOps)
	{
		return setupAuto(
				i_shackSphere2DDataOps->space_res_physical,
				i_shackSphere2DDataOps->space_res_spectral,
				i_shackSphere2DDataOps->reuse_spectral_transformation_plans,
				i_shackSphere2DDataOps->sh_setup_verbosity,
				i_shackSphere2DDataOps->sh_setup_num_threads
			);
	}


public:
	bool setupAuto(Sphere2D::Shack &i_shackSphere2DDataOps)
	{
		return setupAuto(&i_shackSphere2DDataOps);
	}



	void clear(
		bool i_full_reset
	)
	{
		if (shtns == nullptr)
			return;

		fftw_free(lat_cogaussian);
		lat_cogaussian = nullptr;

		fftw_free(lat_gaussian);
		lat_gaussian = nullptr;

		fftw_free(lat);
		lat = nullptr;

		shtns_unset_grid(shtns);
		shtns_destroy(shtns);
		shtns = nullptr;

		if (i_full_reset)
		{
#if SWEET_USE_THREADING
			fftw_cleanup_threads();
#endif
			fftw_cleanup();
		}
	}

	void clear()
	{
		clear(false);
	}


	~Config()
	{
		if (shtns != nullptr)
			clear();
	}
};

}}}

#endif
