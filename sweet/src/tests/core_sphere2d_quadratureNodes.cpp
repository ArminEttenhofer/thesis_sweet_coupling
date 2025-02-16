/*
 * test_sph_quadrature_nodes.hpp
 *
 *  Created on: 3 Feb 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_TESTSPH_QUADRATURE_NODES_HPP_
#define SRC_TESTSPH_QUADRATURE_NODES_HPP_

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/LibMath/BandedMatrixGridReal.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>

#include "../programs/PDE_SWESphere2D/TimeHelpers/SphBandedMatrix_GridReal.hpp"





class Sphere2DDataErrorCheck
{
public:
	static
	bool check(
			const sweet::Data::Sphere2D::DataSpectral &i_lhs,
			const sweet::Data::Sphere2D::DataSpectral &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		const sweet::Data::Sphere2D::DataSpectral lhs = i_lhs;
		const sweet::Data::Sphere2D::DataSpectral rhs = i_rhs;

		sweet::Data::Sphere2D::DataGrid diff = i_lhs.toGrid()-rhs.toGrid();

		double lhs_maxabs = lhs.toGrid().grid_reduce_max_abs();
		double rhs_maxabs = rhs.toGrid().grid_reduce_max_abs();

		double normalize_fac = 1.0;

		if (i_normalization)
		{
			normalize_fac = std::max(lhs_maxabs, rhs_maxabs);

			if (normalize_fac < i_error_threshold)
			{
				std::cout << "Normalization for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				normalize_fac = 1.0;
			}
		}

		double rel_max_abs = diff.grid_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.grid_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\terror threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored" << std::endl;
			}
			else
			{
				lhs.toGrid().grid_file_write("o_error_lhs_values.csv");
				rhs.toGrid().grid_file_write("o_error_rhs_values.csv");
				(lhs-rhs).toGrid().grid_file_write("o_error_lhs_rhs_diff_spectral.csv");
				diff.grid_file_write("o_error_lhs_rhs_diff_physical.csv");

				SWEETErrorFatal("Error too high");
			}

			return true;
		}
		return false;
	}



public:
	static
	bool check(
			const sweet::Data::Sphere2DComplex::DataSpectral &i_lhs,
			const sweet::Data::Sphere2DComplex::DataSpectral &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		sweet::Data::Sphere2DComplex::DataSpectral diff = i_lhs - i_rhs;

		double lhs_maxabs = sweet::Data::Sphere2DComplex::DataSpectral(i_lhs).toGrid().grid_reduce_max_abs();
		double rhs_maxabs = sweet::Data::Sphere2DComplex::DataSpectral(i_rhs).toGrid().grid_reduce_max_abs();

		double normalize_fac = 1.0;

		if (i_normalization)
		{
			normalize_fac = std::max(lhs_maxabs, rhs_maxabs);

			if (normalize_fac < i_error_threshold)
			{
				std::cout << "Normalization for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				normalize_fac = 1.0;
			}
		}

		double rel_max_abs = diff.toGrid().grid_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.toGrid().grid_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\terror threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored" << std::endl;
			}
			else
			{
				sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(i_lhs).toGrid().grid_file_write("o_error_lhs_values.csv");
				sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(i_rhs).toGrid().grid_file_write("o_error_rhs_values.csv");
				sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(i_lhs-i_rhs).toGrid().grid_file_write("o_error_lhs_rhs_diff_spectral.csv");
				sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(diff).toGrid().grid_file_write("o_error_lhs_rhs_diff_physical.csv");

				SWEETErrorFatal("Error too high");
			}

			return true;
		}
		return false;
	}



public:
	static
	bool checkTruncated(
			const sweet::Data::Sphere2D::DataSpectral &i_lhs,
			const sweet::Data::Sphere2D::DataSpectral &i_rhs,
			const sweet::Data::Sphere2D::Config *i_sphere2DDataConfig,
			const std::string &i_id,
			double i_error_threshold,	// = 1.0,
			double i_ignore_error,		// = false,
			bool i_normalization		// = true
	)
	{
		sweet::Data::Sphere2D::DataSpectral lhsr = sweet::Data::Sphere2D::DataSpectral(i_lhs).spectral_returnWithDifferentModes(i_sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral rhsr = sweet::Data::Sphere2D::DataSpectral(i_rhs).spectral_returnWithDifferentModes(i_sphere2DDataConfig);

		sweet::Data::Sphere2D::DataGrid diff = lhsr.toGrid()-rhsr.toGrid();

		double lhs_maxabs = sweet::Data::Sphere2D::DataSpectral(lhsr).toGrid().grid_reduce_max_abs();
		double rhs_maxabs = sweet::Data::Sphere2D::DataSpectral(rhsr).toGrid().grid_reduce_max_abs();

		double normalize_fac = std::min(lhs_maxabs, rhs_maxabs);

		if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
		{
			std::cout << "Error computation for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
			return false;
		}


		double rel_max_abs = diff.grid_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.grid_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tError threshold: " << i_error_threshold << "\tNormalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
				return true;
			}

			lhsr.toGrid().grid_file_write("o_error_lhs.csv");
			rhsr.toGrid().grid_file_write("o_error_rhs.csv");
			(lhsr-rhsr).toGrid().grid_file_write("o_error_lhs_rhs_diff_spectral.csv");
			diff.grid_file_write("o_error_lhs_rhs_diff_physical.csv");

			SWEETErrorFatal("Error too high");
			return true;
		}
		return false;
	}


public:
	static
	bool checkTruncated(
			const sweet::Data::Sphere2DComplex::DataSpectral &i_lhs,
			const sweet::Data::Sphere2DComplex::DataSpectral &i_rhs,
			const sweet::Data::Sphere2D::Config *i_sphere2DDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		sweet::Data::Sphere2DComplex::DataSpectral lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphere2DDataConfig);
		sweet::Data::Sphere2DComplex::DataSpectral rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphere2DDataConfig);

		sweet::Data::Sphere2DComplex::DataSpectral diff = lhsr - rhsr;
//				Convert_Sphere2DDataSpectralComplex_2_Sphere2DDataSpectral::grid_convert_real(lhsr)
//				- Convert_Sphere2DDataSpectralComplex_2_Sphere2DDataSpectral::grid_convert_real(rhsr);

		double normalize_fac;

		if (i_normalization)
		{
			double lhs_maxabs = lhsr.toGrid().grid_reduce_max_abs();
			double rhs_maxabs = rhsr.toGrid().grid_reduce_max_abs();

			if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
			{
				std::cout << "Error for " << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				return false;
			}

			normalize_fac = std::min(lhsr.toGrid().grid_reduce_max_abs(), rhsr.toGrid().grid_reduce_max_abs());

			if (normalize_fac == 0)
			{
				std::cout << "Error for " << i_id << "' ignored since at least one field is Zero" << std::endl;
				return false;
			}
		}
		else
		{
			normalize_fac = 1.0;
		}

		double rel_max_abs = diff.toGrid().grid_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.toGrid().grid_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tError threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
				return false;
			}

			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(lhsr).toGrid().grid_file_write("o_error_lhs.csv");
			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(rhsr).toGrid().grid_file_write("o_error_rhs.csv");
			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(lhsr-rhsr).toGrid().grid_file_write("o_error_lhs_rhs_diff_spectral.csv");
			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(diff).toGrid().grid_file_write("o_error_lhs_rhs_diff_physical.csv");

			SWEETErrorFatal("Error too high");
			return true;
		}
		return false;
	}



public:
	static
	bool checkTruncatedSpectral(
			const sweet::Data::Sphere2DComplex::DataSpectral &i_lhs,
			const sweet::Data::Sphere2DComplex::DataSpectral &i_rhs,
			const sweet::Data::Sphere2D::Config *i_sphere2DDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		sweet::Data::Sphere2DComplex::DataSpectral lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphere2DDataConfig);
		sweet::Data::Sphere2DComplex::DataSpectral rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphere2DDataConfig);

		sweet::Data::Sphere2DComplex::DataSpectral diff = lhsr-rhsr;

		double normalize_fac = 1.0;

		if (i_normalization)
		{
			double lhs_maxabs = lhsr.toGrid().grid_reduce_max_abs();
			double rhs_maxabs = rhsr.toGrid().grid_reduce_max_abs();

			if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
			{
				std::cout << "Error for " << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				return false;
			}

			normalize_fac = std::min(lhsr.toGrid().grid_reduce_max_abs(), rhsr.toGrid().grid_reduce_max_abs());

			if (normalize_fac == 0)
			{
				std::cout << "Error for " << i_id << "' ignored since at least one field is Zero" << std::endl;
				return false;
			}
		}

		double rel_max_abs = diff.toGrid().grid_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.toGrid().grid_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tError threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
				return false;
			}

			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(lhsr).toGrid().grid_file_write("o_error_lhs.csv");
			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(rhsr).toGrid().grid_file_write("o_error_rhs.csv");
			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(lhsr-rhsr).toGrid().grid_file_write("o_error_lhs_rhs_diff_spectral.csv");
			sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(diff).toGrid().grid_file_write("o_error_lhs_rhs_diff_physical.csv");

			SWEETErrorFatal("Error too high");
			return true;
		}
		return false;
	}

};


void run_tests(
		sweet::Data::Sphere2D::Config *sphere2DDataConfig
)
{
	double epsilon = 1e-12;
	epsilon *= (sphere2DDataConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;

	std::cout << std::setprecision(10);

	{
		sweet::Data::Sphere2D::DataGrid physical(sphere2DDataConfig);
		physical.grid_update_lambda_cogaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(2.0)*0.5;
			}
		);

		sweet::Data::Sphere2D::DataSpectral spectral(sphere2DDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 0 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		Sphere2DDataErrorCheck::check(sweet::Data::Sphere2D::DataSpectral(physical), spectral, "n=0, m=0", epsilon, false, true);
	}

	{
		sweet::Data::Sphere2D::DataGrid physical(sphere2DDataConfig);
		physical.grid_update_lambda_gaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(6.0)*mu*0.5;
			}
		);

		sweet::Data::Sphere2D::DataSpectral spectral(sphere2DDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 1 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		sweet::Data::Sphere2D::DataSpectral(physical).spectral_print();
		spectral.spectral_print();

		Sphere2DDataErrorCheck::check(sweet::Data::Sphere2D::DataSpectral(physical), spectral, "n=1, m=0", epsilon, false, true);
	}

	{
		sweet::Data::Sphere2D::DataGrid physical(sphere2DDataConfig);
		physical.grid_update_lambda_gaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(10.0)/4.0 * (3.0*mu*mu - 1.0);
			}
		);

		sweet::Data::Sphere2D::DataSpectral spectral(sphere2DDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 2 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		Sphere2DDataErrorCheck::check(sweet::Data::Sphere2D::DataSpectral(physical), spectral, "n=2, m=0", epsilon, false, true);
	}


	{
		sweet::Data::Sphere2D::DataGrid physical(sphere2DDataConfig);
		physical.grid_update_lambda_gaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(14.0)/4.0*mu * (5.0*mu*mu - 3.0);
			}
		);

		sweet::Data::Sphere2D::DataSpectral spectral(sphere2DDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 3 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		Sphere2DDataErrorCheck::check(sweet::Data::Sphere2D::DataSpectral(physical), spectral, "n=3, m=0", epsilon, false, true);
	}
}




int main(
		int i_argc,
		char *const i_argv[]
)
{

	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	if (shackSphere2DDataOps->space_res_spectral[0] == 0)
		SWEETErrorFatal("Set number of spectral modes to use SPH!");

	sweet::Data::Sphere2D::Config sphere2DDataConfig;
	sphere2DDataConfig.setupAuto(shackSphere2DDataOps);

	run_tests(&sphere2DDataConfig);
}



#endif
