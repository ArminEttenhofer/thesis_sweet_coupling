/*
 * Error.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef INCLUDE_SWEET_XBRAID_ERROR_HPP
#define INCLUDE_SWEET_XBRAID_ERROR_HPP

#if SWEET_GUI
#include<sweet/GUI/VisSweet.hpp>
#endif

#include <algorithm>
#include <map>

#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/XBraid/Shack.hpp>

#include <sweet/Data/GenericContainer/ConfigBase.hpp>

#include "Vector.hpp"

namespace sweet {
namespace XBraid {

class Error
{

	sweet::IO::Shack* shackIOData = nullptr;

public:
	Error()
	{
	}

public:
	~Error()
	{
		clear();
	}

private:
	void clear()
	{
		shackIOData = nullptr;
	}

public:
	bool setup(
			sweet::IO::Shack *i_shackIOData
		)
	{
		shackIOData = i_shackIOData;
		return true;
	}


	/**
	 * Load reference solution
	 */
	void loadRefSolution(
			sweet::XBraid::Vector* pint_data_ref,
			double t,
			std::string path_ref
	)
	{

		pint_data_ref->fileLoad(
					shackIOData->outputFileName.c_str(),
					path_ref,
					shackIOData->outputFileMode,
					t,
					shackIOData->outputFormatTimeScale
		);

	}

	/**
	 * Compute and store parareal errors during simulation
	 */
	void computeStoreError(
			sweet::XBraid::Vector* i_diff,
			sweet::XBraid::Vector* i_data_ref,
			int i_iteration_id,
			int i_time_slice_id,
			double t,
			std::string i_path_ref,
			std::string i_base_solution,	// "ref" or "fine"
			int i_precision = 32
	)
	{

		double err_L1; // physical space
		double err_L2; // physical space
		double err_Linf; // physical space

		std::map<std::size_t, double> err_Linf_spectral;
		std::vector<std::size_t> rnorms = getRnorms();

		for (int ivar = 0; ivar < (int)i_diff->N; ivar++)
		{

			std::string var_name = i_diff->var_names[ivar];

			err_L1 = i_diff->reduceNormL1Grid(ivar, true);
			err_L2 = i_diff->reduceNormL2Grid(ivar, true);
			err_Linf = i_diff->reduceNormLinfGrid(ivar);

			// Spectral space
			double small = 1e-20;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
			{
				double norm_diff = std::sqrt(i_diff->reduceMaxAbs(ivar, *it) );
				double norm_ref = std::sqrt(i_data_ref->reduceMaxAbs(ivar, *it) );
				if (norm_diff < small and norm_ref < small)
					err_Linf_spectral.emplace(std::make_pair(*it, 0.));
				else
					err_Linf_spectral.emplace(std::make_pair(*it, norm_diff / norm_ref ));
			}

			// save physical errors in file
			char buffer_out[1024];
			std::string str = "xbraid_error_%s_%s_t%020.8f_iter%03d.csv";
			const char* filename_template_out = str.c_str();
			sprintf(buffer_out, filename_template_out, i_base_solution.c_str(), var_name.c_str(), t * shackIOData->outputFormatTimeScale, i_iteration_id);

			std::ofstream file(buffer_out, std::ios_base::trunc);
			file << std::setprecision(i_precision);

			file << "#BASESOLUTION " << i_base_solution << " " << i_path_ref << std::endl;
			file << "#VAR " << var_name << std::endl;
			file << "#ITERATION " << i_iteration_id << std::endl;
			file << "#TIMESLICE " << i_time_slice_id << std::endl;
			file << "#TIMEFRAMEEND " << t  * shackIOData->outputFormatTimeScale << std::endl;
			file << "errL1 " << err_L1 << std::endl;
			file << "errL2 " << err_L2 << std::endl;
			file << "errLinf " << err_Linf << std::endl;

			file.close();

			// save spectral errors in file
			char buffer_out_spec[1024];
			std::string str2 = "xbraid_error_spec_%s_%s_t%020.8f_iter%03d.csv";
			const char* filename_template_out_spec = str2.c_str();
			sprintf(buffer_out_spec, filename_template_out_spec, i_base_solution.c_str(), var_name.c_str(), t * shackIOData->outputFormatTimeScale, i_iteration_id);

			std::ofstream file_spec(buffer_out_spec, std::ios_base::trunc);
			file_spec << std::setprecision(i_precision);

			file_spec << "#BASESOLUTION " << i_base_solution << " " << i_path_ref << std::endl;
			file_spec << "#VAR " << var_name << std::endl;
			file_spec << "#ITERATION " << i_iteration_id << std::endl;
			file_spec << "#TIMESLICE " << i_time_slice_id << std::endl;
			file_spec << "#TIMEFRAMEEND " << t  * shackIOData->outputFormatTimeScale << std::endl;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
				file_spec << "errLinf " << *it + 1 << " " << err_Linf_spectral.at(*it) << std::endl;

			file_spec.close();

		}
	}

	virtual
	std::vector<size_t> getRnorms() = 0;
};

}}

#endif
