/*
 * FileOutput.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef INCLUDE_SWEET_XBRAID_FILEOUTPUT_HPP
#define INCLUDE_SWEET_XBRAID_FILEOUTPUT_HPP

////#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
////#include <sweet/_DEPRECATED_pint/PInT_Common.hpp>

///#include <sweet/_DEPRECATED_pint/PInT_Common.hpp>
///#include <xbraid/braid.hpp>

#if SWEET_GUI
#include<sweet/GUI/VisSweet.hpp>
#endif

#include <algorithm>
#include <fstream>

#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/XBraid/Shack.hpp>

#include <sweet/Data/GenericContainer/ConfigBase.hpp>

///#include "GeometryDependentDefinitions.hpp"
#include "Vector.hpp"

namespace sweet {
namespace XBraid {

class FileOutput
{

protected:
	sweet::IO::Shack* shackIOData = nullptr;
	sweet::TimeTree::Shack* shackTimestepControl = nullptr;

public:
	FileOutput()
	{
	}

public:
	~FileOutput()
	{
		clear();
	}

public:
	void clear()
	{
	shackIOData = nullptr;
	shackTimestepControl = nullptr;
	}

public:
	bool setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl
		)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
		return true;
	}


	/**
	 * Output XBraid residual at each iteration
	 */
	void output_residual_file(
			double res,
			int iteration_id
	)
	{

		char buffer[1024];

		const char* filename_template = "residual_iter%03d.csv";
		sprintf(buffer, filename_template, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(16);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE SCALAR" << std::endl;

		file << res;

		file.close();

	}

	/**
	 * Output XBraid solution at given iteration and time
	 */
	void output_data_file(
			sweet::XBraid::Vector* i_vector,
			int i_iteration_id,
			int i_time_slice_id,
			double i_t
	)
	{

		shackTimestepControl->currentSimulationTime = i_t;

		// Store csv files
		// Include iteration in filename
		char char_iter[10];
		std::sprintf(char_iter, "%03d", i_iteration_id);
		std::string str_iter(char_iter);

		std::string filename_template;
		if (shackIOData->outputFileMode == "csv")
			filename_template = "output_%s_t%020.8f_iter" + str_iter + ".csv";
		else if (shackIOData->outputFileMode == "bin")
			filename_template = "output_%s_t%020.8f_iter" + str_iter + ".sweet";
		else
			SWEETErrorFatal("Unknown output file mode '"+shackIOData->outputFileMode+"'");
		shackIOData->outputFileName = filename_template;
		fileSave(i_vector);

	}

	virtual
	void fileSave(
			sweet::XBraid::Vector* i_vector
	) = 0;

};

}}

#endif
