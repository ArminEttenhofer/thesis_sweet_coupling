/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_FILEOUTPUT_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_FILEOUTPUT_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#include <programs/ODE_Generic/FileOutput/Base.hpp>

#include "Shack.hpp"

namespace ODE_Generic {
namespace DE_Dahlquist {

class FileOutput
		: public ODE_Generic::FileOutput::Base
{
public:
	///////typedef std::complex<double> T;

	///////sweet::Error::Base error;

	///////sweet::IO::Shack *shackIOData;
	///////sweet::TimeTree::Shack *shackTimestepControl;

	///////void setup(
	///////		sweet::IO::Shack *i_shackIOData,
	///////		sweet::TimeTree::Shack *i_shackTimestepControl
	///////)
	///////{
	///////	shackIOData = i_shackIOData;
	///////	shackTimestepControl = i_shackTimestepControl;
	///////}

//////	void clear()
//////	{
//////		shackIOData = nullptr;
//////		shackTimestepControl = nullptr;
//////	}
//////
//////	/**
//////	 * Write file to data and return string of file name
//////	 */
//////	std::string write_file_csv(
//////			const T i_data,
//////			const char* i_name		//!< name of output variable
//////	)
//////	{
//////		char buffer[1024];
//////
//////		const char* filename_template = shackIOData->outputFileName.c_str();
//////		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);
//////
//////		std::ofstream file(buffer, std::ios_base::trunc);
//////		file << std::setprecision(16);
//////
//////		file << "#SWEET" << std::endl;
//////		file << "#FORMAT ASCII" << std::endl;
//////		file << "#PRIMITIVE SCALAR" << std::endl;
//////
//////		file << i_data;
//////
//////		file.close();
//////
//////		return buffer;
//////	}

//////	/**
//////	 * Write file to data and return string of file name
//////	 */
//////	std::string write_file_bin(
//////			const T i_data,
//////			const char* i_name
//////	)
//////	{
//////		return "";
//////	}

	bool fileSave(
			const T &i_U,
			///const std::string &i_filename
			const std::string &i_name
	)	const override
	{

		char buffer[1024];
		const char* filename_template = shackIOData->outputFileName.c_str();
		sprintf(buffer, filename_template, i_name.c_str(), shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);

		//SWEET_ASSERT_MSG(i_splitFiles == false, "Splitting of files not supported so far");

		sweet::Dict::Dict dict;

		// Setup Header
		dict.set("sweetMagicCode", "SWEET1505");
		dict.set("dataType", "VectorComplex");

		/*
		 * Instead of storing just a single scalar, we store it as a vector
		 * This makes it more generic
		 */
		sweet::Dict::TypesArrayND<1,std::complex<double>> vector;
		vector.resize(1);
		vector.set(0, i_U);

		// Store data
		dict.set("data", vector);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_BOOLEAN(dict);

		std::cout << "Writing to file '" << buffer << "'" << std::endl;
		// Write to file
		dict.fileSave(buffer);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_BOOLEAN(dict);

#if !SWEET_XBRAID
		if (shackTimestepControl->currentSimulationTime == shackTimestepControl->maxSimulationTime)
			std::cout << "[MULE] reference_filenames: " << buffer << std::endl;
#endif

		return true;
	}


};

}}

#endif
