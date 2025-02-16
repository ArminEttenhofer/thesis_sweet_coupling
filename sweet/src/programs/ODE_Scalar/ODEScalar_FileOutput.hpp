/*
 * Author: JOAO STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_SCALAR_ODESCALAR_FILEOUTPUT_HPP
#define PROGRAMS_ODE_SCALAR_ODESCALAR_FILEOUTPUT_HPP

#include <programs/ODE_Scalar/ShackODEScalar.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <fstream>

// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>



class ODEScalar_FileOutput
{
public:
	sweet::Error::Base error;

	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackODEScalar *shackODEScalar;

	void setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl,
			ShackODEScalar *i_shackODEScalar
	)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
		shackODEScalar = i_shackODEScalar;
	}

	void clear()
	{
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackODEScalar = nullptr;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const double u,
			const char* i_name		//!< name of output variable
	)
	{
		char buffer[1024];

		const char* filename_template = shackIOData->outputFileName.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);

		std::ofstream file(buffer, std::ios_base::trunc);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE SCALAR" << std::endl;

		file << std::setprecision(16);
		file << u;

		file.close();

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output(
			double u
	)
	{
		if (shackIOData->outputFileName.length() == 0)
			return;

		std::cout << "Writing output files at simulation time: " << shackTimestepControl->currentSimulationTime << " secs" << std::endl;

		if (shackIOData->outputFileMode == "csv")
		{
			std::string output_filename;

			output_filename = write_file_csv(u, "prog_u");
			output_reference_filenames += ";"+output_filename;
		}
		else if (shackIOData->outputFileMode == "bin")
		{
			SWEETErrorFatal("Bin output not available for ODEScalar");
		}
		else
		{
			SWEETErrorFatal("Unknown output file mode '"+shackIOData->outputFileMode+"'");
		}
	}
};

#endif
