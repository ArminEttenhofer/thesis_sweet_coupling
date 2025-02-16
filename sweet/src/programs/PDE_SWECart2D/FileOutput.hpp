/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_PDE_SWECART2D_FILEOUTPUT_HPP
#define PROGRAMS_PDE_SWECART2D_FILEOUTPUT_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
/////#include <sweet/Data/Cart2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>
/////#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#include "Shack.hpp"

namespace PDE_SWECart2D {

class FileOutput
{
public:
	sweet::Error::Base error;

	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	Shack *shackPDESWECart2D;
	sweet::Data::Cart2D::Config *cart2DDataConfig;
	sweet::Data::Cart2D::Operators *ops;
	sweet::Data::Cart2DComplex::Operators *opsComplex;

	sweet::Data::Cart2D::GridMapping gridMapping;

	void setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl,
			sweet::Data::Cart2D::Shack *i_shackCart2DDataOps,
			Shack *i_shackPDESWECart2D,
			sweet::Data::Cart2D::Config *i_cart2DDataConfig,
			sweet::Data::Cart2D::Operators *i_ops,
			sweet::Data::Cart2DComplex::Operators *i_opsComplex
	)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
		shackCart2DDataOps = i_shackCart2DDataOps;
		shackPDESWECart2D = i_shackPDESWECart2D;
		cart2DDataConfig = i_cart2DDataConfig;
		ops = i_ops;
		i_opsComplex = opsComplex;

		if (i_shackCart2DDataOps->space_grid_use_c_staggering)
			gridMapping.setup(i_shackCart2DDataOps, cart2DDataConfig);
	}

	void clear()
	{
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackCart2DDataOps = nullptr;
		shackPDESWECart2D = nullptr;
		cart2DDataConfig = nullptr;
		ops = nullptr;
		opsComplex = nullptr;
	}

	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv_spec_evol(
			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
			const char* i_name		//!< name of output variable
	)
	{
		char buffer[1024];
		std::string phase = "_arg";
		std::string ampl = "_amp"; 

		const char* filename_template_ampl = "output_spec_ampl_%s.txt"; //.c_str();
		const char* filename_template_arg = "output_spec_arg_%s.txt"; //.c_str();
		int reduce_mode_factor = 4;

		sprintf(buffer, filename_template_arg, i_name);
		i_cart2DData.spectrum_phase_file_write_line(buffer, 
			i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale,
			20, 10e-20, reduce_mode_factor);

		sprintf(buffer, filename_template_ampl, i_name);
		i_cart2DData.spectrum_abs_file_write_line(buffer, 
			i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale,
			20, 10e-20, reduce_mode_factor);

		return buffer;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
			const char* i_name	//!< name of output variable
		)
	{
		char buffer[1024];

		// TODO: convert spectral datato physical

		const char* filename_template = shackIOData->outputFileName.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);
		i_cart2DData.toGrid().file_grid_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin(
			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
			const char* i_name
	)
	{
		char buffer[1024];

		sweet::Data::Cart2D::DataSpectral cart2DData(i_cart2DData);
		const char* filename_template = shackIOData->outputFileName.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);
		cart2DData.file_write_binary_spectral(buffer);

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output(
			sweet::Data::Cart2D::Operators &i_ops,
			sweet::Data::Cart2D::DataSpectral &i_prog_h_pert,
			sweet::Data::Cart2D::DataSpectral &i_prog_u,
			sweet::Data::Cart2D::DataSpectral &i_prog_v
	)
	{

		if (shackIOData->outputFileName.length() == 0)
			return;

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::Data::Cart2D::DataGrid t_h(i_prog_h_pert.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid t_u(i_prog_h_pert.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid t_v(i_prog_h_pert.cart2DDataConfig);

		if (shackCart2DDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
		{
			t_h = i_prog_h_pert.toGrid();
			gridMapping.mapCtoA_u(i_prog_u.toGrid(), t_u);
			gridMapping.mapCtoA_v(i_prog_v.toGrid(), t_v);
		}
		else
		{
			t_h = i_prog_h_pert.toGrid();
			t_u = i_prog_u.toGrid();
			t_v = i_prog_v.toGrid();
		}

		std::cout << "Writing output files at simulation time: " << shackTimestepControl->currentSimulationTime << " secs" << std::endl;

		std::string output_filenames;
		// Dump  data in csv, if output filename is not empty
		if (shackIOData->outputFileMode == "csv")
		/////if (shackIOData->outputFileName.size() > 0)
		{

			output_reference_filenames = write_file_csv(t_h, "prog_h_pert");
			output_reference_filenames += ";" + write_file_csv(t_u, "prog_u");
			output_reference_filenames += ";" + write_file_csv(t_v, "prog_v");

			output_reference_filenames += ";" + write_file_csv(ops->ke(t_u,t_v),"diag_ke");

#if SWEET_USE_CART2D_SPECTRAL_SPACE
			output_reference_filenames += ";" + write_file_csv_spec_evol(ops->ke(t_u,t_v),"diag_ke_spec");

			output_reference_filenames += ";" + write_file_csv_spec_evol(t_h, "prog_h_pert_spec");
			output_reference_filenames += ";" + write_file_csv_spec_evol(t_u, "prog_u_spec");
			output_reference_filenames += ";" + write_file_csv_spec_evol(t_v, "prog_v_spec");

			output_reference_filenames += ";" + write_file_csv_spec_evol(ops->ke(t_u,t_v).toGrid(), "diag_ke_spec");
#endif

			output_reference_filenames += ";" + write_file_csv(ops->vort(t_u, t_v), "diag_vort");
			output_reference_filenames += ";" + write_file_csv(ops->div(t_u, t_v), "diag_div");

		}
		else if (shackIOData->outputFileMode == "bin")
		{

			output_reference_filenames = write_file_bin(t_h, "prog_h_pert");
			output_reference_filenames += ";" + write_file_bin(t_u, "prog_u");
			output_reference_filenames += ";" + write_file_bin(t_v, "prog_v");

			output_reference_filenames += ";" + write_file_bin(ops->ke(t_u,t_v),"diag_ke");

#if SWEET_USE_CART2D_SPECTRAL_SPACE
			output_reference_filenames += ";" + write_file_csv_spec_evol(ops->ke(t_u,t_v),"diag_ke_spec");

			output_reference_filenames += ";" + write_file_csv_spec_evol(t_h, "prog_h_pert_spec");
			output_reference_filenames += ";" + write_file_csv_spec_evol(t_u, "prog_u_spec");
			output_reference_filenames += ";" + write_file_csv_spec_evol(t_v, "prog_v_spec");

			output_reference_filenames += ";" + write_file_csv_spec_evol(ops->ke(t_u,t_v).toGrid(), "diag_ke_spec");
#endif

			output_reference_filenames += ";" + write_file_bin(ops->vort(t_u, t_v), "diag_vort");
			output_reference_filenames += ";" + write_file_bin(ops->div(t_u, t_v), "diag_div");

		}
		else if (shackIOData->outputFileMode == "csv_spec_evol")
		{
		}
		else
		{
			SWEETErrorFatal("Unknown output file mode '"+shackIOData->outputFileMode+"'");
		}

	}

	bool fileSave(
			sweet::Data::Cart2D::DataSpectral &i_U,
			const std::string &i_name
	)
	{

		// TODO: use dict???????
		sweet::Data::Cart2D::DataGrid t_U(i_U.cart2DDataConfig);
		////sweet::Data::Cart2D::DataSpectral t_U_spec(i_U.cart2DDataConfig);

		if (shackCart2DDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
		{
			if (i_name == "prog_h_pert")
				t_U = i_U.toGrid();
			else if (i_name == "prog_u")
				gridMapping.mapCtoA_u(i_U.toGrid(), t_U);
			else if (i_name == "prog_v")
				gridMapping.mapCtoA_v(i_U.toGrid(), t_U);
			else
				SWEETErrorFatal("Unknown var name '"+i_name+"'");
		}
		else
			t_U = i_U.toGrid();

		/////t_U_spec.loadDataGrid(t_U);

		if (shackIOData->outputFileMode == "csv")
			write_file_csv(t_U, i_name.c_str());
		else if (shackIOData->outputFileMode == "bin")
			write_file_bin(t_U, i_name.c_str());
		else if (shackIOData->outputFileMode == "csv_spec_evol")
			write_file_csv_spec_evol(t_U, i_name.c_str());
		else
			SWEETErrorFatal("Unknown output file mode '"+shackIOData->outputFileMode+"'");

		return true;
	}


	/**
	 * Write current time step info to file
	 */

	std::string write_output_file(
			std::stringstream &buffer
		)
	{
		const char* filename_template = "output_diag_evol.txt";
		std::ofstream file(filename_template, std::ofstream::out | std::ofstream::app);
		file << std::setprecision(12);
			file << buffer.str() << std::endl;

		return buffer.str();
	}


};

}

#endif
