/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_FILEOUTPUT_HPP
#define PROGRAMS_PDE_SWESPHERE2D_FILEOUTPUT_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#include "Shack.hpp"

namespace PDE_SWESphere2D {

class FileOutput
{
public:
	sweet::Error::Base error;

	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	Shack *shackPDESWESphere2D;

	void setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl,
			Shack *i_shackPDESWESphere2D
	)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
		shackPDESWESphere2D = i_shackPDESWESphere2D;
	}

	void clear()
	{
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackPDESWESphere2D = nullptr;
	}

	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv_spec_evol(
			const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
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
		i_sphere2DData.spectrum_phase_file_write_line(buffer, 
			i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale,
			20, 10e-20, reduce_mode_factor);

		sprintf(buffer, filename_template_ampl, i_name);
		i_sphere2DData.spectrum_abs_file_write_line(buffer, 
			i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale,
			20, 10e-20, reduce_mode_factor);

		return buffer;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
			const char* i_name,		//!< name of output variable
			bool i_phi_shifted = false
	)
	{
		char buffer[1024];

		// create copy
		sweet::Data::Sphere2D::DataGrid sphere2DData = i_sphere2DData.toGrid();

		const char* filename_template = shackIOData->outputFileName.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);

		if (i_phi_shifted)
			sphere2DData.grid_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphere2DData.grid_file_write(buffer);

		return buffer;
	}



	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin(
			const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
			const char* i_name
	)
	{
		char buffer[1024];

		sweet::Data::Sphere2D::DataSpectral sphere2DData(i_sphere2DData);
		const char* filename_template = shackIOData->outputFileName.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);
		sphere2DData.file_write_binary_spectral(buffer);

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output(
			sweet::Data::Sphere2D::Operators &i_ops,
			sweet::Data::Sphere2D::DataSpectral &i_prog_phi_pert,
			sweet::Data::Sphere2D::DataSpectral &i_prog_div,
			sweet::Data::Sphere2D::DataSpectral &i_prog_vrt
	)
	{
		if (shackIOData->outputFileName.length() == 0)
			return;

		std::cout << "Writing output files at simulation time: " << shackTimestepControl->currentSimulationTime << " secs" << std::endl;

		if (shackIOData->outputFileMode == "csv")
		{
			std::string output_filename;

			sweet::Data::Sphere2D::DataSpectral h = i_prog_phi_pert*(1.0/shackPDESWESphere2D->gravitation);
			h += shackPDESWESphere2D->h0;

			output_filename = write_file_csv(h, "prog_h");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << " (min: " << h.toGrid().grid_reduce_min() << ", max: " << h.toGrid().grid_reduce_max() << ")" << std::endl;

			output_filename = write_file_csv(i_prog_phi_pert, "prog_phi_pert");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << i_prog_phi_pert.toGrid().grid_reduce_min() << ", max: " << i_prog_phi_pert.toGrid().grid_reduce_max() << ")" << std::endl;

			sweet::Data::Sphere2D::DataGrid phi_phys = h.toGrid() * shackPDESWESphere2D->gravitation;
			sweet::Data::Sphere2D::DataSpectral phi(i_ops.sphere2DDataConfig);
			phi.loadSphere2DDataGrid(phi_phys);
			output_filename = write_file_csv(phi, "prog_phi");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << phi_phys.grid_reduce_min() << ", max: " << phi_phys.grid_reduce_max() << ")" << std::endl;

			sweet::Data::Sphere2D::DataGrid u(i_ops.sphere2DDataConfig);
			sweet::Data::Sphere2D::DataGrid v(i_ops.sphere2DDataConfig);

			i_ops.vrtdiv_2_uv(i_prog_vrt, i_prog_div, u, v);

			output_filename = write_file_csv(u, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(v, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(i_prog_vrt, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(i_prog_div, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			sweet::Data::Sphere2D::DataSpectral potvrt = (i_prog_phi_pert/shackPDESWESphere2D->gravitation)*i_prog_vrt;

			output_filename = write_file_csv(potvrt, "prog_potvrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;
		}
		else if (shackIOData->outputFileMode == "bin")
		{
			std::string output_filename;

            {
                output_filename = write_file_bin(i_prog_phi_pert * (1.0/shackPDESWESphere2D->gravitation), "prog_h_pert");
                output_reference_filenames = output_filename;
                sweet::Data::Sphere2D::DataGrid prog_h = i_prog_phi_pert.toGrid() * (1.0/shackPDESWESphere2D->gravitation);

                std::cout << " + " << output_filename << " (min: " << prog_h.grid_reduce_min() << ", max: " << prog_h.grid_reduce_max() << ")" << std::endl;
            }

			{
				output_filename = write_file_bin(i_prog_phi_pert, "prog_phi_pert");
				output_reference_filenames = output_filename;
				sweet::Data::Sphere2D::DataGrid prog_phys = i_prog_phi_pert.toGrid();

				std::cout << " + " << output_filename << " (min: " << prog_phys.grid_reduce_min() << ", max: " << prog_phys.grid_reduce_max() << ")" << std::endl;
			}


			{
                // Recover U, V output
                sweet::Data::Sphere2D::DataGrid u{i_prog_vrt.sphere2DDataConfig};
                sweet::Data::Sphere2D::DataGrid v{i_prog_vrt.sphere2DDataConfig};

                i_ops.vrtdiv_2_uv(i_prog_vrt, i_prog_div, u, v);

                sweet::Data::Sphere2D::DataSpectral u_spec =  i_ops.scalar_grid_2_spectral(u);
                sweet::Data::Sphere2D::DataSpectral v_spec =   i_ops.scalar_grid_2_spectral(v);

				output_filename = write_file_bin(u_spec, "prog_u");
                output_reference_filenames += ";"+output_filename;
                output_filename = write_file_bin(v_spec, "prog_v");
                output_reference_filenames += ";"+output_filename;
			}
//
//			{
//				output_filename = write_file_bin(i_prog_div, "prog_div");
//				output_reference_filenames += ";"+output_filename;
//				sweet::Data::Sphere2D::DataGrid prog_phys = i_prog_div.toGrid();
//
//				std::cout << " + " << output_filename << " (min: " << prog_phys.grid_reduce_min() << ", max: " << prog_phys.grid_reduce_max() << ")" << std::endl;
//			}
		}
		else if (shackIOData->outputFileMode == "csv_spec_evol"){

			std::string output_filename;

			{ 
				/*
				* Spectral kinetic energy and potential enstrophy calculation and output
				*
				* Details in Jakob-Chien, Ruediger, James J. Hack, and David L. Williamson. 
				* "Spectral transform solutions to the shallow water test set." Journal of Computational Physics 119, no. 1 (1995): 164-187.
				*/
				// Kinetic energy is given in spectral space as
				// KE per mode = a^2/((n(n+1)))*(vrt*conj(vrt))+a^2/((n(n+1)))*(div*conj(div))
				// r = a/(sqrt(n(n+1))) (root_laplace)
				// KE per mode = (r*vrt*conj(r*vrt))+(r*div*conj(r*div))
				sweet::Data::Sphere2D::DataSpectral rlap_vrt = i_ops.inv_root_laplace(i_prog_vrt);
				sweet::Data::Sphere2D::DataSpectral rlap_div = i_ops.inv_root_laplace(i_prog_div);
				sweet::Data::Sphere2D::DataSpectral kin_en = rlap_vrt + rlap_div ;

				output_filename = write_file_csv_spec_evol(kin_en, "kin_en"); 
				std::cout << " + " << output_filename << " (Total Kin Energy : " << 0.25*kin_en.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// For Barotropic vort eq: See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere2D (2009) for eps=0
				// Kinetic energy is given in spectral space as
				// Vortical energy per mode is (0.5 n*(n+1) / a^2) *psi*conj(psi) in spectral space
				//Sphere2DData_Spectral psi = op.inv_laplace(prog_vrt); // 
				// multiply psi by sqrt( n * (n+1))/a (apply root laplacian)
				//Sphere2DData_Spectral psi_root = op.root_laplace(psi);
				//output_filename = write_file_csv_spec_evol(psi_root*std::sqrt(0.5), "spec_energy"); 
				//std::cout << " + " << output_filename << " (Kinetic energy : " << (0.5)*psi_root.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere2D (2009) for eps=0
				// enstrophy per mode is 0.5 vrt*conj(vrt) in spectral space
				// Total enstrophy is the sum of these (counting twice modes with m>0 and once when m=0)
				output_filename = write_file_csv_spec_evol(i_prog_vrt, "enstrophy");
				std::cout << " + " << output_filename << " (Total Enstrophy : " << i_prog_vrt.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

			}
		}
		else
		{
			SWEETErrorFatal("Unknown output file mode '"+shackIOData->outputFileMode+"'");
		}
	}

	bool fileSave(
			sweet::Data::Sphere2D::DataSpectral &i_U,
			const std::string &i_name
	)
	{

		if (shackIOData->outputFileMode == "csv")
			write_file_csv(i_U, i_name.c_str());
		else if (shackIOData->outputFileMode == "bin")
			write_file_bin(i_U, i_name.c_str());
		else if (shackIOData->outputFileMode == "csv_spec_evol")
			write_file_csv_spec_evol(i_U, i_name.c_str());
		else
			SWEETErrorFatal("Unknown output file mode '"+shackIOData->outputFileMode+"'");

		return true;
	}


};

}

#endif
