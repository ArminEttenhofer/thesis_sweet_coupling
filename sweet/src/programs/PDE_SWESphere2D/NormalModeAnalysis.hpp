/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_NORMALMODEANALYSIS_HPP
#define PROGRAMS_PDE_SWESPHERE2D_NORMALMODEANALYSIS_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

namespace PDE_SWESphere2D {

class NormalModeAnalysis
{
public:
	template <typename TCallbackClass>
	static
	void normal_mode_analysis(
			sweet::Data::Sphere2D::DataSpectral &io_prog_phi,
			sweet::Data::Sphere2D::DataSpectral &io_prog_vort,
			sweet::Data::Sphere2D::DataSpectral &io_prog_div,

			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl,

			int mode,


			TCallbackClass *i_class,
			bool(TCallbackClass::* const i_run_timestep_method)(void)
	)
	{
		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_prog_phi.sphere2DDataConfig;

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		char buffer_real[1024];
		const char* filename = i_shackIOData->outputFileName.c_str();

		if (i_shackTimestepControl->maxTimestepsNr > 0)
			sprintf(buffer_real, filename, "normal_modes_physical", i_shackTimestepControl->currentTimestepSize*i_shackTimestepControl->maxTimestepsNr*i_shackIOData->outputFormatTimeScale);
		else
			sprintf(buffer_real, filename, "normal_modes_physical", i_shackTimestepControl->currentTimestepSize*i_shackIOData->outputFormatTimeScale);

		std::ofstream file(buffer_real, std::ios_base::trunc);
		std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

		std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

		// use very high precision
		file << std::setprecision(20);

		sweet::Data::Sphere2D::DataSpectral* prog[3] = {&io_prog_phi, &io_prog_vort, &io_prog_div};

		int max_prog_id = 3;
		prog[0]->spectral_setZero();
		prog[1]->spectral_setZero();
		prog[2]->spectral_setZero();

		int num_timesteps = 1;
		if (mode >= 10)
		{
			if (i_shackTimestepControl->maxTimestepsNr > 0)
				num_timesteps = i_shackTimestepControl->maxTimestepsNr;
		}

		if (i_shackTimestepControl->maxSimulationTime > 0)
			file << "# t " << i_shackTimestepControl->maxSimulationTime << std::endl;
		else
			file << "# t " << (num_timesteps*i_shackTimestepControl->currentTimestepSize) << std::endl;

#if 0
		file << "# g " << i_shackDict.sim.gravitation << std::endl;
		file << "# h " << i_shackDict.sim.h0 << std::endl;
		file << "# r " << i_shackDict.sim.sphere2d_radius << std::endl;
		file << "# f " << i_shackDict.sim.sphere2d_rotating_coriolis_omega << std::endl;
#endif

		// iterate over all prognostic variables
		for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
		{
			std::cout << "normal mode analysis for prog " << outer_prog_id << std::endl;

			if (mode == 1 || mode == 11)
			{
				SWEETErrorFatal("Not supported anymore");
#if 0
				// iterate over physical space
				for (int outer_i = 0; outer_i < sphere2DDataConfig->grid_number_elements; outer_i++)
				{
					// reset time control
					i_shackTimestepControl->currentTimestepNr = 0;
					i_shackTimestepControl->currentSimulationTime = 0;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->spectral_setZero();

					// activate mode
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->grid_space_data[outer_i] = 1;

					/*
					 * RUN timestep
					 */

					// In case of a multi-step scheme, reset it!
					if (i_shackDict.disc.timestepping_method.find("_lf") != std::string::npos)
					{
						SWEETErrorFatal("TODO 01943934");
						//sphere2Ddata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, i_shackDict.disc.timestepping_order, i_shackDict.disc.leapfrog_robert_asselin_filter);

						i_shackTimestepControl->currentTimestepSize = leapfrog_start_timesteps_size;

						for (int i = 0; i < leapfrog_start_num_timesteps; i++)
							(i_class->*i_run_timestep_method)();

						i_shackTimestepControl->currentTimestepSize = leapfrog_end_timestepSize;

						(i_class->*i_run_timestep_method)();

						i_shackTimestepControl->currentTimestepSize = leapfrog_original_timestepSize;
					}
					else
					{
						(i_class->*i_run_timestep_method)();
					}

					for (int i = 1; i < num_timesteps; i++)
					{
						(i_class->*i_run_timestep_method)();
					}


					if (mode == 1)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->grid_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/i_shackTimestepControl->currentTimestepSize);
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_physical();
						for (int k = 0; k < sphere2DDataConfig->grid_number_elements; k++)
						{
							file << prog[inner_prog_id]->grid_space_data[k];
							if (inner_prog_id != max_prog_id-1 || k != sphere2DDataConfig->grid_number_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
#endif
			}
			else if (mode == 2 || mode == 12)
			{

				// iterate over physical space
				for (int outer_i = 0; outer_i < sphere2DDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					for (int imag_i = 0; imag_i < 2; imag_i++)
					{
						// reset time control
						i_shackTimestepControl->currentTimestepNr = 0;
						i_shackTimestepControl->currentSimulationTime = 0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->spectral_setZero();

						// activate mode
						if (imag_i)
							prog[outer_prog_id]->spectral_space_data[outer_i].imag(1);
						else
							prog[outer_prog_id]->spectral_space_data[outer_i].real(1);

						for (int i = 0; i < num_timesteps; i++)
							(i_class->*i_run_timestep_method)();


						if (mode == 2)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								prog[inner_prog_id]->operator*=(1.0/i_shackTimestepControl->currentTimestepSize);
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							for (int k = 0; k < sphere2DDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].real();
								file << "\t";
							}
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							for (int k = 0; k < sphere2DDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].imag();
								if (inner_prog_id != max_prog_id-1 || k != sphere2DDataConfig->spectral_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
			}
			else if (mode == 3 || mode == 13)
			{
				// iterate over spectral space
				for (int outer_i = 0; outer_i < sphere2DDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					i_shackTimestepControl->currentTimestepNr = 0;
					i_shackTimestepControl->currentSimulationTime = 0;

					std::cout << "." << std::flush;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->spectral_setZero();

					// activate mode via real coefficient
					prog[outer_prog_id]->spectral_space_data[outer_i].real(1);


					(i_class->*i_run_timestep_method)();

					for (int i = 0; i < num_timesteps; i++)
					{
						(i_class->*i_run_timestep_method)();
					}

					if (mode == 3)
					{
						/*
						 * Compute
						 *    1/dt * (U(t+1) - U(t))
						 * for linearization
						 */
						prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/i_shackTimestepControl->currentTimestepSize);
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						for (int k = 0; k < sphere2DDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].real();
							file << "\t";
						}
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						for (int k = 0; k < sphere2DDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].imag();

							if (inner_prog_id != max_prog_id-1 || k != sphere2DDataConfig->spectral_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
		}
	}
};

}

#endif
