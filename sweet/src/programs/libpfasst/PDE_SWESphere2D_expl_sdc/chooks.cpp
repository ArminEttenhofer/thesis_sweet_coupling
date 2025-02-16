#include <mpi.h>
#include "../interface/Sphere2DDataVars.hpp"
#include "Sphere2DDataCtxSDC.hpp"
#include "ceval.hpp"

extern "C"
{
	void cecho_error(sweet::Data::Sphere2D::DataSpectral *sd,
					 int step)
	{
		// not implemented
	}

	void cecho_residual(Sphere2DDataCtxSDC *i_ctx,
						double i_norm,
						int i_current_proc)
	{
		// get the residual vector
		std::vector<std::vector<double>> &residuals = i_ctx->get_residuals();

		// save the residual
		residuals[i_current_proc].push_back(i_norm);
	}

	void cecho_output_invariants(Sphere2DDataCtxSDC *i_ctx,
								 Sphere2DDataVars *i_Y,
								 int i_current_proc,
								 int i_current_step,
								 int i_current_iter,
								 int i_nnodes,
								 int i_niters)
	{
		const sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = i_Y->get_phi_pert();
		const sweet::Data::Sphere2D::DataSpectral &vrt_Y = i_Y->get_vrt();
		const sweet::Data::Sphere2D::DataSpectral &div_Y = i_Y->get_div();

		// get the sweet::Data::Sphere2D::Sphere2DOperators object from context
		sweet::Data::Sphere2D::Operators *sphere2DOperators = i_ctx->get_sphere2d_operators();

		// compute the invariants
		i_ctx->diagnostics.update_phi_vrt_div_2_mass_energy_enstrophy(
			sphere2DOperators,
			phi_pert_Y,
			vrt_Y,
			div_Y,
			i_ctx->shackSphere2DDataOps->sphere2d_radius,
			i_ctx->shackPDESWESphere2D->gravitation);

		std::cout << "[MULE] libpfasst.mass_s" << std::setfill('0') << std::setw(5) << i_current_step;
		std::cout << " = " << std::setprecision(20) << i_ctx->diagnostics.total_mass << std::endl;
		std::cout << "[MULE] libpfasst.energy_s" << std::setfill('0') << std::setw(5) << i_current_step;
		std::cout << " = " << std::setprecision(20) << i_ctx->diagnostics.total_energy << std::endl;
		std::cout << "[MULE] libpfasst.potential_enstrophy_s" << std::setfill('0') << std::setw(5) << i_current_step;
		std::cout << " = " << std::setprecision(20) << i_ctx->diagnostics.total_potential_enstrophy << std::endl;

		// save the invariants for plotting at the end
		i_ctx->save_grid_invariants(i_current_step);
	}

	void cecho_output_solution(Sphere2DDataCtxSDC *i_ctx,
							   Sphere2DDataVars *i_Y,
							   int i_current_proc,
							   int i_current_step,
							   int i_current_iter,
							   int i_nnodes,
							   int i_niters)
	{
		// get the pointer to the Simulation Variables object
		// sweet::Shacks::ShackDictionary *shackDict = i_ctx->get_simulation_variables();

		// update timecontrol information
		i_ctx->shackTimestepControl->currentTimestepNr = i_current_step + 1;
		auto current_dt = i_ctx->shackTimestepControl->currentTimestepSize;
		i_ctx->shackTimestepControl->currentSimulationTime = (i_current_step + 1) * current_dt;

		// check if we should write output
		if (!timestep_check_output(i_ctx, i_current_iter, i_niters))
		{
			return;
		}

		// update when to write output the next time
		i_ctx->shackIOData->_outputNextSimTime += i_ctx->shackIOData->outputEachSimTime;

		const sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = i_Y->get_phi_pert();
		const sweet::Data::Sphere2D::DataSpectral &vrt_Y = i_Y->get_vrt();
		const sweet::Data::Sphere2D::DataSpectral &div_Y = i_Y->get_div();

		// write the data to file
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank != 0)
		{
			return;
		}
		std::string filename = "prog_phi_pert";
		write_file(i_ctx, phi_pert_Y, filename.c_str());

		filename = "prog_vrt";
		write_file(i_ctx, vrt_Y, filename.c_str());

		filename = "prog_div";
		write_file(i_ctx, div_Y, filename.c_str());
	}
}
