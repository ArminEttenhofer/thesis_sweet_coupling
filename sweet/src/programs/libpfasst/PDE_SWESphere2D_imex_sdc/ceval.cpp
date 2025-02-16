#include <iomanip>
#include <math.h>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <string>

#include "../../PDE_SWESphere2D/BenchmarksCombined.hpp"
#include "../../PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"
#include "../../PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_lg_irk.hpp"

#include "ceval.hpp"
#include "cencap.hpp"

extern "C"
{
#include <mpi.h>
}

bool timestep_check_output(Sphere2DDataCtxSDC *i_ctx,
						   int i_current_iter,
						   int i_niters)
{
	if (i_current_iter < i_niters)
	{
		// TODO: make this controllable via command line argument
		return false;
	}

	// get the simulation variables
	// sweet::Shacks::ShackDictionary *shackDict = i_ctx->get_simulation_variables();

	if (i_ctx->shackIOData->outputEachSimTime < 0)
	{
		// write no output between start and end of simulation
		return false;
	}

	if (i_ctx->shackIOData->outputEachSimTime == 0)
	{
		// write output at every time step
		return true;
	}

	if (i_ctx->shackTimestepControl->currentSimulationTime < i_ctx->shackIOData->_outputNextSimTime)
	{
		// we have not reached the next output time step
		return false;
	}

	if (i_ctx->shackTimestepControl->maxSimulationTime - i_ctx->shackTimestepControl->currentSimulationTime < 1e-3)
	{
		// do not write output if final time step is reached
		// (output will be written in cfinal anyways)
		return false;
	}

	// we have reached the next output time step
	return true;
}

/**
 * Write data to file and return string of file name
 */
std::string write_file(
	Sphere2DDataCtxSDC *i_ctx,
	const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
	const char *i_name //!< name of output variable
)
{
	char buffer[1024];

	// get the pointer to the Simulation Variables object
	// sweet::Shacks::ShackDictionary *shackDict = i_ctx.get_simulation_variables();

	// create copy
	sweet::Data::Sphere2D::DataSpectral sphere2DData(i_sphere2DData);

	// Write the data into the file
	const char *filename_template = i_ctx->shackIOData->outputFileName.c_str();
	sprintf(buffer,
			filename_template,
			i_name,
			i_ctx->shackTimestepControl->currentSimulationTime * i_ctx->shackIOData->outputFormatTimeScale);
	sphere2DData.file_write_binary_spectral(buffer);

	return buffer;
}

extern "C"
{
	// initialization of the variables (initial condition)
	void cinitial(
		Sphere2DDataCtxSDC *i_ctx,
		double i_t,
		double i_dt,
		Sphere2DDataVars *o_Y)
	{
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = o_Y->get_phi_pert();
		sweet::Data::Sphere2D::DataSpectral &vrt_Y = o_Y->get_vrt();
		sweet::Data::Sphere2D::DataSpectral &div_Y = o_Y->get_div();

		// get the sweet::Shacks::ShackDictionary object from context
		// sweet::Shacks::ShackDictionary *shackDict(i_ctx->get_simulation_variables());

		// if (shackDict->benchmark.use_topography)
		//	write_file(*i_ctx, shackDict->benchmark.h_topo, "prog_h_topo");

		PDE_SWESphere2D::Benchmarks::BenchmarksCombined *benchmarks = i_ctx->get_swe_benchmark();

		// use dealiased physical space for setup
		// get operator for this level
		sweet::Data::Sphere2D::Operators *op = i_ctx->get_sphere2d_operators();
		benchmarks->setup_1_registerAllBenchmark();
		benchmarks->setup_2_shackRegistration(i_ctx->shackDict);
		benchmarks->setup_3_benchmarkDetection();
		benchmarks->setup_4_benchmarkSetup_1_withoutOps();
		benchmarks->setup_5_benchmarkSetup_2_withOps(op);

		benchmarks->benchmark->getInitialState(phi_pert_Y, vrt_Y, div_Y);
		// output the configuration
		i_ctx->shackDict->printShackData();

		if (rank == 0)
		{
			if (i_ctx->shackIOData->outputEachSimTime >= 0)
			{
				write_file(i_ctx, phi_pert_Y, "prog_phi_pert");
				write_file(i_ctx, vrt_Y, "prog_vrt");
				write_file(i_ctx, div_Y, "prog_div");
			}
			if (i_ctx->shackIOData->outputEachSimTime < 0)
			{
				// do not write output
				i_ctx->shackIOData->_outputNextSimTime = i_ctx->shackTimestepControl->maxSimulationTime;
			}
			else if (i_ctx->shackIOData->outputEachSimTime > 0)
			{
				// write output every output_each_sim_seconds
				// next output time is thus equal to output_each_sim_seconds
				i_ctx->shackIOData->_outputNextSimTime = i_ctx->shackIOData->outputEachSimTime;
			}
			else
			{
				// output at every time step
				i_ctx->shackIOData->_outputNextSimTime = i_ctx->shackTimestepControl->currentTimestepSize;
			}
		}

		sweet::Data::Sphere2D::DataSpectral phi_pert_Y_init(phi_pert_Y);
		sweet::Data::Sphere2D::DataSpectral phi_pert_Y_final(phi_pert_Y);
		phi_pert_Y_init -= phi_pert_Y_final;

		sweet::Data::Sphere2D::DataSpectral div_Y_init(div_Y);
		sweet::Data::Sphere2D::DataSpectral div_Y_final(div_Y);
		div_Y_init -= div_Y_final;

		sweet::Data::Sphere2D::DataSpectral vrt_Y_init(vrt_Y);
		sweet::Data::Sphere2D::DataSpectral vrt_Y_final(vrt_Y);
		vrt_Y_init -= vrt_Y_final;
	}

	// finalizes the time step when libpfasst is done
	// currently does nothing else than outputting the solution
	void cfinal(
		Sphere2DDataCtxSDC *i_ctx,
		Sphere2DDataVars *i_Y,
		int i_nnodes,
		int i_niters)
	{
		int rank = 0;
		int nprocs = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		const sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = i_Y->get_phi_pert();
		const sweet::Data::Sphere2D::DataSpectral &vrt_Y = i_Y->get_vrt();
		const sweet::Data::Sphere2D::DataSpectral &div_Y = i_Y->get_div();

		// const int& level_id = i_Y->get_level();

		// get the sweet::Shacks::ShackDictionary object from context
		//sweet::Shacks::ShackDictionary *shackDict(i_ctx->get_simulation_variables());

		sweet::Data::Sphere2D::DataSpectral phi_pert_Y_init(phi_pert_Y);
		sweet::Data::Sphere2D::DataSpectral phi_pert_Y_final(phi_pert_Y);
		phi_pert_Y_init -= phi_pert_Y_final;

		sweet::Data::Sphere2D::DataSpectral div_Y_init(div_Y);
		sweet::Data::Sphere2D::DataSpectral div_Y_final(div_Y);
		div_Y_init -= div_Y_final;

		sweet::Data::Sphere2D::DataSpectral vrt_Y_init(vrt_Y);
		sweet::Data::Sphere2D::DataSpectral vrt_Y_final(vrt_Y);
		vrt_Y_init -= vrt_Y_final;

		if (i_ctx->shackIOData->outputEachSimTime < 0)
		{
			// do not write output
			return;
		}

		if (rank == 0)
		{
			std::string filename = "prog_phi_pert";
			write_file(i_ctx, phi_pert_Y, filename.c_str());

			filename = "prog_vrt";
			write_file(i_ctx, vrt_Y, filename.c_str());

			filename = "prog_div";
			write_file(i_ctx, div_Y, filename.c_str());
		}
	}

	// evaluates the explicit (nonlinear) piece
	void ceval_f1(Sphere2DDataVars *i_Y,
				  double i_t,
				  Sphere2DDataCtxSDC *i_ctx,
				  Sphere2DDataVars *o_F1)
	{
		const sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = i_Y->get_phi_pert();
		const sweet::Data::Sphere2D::DataSpectral &vrt_Y = i_Y->get_vrt();
		const sweet::Data::Sphere2D::DataSpectral &div_Y = i_Y->get_div();

		sweet::Data::Sphere2D::DataSpectral &phi_pert_F1 = o_F1->get_phi_pert();
		sweet::Data::Sphere2D::DataSpectral &vrt_F1 = o_F1->get_vrt();
		sweet::Data::Sphere2D::DataSpectral &div_F1 = o_F1->get_div();

		// get the time step parameters
		//sweet::Shacks::ShackDictionary *shackDict = i_ctx->get_simulation_variables();

		PDESWESphere2DTS_lg_erk_lc_n_erk *timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper();
		// compute the explicit nonlinear right-hand side
		timestepper->euler_timestep_update_lc_n(
			phi_pert_Y,
			vrt_Y,
			div_Y,
			phi_pert_F1,
			vrt_F1,
			div_F1,
			i_ctx->shackTimestepControl->currentSimulationTime);
	}

	// evaluates the implicit (linear) piece
	void ceval_f2(Sphere2DDataVars *i_Y,
				  double i_t,
				  Sphere2DDataCtxSDC *i_ctx,
				  Sphere2DDataVars *o_F2)
	{
		const sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = i_Y->get_phi_pert();
		const sweet::Data::Sphere2D::DataSpectral &vrt_Y = i_Y->get_vrt();
		const sweet::Data::Sphere2D::DataSpectral &div_Y = i_Y->get_div();

		sweet::Data::Sphere2D::DataSpectral &phi_pert_F2 = o_F2->get_phi_pert();
		sweet::Data::Sphere2D::DataSpectral &vrt_F2 = o_F2->get_vrt();
		sweet::Data::Sphere2D::DataSpectral &div_F2 = o_F2->get_div();

		// get the time step parameters
		//sweet::Shacks::ShackDictionary *shackDict = i_ctx->get_simulation_variables();

		PDESWESphere2DTS_lg_erk_lc_n_erk *timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper();
		// compute the linear right-hand side
		timestepper->euler_timestep_update_linear(
			phi_pert_Y,
			vrt_Y,
			div_Y,
			phi_pert_F2,
			vrt_F2,
			div_F2,
			i_ctx->shackTimestepControl->currentSimulationTime);
	}

	// solves the first implicit system for io_Y
	// then updates o_F2 with the new value of F2(io_Y)
	void ccomp_f2(
		Sphere2DDataVars *io_Y,
		double i_t,
		double i_dtq,
		Sphere2DDataVars *i_Rhs,
		Sphere2DDataCtxSDC *i_ctx,
		Sphere2DDataVars *o_F2)
	{
		// get the time step parameters
		//sweet::Shacks::ShackDictionary *shackDict = i_ctx->get_simulation_variables();

		sweet::Data::Sphere2D::DataSpectral &phi_pert_Y = io_Y->get_phi_pert();
		sweet::Data::Sphere2D::DataSpectral &vrt_Y = io_Y->get_vrt();
		sweet::Data::Sphere2D::DataSpectral &div_Y = io_Y->get_div();

		const sweet::Data::Sphere2D::DataSpectral &phi_pert_Rhs = i_Rhs->get_phi_pert();
		const sweet::Data::Sphere2D::DataSpectral &vrt_Rhs = i_Rhs->get_vrt();
		const sweet::Data::Sphere2D::DataSpectral &div_Rhs = i_Rhs->get_div();

		// first copy the rhs into the solution vector
		// this is needed to call the SWEET function run_timestep
		phi_pert_Y = phi_pert_Rhs;
		vrt_Y = vrt_Rhs;
		div_Y = div_Rhs;

		if (i_dtq == 0)
		{
			// quadrature weight is zero -> return trivial solution
			// y = rhs (already done), f = 0.0
			c_sweet_data_setval(o_F2, 0.0);
			return;
		}

		PDESWESphere2DTS_lg_irk *timestepper = i_ctx->get_lg_irk_timestepper();
		// solve the implicit system using the Helmholtz solver
		timestepper->runTimestep(
			phi_pert_Y,
			vrt_Y,
			div_Y,
			i_dtq,
			i_ctx->shackTimestepControl->maxSimulationTime);

		sweet::Data::Sphere2D::DataSpectral &phi_pert_F2 = o_F2->get_phi_pert();
		sweet::Data::Sphere2D::DataSpectral &vrt_F2 = o_F2->get_vrt();
		sweet::Data::Sphere2D::DataSpectral &div_F2 = o_F2->get_div();

		phi_pert_F2 = (phi_pert_Y - phi_pert_Rhs) / i_dtq;
		vrt_F2 = (vrt_Y - vrt_Rhs) / i_dtq;
		div_F2 = (div_Y - div_Rhs) / i_dtq;

		return;
	}
}
