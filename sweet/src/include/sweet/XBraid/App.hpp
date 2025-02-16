/*
 * App.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef INCLUDE_SWEET_XBRAID_APP_HPP
#define INCLUDE_SWEET_XBRAID_APP_HPP

#include <xbraid/braid.hpp>

#if SWEET_GUI
#include<sweet/GUI/VisSweet.hpp>
#endif

#include <algorithm>

#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/XBraid/Shack.hpp>

#include <sweet/Data/GenericContainer/ConfigBase.hpp>

#include "Vector.hpp"
#include "TimeTree.hpp"
#include "FileOutput.hpp"
#include "Error.hpp"

namespace sweet {
namespace XBraid {

// Wrapper for BRAID's App object·
// --> Put all time INDEPENDENT information here
class App
			: public BraidApp
{

public:

	sweet::Error::Base error;

	// General Shacks
	sweet::Shacks::Dictionary*			shackDict;
	sweet::XBraid::Shack*				shackXBraid;
	sweet::Parallelization::Shack*			shackParallelization;
	sweet::TimeTree::Shack*				shackTimestepControl;
	sweet::IO::Shack*				shackIOData;

	// Level-defined shacks
	std::vector<sweet::Shacks::Dictionary*>		shacksDict_levels;
	std::vector<sweet::TimeTree::Shack*>		shacksTimestepControl_levels;

	// Time step size on the finest level
	double dt;

	// Time stepper on each level
	std::vector<std::string> tsms;
	std::vector<int> tsos;
	std::vector<int> tsos2;
	std::vector<bool> is_SL;
	bool contains_SL = false;
	std::vector<sweet::XBraid::TimeTree*>	timeSteppers_NewTS;

	// Initial solution and initial guess
	sweet::XBraid::Vector*		prog_initial_solution = nullptr;
	sweet::XBraid::Vector*		prog_initial_guess = nullptr;

	// Solution on each level
	std::vector<sweet::XBraid::Vector*>	data_container;

	// Viscosity on each level
	std::vector<int>	viscosity_orders;
	std::vector<double>	viscosity_coefficients;


#if SWEET_GUI
	std::vector<SimulationGuiCallBacks*> levels_simulations;
#endif

	// Overestimated buffer size for MPI communications
	int size_buffer;		// overestimated

	// MPI rank
	int rank;

	// Reference solutions used for online error computation
	// REF: given reference solution
	// FINE: solution on the finest level
	std::vector<sweet::XBraid::Vector*> xbraid_data_ref_exact;
	std::vector<sweet::XBraid::Vector*> xbraid_data_fine_exact;

	// Solution from previous timestep (for SL)
	std::vector<std::vector<sweet::XBraid::Vector*>> 	sol_prev; // sol_prev[level][timestep]
	std::vector<std::vector<int>>				sol_prev_iter; // store iteration in which the solution has been stored
	std::vector<int>					first_timeid_level; // store first timestep (time id) in this level and processor
	std::vector<int>					last_timeid_level; // store last timestep (time id) in this level and processor

	// Custom time grid
	std::vector<double> custom_time_steps = {};

	// Effective number of levels ( <= xbraid_max_levels)
	int nlevels = -1;

	// Use new timesteppers (...,...,...)
	bool useNewTimeSteppers = true;

	// Output solution and errors
	sweet::XBraid::FileOutput*	file_output = nullptr;
	sweet::XBraid::Error*		error_xbraid = nullptr;

public:

	// Constructor·
	App(		MPI_Comm		i_comm_t,
			int			i_rank,
			double			i_tstart,
			double			i_tstop,
			int 			i_ntime
	)
		:
			BraidApp(i_comm_t, i_tstart, i_tstop, i_ntime),
			rank(i_rank)
	{
	}


	virtual ~App()
	{

		for (std::vector<sweet::XBraid::TimeTree*>::iterator it = timeSteppers_NewTS.begin();
							it != timeSteppers_NewTS.end();
							it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<sweet::XBraid::Vector*>::iterator it = data_container.begin();
							it != data_container.end();
							it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<sweet::XBraid::Vector*>::iterator it = xbraid_data_ref_exact.begin();
								it != xbraid_data_ref_exact.end();
								it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<sweet::XBraid::Vector*>::iterator it = xbraid_data_fine_exact.begin();
								it != xbraid_data_fine_exact.end();
								it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<std::vector<sweet::XBraid::Vector*>>::iterator it = sol_prev.begin();
										it != sol_prev.end();
										it++)
			for (std::vector<sweet::XBraid::Vector*>::iterator it2 = it->begin();
										it2 != it->end();
										it2++)
				if (*it2)
				{
					delete *it2;
					*it2 = nullptr;
				}

#if SWEET_GUI
		for (std::vector<sweet::SimulationGuiCallbacks*>::iterator it =
				levels_simulations.begin();
				it != levels_simulations.end();
			it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}
#endif

		if (prog_initial_solution)
		{
			delete prog_initial_solution;
			prog_initial_solution = nullptr;
		}

		if (prog_initial_guess)
		{
			delete prog_initial_guess;
			prog_initial_guess = nullptr;
		}

	}

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		std::cout << "GLOBAL REGISTRATION" << std::endl;
		shackDict = &io_shackDict;
		shackIOData = io_shackDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimestepControl = io_shackDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackXBraid = io_shackDict.getAutoRegistration<sweet::XBraid::Shack>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(io_shackDict);

		return true;
	}

public:
	bool shackRegistrationLevels(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{

		std::cout << "LOCAL REGISTRATION" << std::endl;
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			std::cout << " --> REGISTRATION LEVEL " << level << std::endl;
			sweet::Shacks::Dictionary* shackDict_level = new sweet::Shacks::Dictionary;
			*shackDict_level = *io_shackDict;

			sweet::TimeTree::Shack* shackTimestepControl_level = nullptr;

			shackTimestepControl_level = shackDict_level->getAutoRegistration<sweet::TimeTree::Shack>();

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackDict_level);

			shacksDict_levels.push_back(shackDict_level);
			shacksTimestepControl_levels.push_back(shackTimestepControl_level);
		}

		return true;
	}


public:
	void setup(
			BraidCore& i_core,
			sweet::XBraid::Vector& i_prog_initial_solution,
			sweet::XBraid::Vector& i_prog_initial_guess,
			sweet::XBraid::FileOutput* i_file_output,
			sweet::XBraid::Error* i_error_xbraid
	)
	{

		/////////////////////////////////////////////////
		// get parameters from simVars and set to Core //
		/////////////////////////////////////////////////

		i_core.SetMaxLevels(shackXBraid->xbraid_max_levels);
		/////i_core.SetIncrMaxLevels();

		i_core.SetSkip(shackXBraid->xbraid_skip);

		i_core.SetMinCoarse(shackXBraid->xbraid_min_coarse);

		///i_core.SetRelaxOnlyCG(shackXBraid->xbraid_relax_only_cg);

		i_core.SetNRelax(-1, shackXBraid->xbraid_nrelax);
		if (shackXBraid->xbraid_nrelax0 > -1)
			i_core.SetNRelax(0, shackXBraid->xbraid_nrelax0);

		i_core.SetAbsTol(shackXBraid->xbraid_tol);
		i_core.SetRelTol(shackXBraid->xbraid_tol);

		i_core.SetTemporalNorm(shackXBraid->xbraid_tnorm);

		i_core.SetCFactor(-1, shackXBraid->xbraid_cfactor);
		if (shackXBraid->xbraid_cfactor0 > -1)
			i_core.SetCFactor(0, shackXBraid->xbraid_cfactor0);

		///i_core.SetPeriodic(shackXBraid->xbraid_periodic);

		////i_core.SetResidual();

		i_core.SetMaxIter(shackXBraid->xbraid_max_iter);

		i_core.SetPrintLevel(shackXBraid->xbraid_print_level);

		i_core.SetSeqSoln(shackXBraid->xbraid_use_seq_soln);

		i_core.SetAccessLevel(shackXBraid->xbraid_access_level);

		i_core.SetNFMG(shackXBraid->xbraid_fmg);
		if (shackXBraid->xbraid_fmg)
			i_core.SetFMG();

		i_core.SetNFMGVcyc(shackXBraid->xbraid_fmg_vcyc);

		i_core.SetStorage(shackXBraid->xbraid_storage);

		//i_core.SetRevertedRanks(shackXBraid->xbraid_reverted_ranks);

		////i_core.SetRefine(shackXBraid->xbraid_refine);
		///i_core.SetMaxRefinements(shackXBraid->xbraid_max_Refinements);

		i_core.SetTimeGrid(App::sweet_TimeGrid);

		///setup_timesteppers();
		setup(i_prog_initial_solution, i_prog_initial_guess, i_file_output, i_error_xbraid);
	}

public:
	/*
	 * Setup timesteppers for each level.
	 * IMPORTANT: this function must be called after setting up initial conditions, since the benchmark may modify simulation parameters
	 *            call it in the first call of Step function.
	 */
	void setup_timesteppers()
	{

		////////////////////////////////////
		// get tsm and tso for each level //
		////////////////////////////////////
		tsms = getLevelParameterFromParameters<std::string>("timestepping_method");

		// create a timeSteppers instance for each level
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{

			// Configure timesteppers with the correct timestep for this level
			shacksTimestepControl_levels[level]->printShack();
			shacksTimestepControl_levels[level]->currentTimestepSize = shacksTimestepControl_levels[level]->currentTimestepSize *
												std::pow(shackXBraid->xbraid_cfactor, level);
			shacksTimestepControl_levels[level]->printShack();

			if (rank == 0)
				std::cout << "Timestep size at level " << level << " : " << shacksTimestepControl_levels[level]->currentTimestepSize << std::endl;

			sweet::XBraid::TimeTree* tsm = createConfigureNewTimestepper(level);
			this->timeSteppers_NewTS.push_back(tsm);

			// check if tsm is SL
			/////if ( std::find(SL_tsm.begin(), SL_tsm.end(), tsms[level]) == SL_tsm.end())
			if (tsms[level].find("SL") == std::string::npos && tsms[level].find("SETTLS") == std::string::npos)
				is_SL.push_back(false);
			else
			{
				is_SL.push_back(true);
				contains_SL = true; // requires extra communication
			}
		}

	}

public:
	virtual
	bool setupDataConfigOps() = 0;

public:
	void setup(
			sweet::XBraid::Vector& i_prog_initial_solution,
			sweet::XBraid::Vector& i_prog_initial_guess,
			sweet::XBraid::FileOutput* i_file_output,
			sweet::XBraid::Error* i_error_xbraid
	)
	{


		file_output = i_file_output;
		error_xbraid = i_error_xbraid;

		// create vectors for storing solutions from previous timestep (SL)
		for (int i = 0; i < shackXBraid->xbraid_max_levels; i++)
		{
			std::vector<sweet::XBraid::Vector*> v = {};
			sol_prev.push_back(v);

			std::vector<int> w = {};
			sol_prev_iter.push_back(w);

			first_timeid_level.push_back(INT_MAX);
			last_timeid_level.push_back(-1);
		}

		// set custom time grid
		double t = 0;
		while (t < shackTimestepControl->maxSimulationTime - 1e-10)
		{
			double dt = shackTimestepControl->currentTimestepSize;
			double dt2;
			if ( t + dt < shackTimestepControl->maxSimulationTime - 1e-10)
				dt2 = dt;
			else
				dt2 = shackTimestepControl->maxSimulationTime - t;
			custom_time_steps.push_back(dt2);
			t += dt2;
		}

		setupDataConfigOps();

		prog_initial_solution = createNewVector(0);
		prog_initial_guess = createNewVector(0);
		prog_initial_solution->op_setVector(i_prog_initial_solution);
		prog_initial_guess->op_setVector(i_prog_initial_guess);

	}


	/*
	 * Get specific parameters for each of the N levels
	 * Input string must contain 1, 2 or N orders separated by comma:
	 *  - 1 tso: same tso for all levels;
	 *  - 2 tso: first one is for level 0 (finest); all coarse levels use the second one
	 *  - N tso: i-th level uses the i-th tso.
	 */
	template<typename T>
	std::vector<T> getLevelParameterFromParameters(std::string i_param_name, int i_order = 0)
	{

		if ( ! (
			i_param_name == "timestepping_method" ||
			i_param_name == "timestepping_order" ||
			i_param_name == "viscosity_order" ||
			i_param_name == "viscosity_coefficient"
			))
			SWEETErrorFatal("Wrong param_name " + i_param_name);

		if (i_param_name == "timestepping_order")
			SWEET_ASSERT(i_order == 1 || i_order == 2);

		std::vector<T> out = {};
		std::stringstream all_param;

		if (i_param_name == "timestepping_method")
		{
			all_param = std::stringstream(this->shackXBraid->xbraid_timestepping_method);
			if (this->shackXBraid->xbraid_timestepping_method.find('(') != std::string::npos)
				this->useNewTimeSteppers = true;
			else
				this->useNewTimeSteppers = false;
		}
		else if (i_param_name == "timestepping_order")
		{
			if (i_order == 1)
				all_param = std::stringstream(shackXBraid->xbraid_timestepping_order);
			else if (i_order == 2)
				all_param = std::stringstream(shackXBraid->xbraid_timestepping_order2);
		}
		else if (i_param_name == "viscosity_order")
			all_param = std::stringstream(shackXBraid->xbraid_viscosity_order);
		else if (i_param_name == "viscosity_coefficient")
			all_param = std::stringstream(shackXBraid->xbraid_viscosity_coefficient);

		if ( ( i_param_name != "timestepping_method" ) || !this->useNewTimeSteppers  )
		{
			/////SWEETErrorFatal("DEPRECATED");
			while (all_param.good())
			{
				std::string str;
				getline(all_param, str, ',');
				std::stringstream ss(str);
				T conv;
				if (ss >> conv)
					out.push_back(conv);
				else
					SWEETErrorFatal("Unable to convert parameter: " + str);
			}
		}
		else
		{
			// string should be: "(...,...,...);(...,...,...);..."
			while (all_param.good())
			{
				std::string str;
				getline(all_param, str, ';');
				std::stringstream ss(str);
				T conv;
				if (ss >> conv)
					out.push_back(conv);
				else
					SWEETErrorFatal("Unable to convert parameter: " + str);
			}
		}

		if ( ! (out.size() == 1 || out.size() == 2 || (int)out.size() == shackXBraid->xbraid_max_levels ) )
			SWEETErrorFatal("xbraid_" + i_param_name +  "must contain 1, 2 or N timestepping orders.");

		// all levels use same param
		if (out.size() == 1)
			for (int level = 1; level < shackXBraid->xbraid_max_levels; level++)
				out.push_back(out[0]);

		// all coarse levels use same tso
		if (out.size() == 2)
			for (int level = 2; level < shackXBraid->xbraid_max_levels; level++)
				out.push_back(out[1]);

		return out;
	}

public:
	virtual
	sweet::XBraid::TimeTree* createConfigureNewTimestepper(int i_level) = 0;

public:
	virtual
	sweet::XBraid::Vector* createNewVector(int i_level) = 0;

private:
	void store_prev_solution(
					sweet::XBraid::Vector* i_U,
					int i_time_id,
					int i_level,
					int iter
				)
	{
		// if not SL scheme: nothing to do
		//if ( std::find(SL_tsm.begin(), SL_tsm.end(), tsms[i_level]) == SL_tsm.end())
		if ( ! is_SL[i_level] )
			return;

		// if solution has already been stored in this iteration: nothing to do
		if ( sol_prev_iter[i_level][i_time_id] == iter )
			return;
		///SWEET_ASSERT(sol_prev_iter[i_level][i_time_id] == iter - 1);

		// create vector if necessary
		if ( ! sol_prev[i_level][i_time_id] )
			sol_prev[i_level][i_time_id] = createNewVector(i_level);

		// set solution
		sol_prev[i_level][i_time_id]->op_setVector(*i_U);
		sol_prev_iter[i_level][i_time_id] = iter;
		first_timeid_level[i_level] = std::min(first_timeid_level[i_level], i_time_id);
		last_timeid_level[i_level] = std::max(last_timeid_level[i_level], i_time_id);
	}

	void set_prev_solution(
					sweet::XBraid::Vector* i_U,
					int i_time_id,
					int i_level
				)
	{

		/////std::cout << "WARNING: not yet implemented" << std::endl;
		/////return;

		// if not SL scheme: nothing to do
		///if (  std::find(SL_tsm.begin(), SL_tsm.end(), tsms[i_level]) == SL_tsm.end())
		if ( ! is_SL[i_level] )
			return;

		// if t == 0 or prev solution not available
		// then: prev_solution = solution
		bool prev_sol_exists = true;

		if ( i_time_id == 0 )
			prev_sol_exists = false;
		if ( prev_sol_exists && (!sol_prev[i_level][i_time_id - 1]) )
			prev_sol_exists = false;

		// only store prev solution if it is not the first time step inside a coarse slice
		//if (i_time_id % shackXBraid->xbraid_cfactor == 0)
		if (i_level < nlevels - 1)
			prev_sol_exists = false;

		//prev_sol_exists = false;

		if (prev_sol_exists)
		{
			this->timeSteppers_NewTS[i_level]->storePreviousSolution(sol_prev[i_level][i_time_id - 1]);
		}
		else
		{
			this->timeSteppers_NewTS[i_level]->storePreviousSolution(i_U);
		}
	}

public:
	/* --------------------------------------------------------------------
	 * Time integrator routine that performs the update
	 *   u_i = Phi_i(u_{i-1}) + g_i 
	 * 
	 * When Phi is called, u is u_{i-1}.
	 * The return value is that u is set to u_i upon completion
	 *
	 * Always receives and returns a solution defined on the finest spatial grid
	 * Spatial interpolation is performed if necessary
	 * -------------------------------------------------------------------- */
	braid_Int
	Step(
			braid_Vector		io_U,
			braid_Vector		i_ustop,
			braid_Vector		i_fstop,
			BraidStepStatus&	io_status
			)
	{

		// First call of this function: create timesteppers
		// Benchmark::setup_initial_conditions has already been called
		if (timeSteppers_NewTS.size() == 0)
		{

			// if rank > 0, call setup_init_conditions to ensure that parameters from benchmark are set
			if (rank > 0)
			{
				braid_Vector dummy;
				Init(0., &dummy);
			}

			setup_timesteppers();

			io_status.GetNLevels(&nlevels);
		}


		// Vector defined in the finest level
		sweet::XBraid::Vector* U = (sweet::XBraid::Vector*) io_U;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;
		int time_id;
		int iter;

		/* Grab status of current time step */
		io_status.GetTstartTstop(&tstart, &tstop);
		io_status.GetLevel(&level);
		io_status.GetTIndex(&time_id);
		io_status.GetIter(&iter);

		// Vector defined in the current level (defined via interpolation if necessary)
		sweet::XBraid::Vector* U_level = createNewVector(level);

		// Interpolate to coarser grid in space if necessary
		if (shackXBraid->xbraid_spatial_coarsening && level > 0)
			U_level->restrict(U);
		else
			U_level->op_setVector(U);


		// create containers for prev solution
		if (sol_prev[level].size() == 0)
		{
			// store nt solutions (overestimated for coarse levels!)
			int nt;
			io_status.GetNTPoints(&nt);
			for (int i = 0; i < nt + 1; i++)
			{
				sol_prev[level].push_back(nullptr);
				sol_prev_iter[level].push_back(-1);
			}
				///sol_prev[level].push_back(createNewVector());
		}

		// store solution for SL
		if (time_id == 0)
			store_prev_solution(U_level, time_id, level, iter);

		// set prev solution for SL
		set_prev_solution(U_level, time_id, level);

		// TODO: check if this is thread safe
		/////simVars->timecontrol.current_simulation_time = tstart;
		/////simVars->timecontrol.current_timestepSize = tstop - tstart;
		// TODO

		//////std::cout << rank << " " << iter << " " << level << " " << tstart << " " << tstop << std::endl;

		sweet::XBraid::Vector* U_level_copy = createNewVector(level);
		U_level_copy->op_setVector(*U_level);

		timeSteppers_NewTS[level]->runIntegration(
					U_level_copy,
					U_level,
					tstart
		);

		delete U_level_copy;


		store_prev_solution(U_level, time_id + 1, level, iter);

		// Interpolate to finest grid in space if necessary
		if (this->useNewTimeSteppers)
		{
			if (this->shackXBraid->xbraid_spatial_coarsening && level > 0)
				U->pad_zeros(U_level);
			else
				U->op_setVector(U_level);
		}
		else
			SWEETErrorFatal("DEPRECATED");

		/* Tell XBraid no refinement */
		io_status.SetRFactor(1);

		delete U_level;

		return 0;
	}

		/* --------------------------------------------------------------------
		 * -------------------------------------------------------------------- */

	virtual braid_Int
	Residual(
				braid_Vector			i_ustop,
				braid_Vector			o_r,
				BraidStepStatus&		io_status
		)
	{

		//braid_StepStatus& status = (braid_StepStatus&) io_status;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;

		/* Grab status of current time step */
		io_status.GetTstartTstop(&tstart, &tstop);
	
		/* Grab level */
		io_status.GetLevel(&level);
	
		/* Set the new dt in the user's manager*/
		dt = tstop - tstart;


		///sweet::XBraid::Vector* u = (sweet::XBraid::Vector*) i_ustop;
		///sweet::XBraid::Vector* r = (sweet::XBraid::Vector*) o_r;

		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Create a vector object for a given time point.
	 * This function is only called on the finest level.
	 * -------------------------------------------------------------------- */
	braid_Int
	Init(
			double		i_t,
			braid_Vector*	o_U
			)
	{

		sweet::XBraid::Vector* U = createNewVector(0);

		if (i_t == tstart)
			U->op_setVector(prog_initial_solution);
		else
			U->op_setVector(prog_initial_guess);

		*o_U = (braid_Vector) U;

		return 0;
	}

	/* --------------------------------------------------------------------
	 * Create a copy of a vector object.
	 * -------------------------------------------------------------------- */
	virtual
	braid_Int
	Clone(
			braid_Vector	i_U,
			braid_Vector*	o_V
		) = 0;

	/* --------------------------------------------------------------------
	 * Destroy vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Free(
			braid_Vector	i_U)
	{
		sweet::XBraid::Vector* U = (sweet::XBraid::Vector*) i_U;
		delete U;

		return 0;
	}

	/* --------------------------------------------------------------------
	 * Compute vector sum y = alpha*x + beta*y.
	 * -------------------------------------------------------------------- */
	braid_Int
	Sum(
			double			i_alpha,
			braid_Vector		i_X,
			double			i_beta,
			braid_Vector		io_Y
		)
	{

		// copy i_X
		braid_Vector i_X_copy;
		Clone(i_X, &i_X_copy);
		sweet::XBraid::Vector* X = (sweet::XBraid::Vector*) i_X_copy;

		////sweet::XBraid::Vector* X = (sweet::XBraid::Vector*) i_X;
		sweet::XBraid::Vector* Y = (sweet::XBraid::Vector*) io_Y;

		///*Y = *X * i_alpha + *Y * i_beta;

		X->op_mulScalar(i_alpha);
		Y->op_mulScalar(i_beta);
		Y->op_addVector(*X);

		// Delete copy
		Free(i_X_copy);

		return 0;
	}


	/* --------------------------------------------------------------------
	 * User access routine to spatial solution vectors and allows for user
	 * output.  The default XBraid parameter of access_level=1, calls 
	 * my_Access only after convergence and at every time point.
	 * -------------------------------------------------------------------- */
	braid_Int
	Access(
				braid_Vector		i_U,
				BraidAccessStatus&	io_astatus
			)
	{
		/////double     tstart         = (tstart);
		/////double     tstop          = (tstop);
		////int        nt             = (nt);

		///double     rnorm, disc_err, t;
		///int        iter, level, done, index, myid, it;
		///char       filename[255], filename_mesh[255], filename_err[255], filename_sol[255];
		double     t;
		int        iter, level, it;


		sweet::XBraid::Vector *U = (sweet::XBraid::Vector*) i_U;

		/* Retrieve current time from Status Object */
		////braid_AccessStatusGetT(astatus, &t);
		io_astatus.GetT(&t);
		io_astatus.GetTIndex(&it);
		io_astatus.GetIter(&iter);
		io_astatus.GetLevel(&level);

		/* Retrieve XBraid State Information from Status Object */
		///////////MPI_Comm_rank(app->comm_x, &myid);
		///////////braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
		///////////braid_AccessStatusGetResidual(astatus, &rnorm);

		////std::cout << "ACCESS " << rank << " " << level << " " << t << std::endl;

		if (level == 0 /*&& rank == 0*/)
		{


			// Decide whether to output or not
			bool do_output = false;
			double small = 1e-10;

			// output each time step if:
			// output_timestep < 0 (i.e. output every timestep)
			// t == 0
			// t == Tmax
			// t is a multiple of dt_output
			if (
				shackIOData->outputEachSimTime < 0 ||
				std::abs(t) < small ||
				std::abs(t - shackTimestepControl->maxSimulationTime) < small ||
				fmod(t, shackIOData->outputEachSimTime) < small ||
				std::abs(fmod(t, shackIOData->outputEachSimTime) - shackIOData->outputEachSimTime) < small
			)
				do_output = true;

			if (shackXBraid->xbraid_no_output)
				do_output = false;

			if (do_output)
			{

				// Load reference solution
				if (iter == 0)
				{

					if (shackXBraid->xbraid_load_ref_csv_files)
					{

						if (xbraid_data_ref_exact.size() == 0)
						{
							int nt;
							io_astatus.GetNTPoints(&nt);
							for (int i = 0; i < nt + 1; i++)
								xbraid_data_ref_exact.push_back(createNewVector(0));
						}

						error_xbraid->loadRefSolution(
										this->xbraid_data_ref_exact[it],
										t,
										this->shackXBraid->xbraid_path_ref_csv_files
						);
					}

					if (shackXBraid->xbraid_load_fine_csv_files)
					{

						if(xbraid_data_fine_exact.size() == 0)
						{
							int nt;
							io_astatus.GetNTPoints(&nt);
							for (int i = 0; i < nt + 1; i++)
								xbraid_data_fine_exact.push_back(createNewVector(0));
						}

						error_xbraid->loadRefSolution(
										this->xbraid_data_fine_exact[it],
										t,
										this->shackXBraid->xbraid_path_fine_csv_files
						);
					}

				}

				// Output physical solution to file
				if (shackXBraid->xbraid_store_iterations)
				{
					file_output->output_data_file(
								U,
								iter,
								it,
								t
					);
				}

				// Compute and store errors w.r.t. ref solution
				if (shackXBraid->xbraid_load_ref_csv_files)
				{

					if (it >= 0)
					{

						sweet::XBraid::Vector* diff = createNewVector(0);
						diff->op_setVector(*U);
						diff->op_subVector(*xbraid_data_ref_exact[it]);

						error_xbraid->computeStoreError(
										diff,
										this->xbraid_data_ref_exact[it],
										iter /* + 1 */,
										it,
										t,
										this->shackXBraid->xbraid_path_ref_csv_files,
										"ref"
						);

						diff->clear();
					}
				}
				// Compute and store errors w.r.t. fine (serial) solution
				if (shackXBraid->xbraid_load_fine_csv_files)
				{

					if (it >= 0)
					{

						sweet::XBraid::Vector* diff = createNewVector(0);
						diff->op_setVector(*U);
						diff->op_subVector(*xbraid_data_fine_exact[it]);

						error_xbraid->computeStoreError(
										diff,
										this->xbraid_data_fine_exact[it],
										iter /* + 1 */,
										it,
										t,
										this->shackXBraid->xbraid_path_fine_csv_files,
										"fine"
						);

						diff->clear();
					}
				}
			}

			// Store residual (residual per iteration)
			if (it == 0) {
				double res;
				io_astatus.GetResidual(&res);
				file_output->output_residual_file(
									res,
									iter
				);
			}

		}

		// TODO: verify if convergence stagnates and stop simulation

		return 0;
	}

	virtual
	int getMaxLevel() = 0;

	/* --------------------------------------------------------------------
	 * Compute norm of a spatial vector 
	 * -------------------------------------------------------------------- */
	braid_Int
	SpatialNorm(
			braid_Vector  i_U,
			double*       o_norm)
	{
		sweet::XBraid::Vector* U = (sweet::XBraid::Vector*) i_U;

		// Compute residual on the coarsest level
		// Restrict solution in spectral space then compute residual in physical space

		//int max_level = (int)this->config_levels.size() - 1;
		int max_level = getMaxLevel();
		sweet::XBraid::Vector* U_level = this->createNewVector(max_level);
		U_level->restrict(U);
		*o_norm = U_level->reduceNormLinfGrid();
		delete U_level;

		return 0;
	}

	/* --------------------------------------------------------------------
	 * Return buffer size needed to pack one spatial braid_Vector.  Here the
	 * vector contains one double at every grid point and thus, the buffer 
	 * size is the number of grid points.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufSize(
			int*			o_size,
			BraidBufferStatus&	o_status)
	{
		*o_size = size_buffer;
		return 0;
	}

	virtual
	int getBufferSize() = 0;

	/* --------------------------------------------------------------------
	 * Pack a braid_Vector into a buffer.
	 *
	 * Issue concerning SL: the first time step in each processor needs to receive
	 *                      the penult time step from the previous processor;
	 *                      However, communication is made only in level 0
	 * Solution (possibly not optimal): if SL is used in at least one level,
	 *                                  each communication includes the previous
	 *                                  time step of all levels.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufPack(
			braid_Vector		i_U,
			void*			o_buffer,
			BraidBufferStatus&	o_status
		)
	{

		sweet::XBraid::Vector* U = (sweet::XBraid::Vector*) i_U;

		std::complex<double>* dbuffer = (std::complex<double>*) o_buffer;

		// get buffer size
		int actual_size_buffer = getBufferSize();

		U->serialize(dbuffer);

		o_status.SetSize( actual_size_buffer );
		return 0;
	}

	/* --------------------------------------------------------------------
	 * Unpack a buffer and place into a braid_Vector
	 * -------------------------------------------------------------------- */
	braid_Int
	BufUnpack(
			void*			i_buffer,
			braid_Vector*		o_U,
			BraidBufferStatus&	io_status
		)
	{

		int level = 0;

		sweet::XBraid::Vector* U = createNewVector(level);

		std::complex<double>* dbuffer = (std::complex<double>*) i_buffer;

		U->deserialize(dbuffer);

		*o_U = (braid_Vector) U;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Define time grid
	 * -------------------------------------------------------------------- */
	static braid_Int
	sweet_TimeGrid(
				_braid_App_struct* i_app,
				braid_Real* i_ta,
				braid_Int* i_ilower,
				braid_Int* i_iupper
			)
	{

		App* app =  (App*) i_app;

		double tstart;
		int lower = *i_ilower;
		int upper = *i_iupper;

		/* Start from the global tstart to compute the local tstart */
		tstart = app->tstart;
		for (int i = 0; i < lower; i++)
			tstart += app->custom_time_steps[i];

		/* Assign time point values for local time point index values lower:upper */
		for (int i = lower; i <= upper; i++)
		{
			i_ta[i - lower] = tstart;
			tstart += app->custom_time_steps[i];
		}

		return 0;
	}

};

}}

#endif
