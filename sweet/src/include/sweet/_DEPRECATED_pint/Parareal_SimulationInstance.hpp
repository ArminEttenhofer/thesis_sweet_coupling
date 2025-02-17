/*
 * PararealSimulation.hpp
 *
 *  Created on: 25 Feb 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: programs/swe_sphere2d_benchmarks/BenchmarksSphere2DSWE.cpp
 * MULE_COMPILE_FILES_AND_DIRS: programs/swe_sphere2d_benchmarks/
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_SIMULATIONINSTANCE_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_SIMULATIONINSTANCE_HPP


#include <sweet/core/SimulationVariables.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
#include <sweet/_DEPRECATED_pint/PInT_Common.hpp>


#if SWEET_PARAREAL_SCALAR
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Scalar.hpp>

#elif SWEET_PARAREAL_CART2D
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Cart2DData_Spectral.hpp>
#include <sweet/Data/Cart2D/Cart2DData_Spectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Data/Cart2D/GridMapping.hpp>
#include "../../programs/swe_cart2d_timeintegrators/SWE_Cart2D_TimeSteppers.hpp"
#include "../../programs/burgers_timeintegrators/Burgers_Cart2D_TimeSteppers.hpp"
#include "../../programs/swe_cart2d_benchmarks/SWECart2DBenchmarksCombined.hpp"

#elif SWEET_PARAREAL_SPHERE2D
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Sphere2DData_Spectral.hpp>
#include <sweet/Data/Sphere2D/Sphere2DData_Spectral.hpp>
#include <sweet/core/sphere2D/Sphere2DOperators.hpp>
#include "../../programs/swe_sphere2d_timeintegrators/SWE_Sphere2D_TimeSteppers.hpp"
#include "../../programs/swe_sphere2d_benchmarks/BenchmarksSphere2DSWE.hpp"
#endif
//#include "../../programs/swe_sphere2d_benchmarks/BenchmarksSphere2DSWE.cpp"


namespace sweet {
namespace DEPRECATED_pint {

/**
 * Interface descriptions which are required
 * to run Parareal Simulations in SWEET.
 *
 * Note that these interface descriptions have
 * to be implemented by the simulation!
 *
 * These interfaces were ported from the Python implementation.
 */
template <class t_tsmType, int N>
class Parareal_SimulationInstance
			: public PInT_Common
{
public:

	// Simulation variables coarse level
	SimulationVariables* simVars_coarse = nullptr; // coarse level (dt, tsm, tso provided by specific parareal parameters)

	// Time slice
	double timeframe_start;
	double timeframe_end;
	int nb_timesteps_fine;
	int nb_timesteps_coarse;
	double dt_fine;
	double dt_coarse;

	// Data containers
	Parareal_GenericData* parareal_data_start = nullptr;
	Parareal_GenericData* parareal_data_fine = nullptr;
	Parareal_GenericData* parareal_data_coarse = nullptr; // defined on fine spatial mesh; used almost everywhere
	Parareal_GenericData* parareal_data_coarse_coarse_mesh = nullptr; // defined on coarse spatial mesh; used only on run_timestep_coarse
	Parareal_GenericData* parareal_data_output = nullptr;
	Parareal_GenericData* parareal_data_error = nullptr;
	Parareal_GenericData* parareal_data_coarse_previous_timestep = nullptr;
	Parareal_GenericData* parareal_data_coarse_previous_time_slice = nullptr;
	Parareal_GenericData* parareal_data_fine_previous_timestep = nullptr;
	Parareal_GenericData* parareal_data_fine_previous_time_slice = nullptr;
	Parareal_GenericData* parareal_data_ref_exact = nullptr;
	Parareal_GenericData* parareal_data_fine_exact = nullptr;
#if SWEET_DEBUG
	Parareal_GenericData* parareal_data_fine_exact_debug = nullptr;
#endif

	Parareal_GenericData* parareal_data_debug = nullptr;
	bool debug_contains_data = false;

	// Fine and coarse timesteppers
	t_tsmType* timeSteppersFine = nullptr;
	t_tsmType* timeSteppersCoarse = nullptr;

	////! list of SL schemes
	///std::vector<std::string> SL_tsm = {};


	bool compute_normal_modes = false;


	// For Burgers
	class BenchmarkErrors
	{
	public:
		// Max difference to initial conditions
		double benchmark_diff_u;
		double benchmark_diff_v;

		// Error measures L2 norm
		double benchmark_analytical_error_rms_u;
		double benchmark_analytical_error_rms_v;

		// Error measures max norm
		double benchmark_analytical_error_maxabs_u;
		double benchmark_analytical_error_maxabs_v;
	};

	BenchmarkErrors benchmark;



public:

	Parareal_SimulationInstance()
	{
	};

#if SWEET_PARAREAL_CART2D
	// Cart2D
	void setup(	SimulationVariables* i_simVars,
			std::vector<Cart2DDataConfig*> i_cart2DDataConfig,
			std::vector<Cart2DOperators*> i_op_cart2d,
			t_tsmType* i_timeSteppersFine,
			t_tsmType* i_timeSteppersCoarse)
	{
		cart2DDataConfig = i_cart2DDataConfig;
		op_cart2d = i_op_cart2d;
		setup(i_simVars, 
				i_timeSteppersFine, i_timeSteppersCoarse);

		// IMPORTANT: setup initial conditions (inside setup) before setting up timesteppers
		// because simulation parameters may change
		timeSteppersFine->setup(
				simVars->disc.timestepping_method,
				simVars->disc.timestepping_order,
				simVars->disc.timestepping_order2,
				*op_cart2d[0],
				*simVars
			);

		timeSteppersCoarse->setup(
				simVars_coarse->disc.timestepping_method,
				simVars_coarse->disc.timestepping_order,
				simVars_coarse->disc.timestepping_order2,
				*op_cart2d[1],
				*simVars_coarse
			);

		PInT_Common::setup();

	#if SWEET_PARAREAL_CART2D_BURGERS
		PInT_Common::set_tsm_burgers(i_timeSteppersFine);
	#endif


		if (simVars->benchmark.benchmark_name == "normalmodes" )
			compute_normal_modes = true;

	}

#elif SWEET_PARAREAL_SPHERE2D
	// Sphere2D
	void setup(	SimulationVariables* i_simVars,
			std::vector<Sphere2DData_Config*> i_sphere2DDataConfig,
			std::vector<Sphere2DOperators*> i_op_sphere2D,
			std::vector<Sphere2DOperators*> i_op_sphere2d_nodealiasing,
			t_tsmType* i_timeSteppersFine,
			t_tsmType* i_timeSteppersCoarse)
	{
		sphere2DDataConfig = i_sphere2DDataConfig;
		op_sphere2D = i_op_sphere2D;
		op_sphere2d_nodealiasing = i_op_sphere2d_nodealiasing;
		setup(i_simVars,
				i_timeSteppersFine, i_timeSteppersCoarse);

		// IMPORTANT: setup initial conditions (inside setup) before setting up timesteppers
		// because simulation parameters may change
		timeSteppersFine->setup(
					simVars->disc.timestepping_method,
					////simVars->disc.timestepping_order,
					////simVars->disc.timestepping_order2,
					*op_sphere2D[0],
					*simVars
				);
		timeSteppersCoarse->setup(
					simVars_coarse->disc.timestepping_method,
					////simVars_coarse->disc.timestepping_order,
					////simVars_coarse->disc.timestepping_order2,
					*op_sphere2D[1],
					*simVars_coarse
				);

		PInT_Common::setup();
	}
#endif

	void setup(SimulationVariables* i_simVars,
			t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{

		simVars = i_simVars;
		simVars_coarse = new SimulationVariables;
		*simVars_coarse = *simVars;
		simVars_coarse->disc.timestepping_method = simVars->parareal.coarse_timestepping_method;
		simVars_coarse->disc.timestepping_order = simVars->parareal.coarse_timestepping_order;
		simVars_coarse->disc.timestepping_order2 = simVars->parareal.coarse_timestepping_order2;
		simVars_coarse->timecontrol.current_timestepSize = simVars->parareal.coarse_timestepSize;

		timeSteppersFine = i_timeSteppersFine;
		timeSteppersCoarse = i_timeSteppersCoarse;


#if SWEET_PARAREAL_SCALAR
		timeSteppersFine->setup(*simVars);
		timeSteppersCoarse->setup(*simVars);
#endif

		// All containers contain data defined on the fine spatial grid
		// (including the coarse_data container, with interpolation performed in run_timestep_coarse)
		// Exception: the auxiliary coarse solutions for SL tsm, since they are not summed up with fine solutions
		parareal_data_start                      = create_new_data_container("fine");
		parareal_data_fine                       = create_new_data_container("fine");
		parareal_data_coarse                     = create_new_data_container("fine");
		parareal_data_coarse_coarse_mesh         = create_new_data_container("coarse");
		parareal_data_output                     = create_new_data_container("fine");
		parareal_data_error                      = create_new_data_container("fine");
		parareal_data_coarse_previous_timestep   = create_new_data_container("coarse");
		parareal_data_coarse_previous_time_slice = create_new_data_container("coarse");
		parareal_data_fine_previous_timestep     = create_new_data_container("fine");
		parareal_data_fine_previous_time_slice   = create_new_data_container("fine");
		parareal_data_ref_exact                  = create_new_data_container("fine");
		parareal_data_fine_exact                 = create_new_data_container("fine");
#if SWEET_DEBUG
		parareal_data_fine_exact_debug           = create_new_data_container("fine");
#endif

		parareal_data_debug           = create_new_data_container("coarse");

		sim_setup_initial_data();
	};


public:

	Parareal_GenericData* create_new_data_container(std::string level)
	{
#if SWEET_PARAREAL_SCALAR
		{
			Parareal_GenericData_Scalar<N>* out = new Parareal_GenericData_Scalar<N>;
			out->allocate_data();
			out->set_time(timeframe_end);
			return out;
		}

#elif SWEET_PARAREAL_CART2D
		{
			Parareal_GenericData_Cart2DData_Spectral<N>* out = new Parareal_GenericData_Cart2DData_Spectral<N>;
			if (level == "fine")
				out->setup_data_config(cart2DDataConfig[0]);
			else if (level == "coarse")
				out->setup_data_config(cart2DDataConfig[1]);
			else
				SWEETErrorFatal("Wrong level.");
			out->allocate_data();
			return out;
		}

#elif SWEET_PARAREAL_SPHERE2D
		{
			Parareal_GenericData_Sphere2DData_Spectral<N>* out = new Parareal_GenericData_Sphere2DData_Spectral<N>;
			if (level == "fine")
				out->setup_data_config(sphere2DDataConfig[0]);
			else if (level == "coarse")
				out->setup_data_config(sphere2DDataConfig[1]);
			else
				SWEETErrorFatal("Wrong level.");
			out->allocate_data();
			return out;
		}
#endif

		// default time set
	}

public:

	/**
	 * Check if the time slice contains an integer number of coarse and fine time steps
	 */
	void sim_check_timesteps(
			double time_slice_size
	)
	{

		// check if seup has been called
		SWEET_ASSERT(simVars_coarse->timecontrol.current_timestepSize == simVars->parareal.coarse_timestepSize);

		// check if each time slice contains an integer number of fine and coarse time steps
		double eps = 1e-12;
		double mod_coarse = fmod(time_slice_size, simVars_coarse->timecontrol.current_timestepSize);
		double mod_fine = fmod(time_slice_size, simVars->timecontrol.current_timestepSize);
				if ( std::abs(mod_coarse) > eps && std::abs(mod_coarse - time_slice_size) > eps )
			SWEETErrorFatal("Time slice length must be an integer multiple of the coarse time step! (" + std::to_string(simVars_coarse->timecontrol.current_timestepSize) + ", " + std::to_string(time_slice_size) + ")");
				if ( std::abs(mod_fine) > eps && std::abs(mod_fine - time_slice_size) > eps )
		{
			std::cout << "Number of timesteps: " << nb_timesteps_fine << std::endl;
			std::cout << "Mod(dt): " << std::abs(mod_fine) << std::endl;
			std::cout << "Mod(dt) - Dt: " << std::abs(mod_fine - time_slice_size) << std::endl;
			SWEETErrorFatal("Time slice length must be an integer multiple of the fine time step! (" + std::to_string(simVars->timecontrol.current_timestepSize) + ", " + std::to_string(time_slice_size) + ")");
		}
	};


	/**
	 * Set the start and end of the coarse time step
	 */
	void sim_set_timeframe(
			double i_timeframe_start,	//!< start timestamp of coarse time step
			double i_timeframe_end		//!< end time stamp of coarse time step
	){
		// check if seup has been called
		SWEET_ASSERT(simVars_coarse->timecontrol.current_timestepSize == simVars->parareal.coarse_timestepSize);

		if (simVars->parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;
		timeframe_start = i_timeframe_start;
		timeframe_end = i_timeframe_end;

		nb_timesteps_fine = (int)((timeframe_end - timeframe_start) / simVars->timecontrol.current_timestepSize);
		nb_timesteps_coarse = (int)((timeframe_end - timeframe_start) / simVars_coarse->timecontrol.current_timestepSize);
		if (timeframe_start + nb_timesteps_fine * simVars->timecontrol.current_timestepSize < timeframe_end - 1e-14)
			nb_timesteps_fine++;
		if (timeframe_start + nb_timesteps_coarse * simVars_coarse->timecontrol.current_timestepSize < timeframe_end - 1e-14)
			nb_timesteps_coarse++;
		SWEET_ASSERT( std::abs(timeframe_start + nb_timesteps_fine * simVars->timecontrol.current_timestepSize - timeframe_end) < 1e-14);
		SWEET_ASSERT( std::abs(timeframe_start + nb_timesteps_coarse * simVars_coarse->timecontrol.current_timestepSize - timeframe_end) < 1e-14);

		dt_fine = simVars->timecontrol.current_timestepSize;
		dt_coarse = simVars->parareal.coarse_timestepSize;

		// set time to parareal_genericdata instances
		parareal_data_start->set_time(i_timeframe_end);
		parareal_data_fine->set_time(i_timeframe_end);
		parareal_data_coarse->set_time(i_timeframe_end);
		parareal_data_coarse_coarse_mesh->set_time(i_timeframe_end);
		parareal_data_output->set_time(i_timeframe_end);
		parareal_data_error->set_time(i_timeframe_end);
		parareal_data_coarse_previous_timestep->set_time(i_timeframe_end);
		parareal_data_coarse_previous_time_slice->set_time(i_timeframe_end);
		parareal_data_fine_previous_timestep->set_time(i_timeframe_end);
		parareal_data_fine_previous_time_slice->set_time(i_timeframe_end);
		parareal_data_ref_exact->set_time(i_timeframe_end);
		parareal_data_fine_exact->set_time(i_timeframe_end);
#if SWEET_DEBUG
		parareal_data_fine_exact_debug->set_time(i_timeframe_end);
#endif

	};

	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		///reset();

#if SWEET_PARAREAL_SCALAR
		double u0 = atof(simVars->user_defined.var[1].c_str());
		parareal_data_start->dataArrays_2_GenericData_Scalar(u0);

#elif SWEET_PARAREAL_CART2D
		Cart2DData_Spectral t0_prog_h_pert(cart2DDataConfig[0]);
		Cart2DData_Spectral t0_prog_u(cart2DDataConfig[0]);
		Cart2DData_Spectral t0_prog_v(cart2DDataConfig[0]);

	#if SWEET_PARAREAL_CART2D_SWE
		SWECart2DBenchmarksCombined sweCart2DBenchmarks;
		sweCart2DBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, *simVars, *op_cart2d[0]);
	#elif SWEET_PARAREAL_CART2D_BURGERS
		Cart2D_DataGrid t0_prog_u_phys = t0_prog_u.toGrid();
		Cart2D_DataGrid t0_prog_v_phys = t0_prog_v.toGrid();
		if (simVars->disc.space_grid_use_c_staggering)
		{
			t0_prog_u_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.cart2d_domain_size[0];
					double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.cart2d_domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
				}
			);
			t0_prog_v_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
				io_data = 0.0;
				}
			);
		}
		else
		{
			t0_prog_u_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.cart2d_domain_size[0];
					double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.cart2d_domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
				}
			);
	
			t0_prog_v_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
				io_data = 0.0;
				}
			);
		}
		t0_prog_u.loadCart2DDataGrid(t0_prog_u_phys);
		t0_prog_v.loadCart2DDataGrid(t0_prog_v_phys);
	#endif


		parareal_data_start->dataArrays_2_GenericData_Cart2DData_Spectral(
	#if SWEET_PARAREAL_CART2D_SWE
											t0_prog_h_pert,
	#endif
											t0_prog_u,
											t0_prog_v);


		if (simVars->parareal.spatial_coarsening)
			parareal_data_coarse_previous_time_slice->restrict(*parareal_data_start);
		else
			parareal_data_coarse_previous_time_slice->dataArrays_2_GenericData_Cart2DData_Spectral(
		#if SWEET_PARAREAL_CART2D_SWE
												t0_prog_h_pert,
		#endif
												t0_prog_u,
												t0_prog_v);

			parareal_data_fine_previous_time_slice->dataArrays_2_GenericData_Cart2DData_Spectral(
	#if SWEET_PARAREAL_CART2D_SWE
											t0_prog_h_pert,
	#endif
											t0_prog_u,
											t0_prog_v);

#elif SWEET_PARAREAL_SPHERE2D
		sweet::Data::Sphere2D::DataSpectral t0_prog_phi_pert(sphere2DDataConfig[0]);
		sweet::Data::Sphere2D::DataSpectral t0_prog_vrt(sphere2DDataConfig[0]);
		sweet::Data::Sphere2D::DataSpectral t0_prog_div(sphere2DDataConfig[0]);

		BenchmarksSphere2DSWE sphere2DBenchmarks;
		sphere2DBenchmarks.setup(*simVars, *op_sphere2D[0]);
		sphere2DBenchmarks.master->get_initial_state(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);

		parareal_data_start->dataArrays_2_GenericData_Sphere2DData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
		if (simVars->parareal.spatial_coarsening)
			parareal_data_coarse_previous_time_slice->restrict(*parareal_data_start);
		else
			parareal_data_coarse_previous_time_slice->dataArrays_2_GenericData_Sphere2DData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
		parareal_data_fine_previous_time_slice->dataArrays_2_GenericData_Sphere2DData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
#endif
	};


	/**
	 * Set simulation data to data given in i_sim_data.
	 * This can be data which is computed by another simulation.
	 * Y^S := i_sim_data
	 */
	void sim_set_data(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data()" << std::endl;

		// copy to buffers
		*parareal_data_start = i_pararealData;

		parareal_data_start->set_time(timeframe_end);
	};

	/**
	 * Set solution of penult coarse timestep of previous time slice
	 */
	void sim_set_data_coarse_previous_time_slice(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data_coarse_previous_time_slice()" << std::endl;

		// copy to buffers
		*parareal_data_coarse_previous_time_slice = i_pararealData;

		parareal_data_coarse_previous_time_slice->set_time(timeframe_end);
	};

	/**
	 * Set solution of penult fine timestep of previous time slice
	 */
	void sim_set_data_fine_previous_time_slice(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data_fine_previous_time_slice()" << std::endl;

		// copy to buffers
		*parareal_data_fine_previous_time_slice = i_pararealData;

		parareal_data_fine_previous_time_slice->set_time(timeframe_end);
	};

#if SWEET_PARAREAL == 2
	void sim_set_data_coarse(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data_coarse()" << std::endl;

		// copy to buffers
		*parareal_data_coarse = i_pararealData;

		parareal_data_coarse->set_time(timeframe_end);
	};

#endif

#if SWEET_DEBUG
	/**
	* Store exact solution (full fine simulation) at the end of the time slice
	*/
	void sim_set_data_fine_exact(
			Parareal_GenericData &i_pararealData
	)
	{
		*parareal_data_fine_exact_debug = i_pararealData;

		parareal_data_fine_exact_debug->set_time(timeframe_end);
	};

	/**
	* Check if solution at time k (end of time slice k-1) is exact (= fine) at iteration k
	*/
	void compare_2_fine_exact(
	)
	{
		Parareal_GenericData* diff = create_new_data_container("fine");
		*diff = *parareal_data_fine_exact_debug;
		*diff -= *parareal_data_output;
		///parareal_data_fine_exact_debug->grid_print();
		///parareal_data_output->grid_print();
		std::cout << "Max serial:" << parareal_data_fine_exact_debug->spectral_reduce_maxAbs() << std::endl;
		std::cout << "Max parareal:" << parareal_data_output->spectral_reduce_maxAbs() << std::endl;
		std::cout << "DIFF: " << diff->spectral_reduce_maxAbs() << std::endl;
		SWEET_ASSERT(diff->spectral_reduce_maxAbs() < 1e-10);
		delete diff;
	};

#endif

	void set_previous_solution(
				std::string tsm_level
				)
	{

		if (tsm_level == "fine")
			timeSteppersFine->master->set_previous_solution(parareal_data_fine_previous_time_slice);
		else if (tsm_level == "coarse")
			timeSteppersCoarse->master->set_previous_solution(parareal_data_coarse_previous_time_slice);
		else
			SWEETErrorFatal("Wrong tsm_level (should be 'fine' or 'coarse')");

	};

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	void sim_set_mpi_comm(
			int i_mpi_comm
	)
	{
	};


	void runTimestep(
			Parareal_GenericData* io_data,
			std::string tsm_level
	)
	{

		if (tsm_level == "fine")
			timeSteppersFine->master->runTimestep(
						io_data,
						simVars->timecontrol.current_timestepSize,
						simVars->timecontrol.current_simulation_time
					);
		else if (tsm_level == "coarse")
			timeSteppersCoarse->master->runTimestep(
						io_data,
						simVars_coarse->timecontrol.current_timestepSize,
						simVars_coarse->timecontrol.current_simulation_time
					);
		else
			SWEETErrorFatal("Wrong tsm_level (should be 'fine' or 'coarse')");

	}



	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	void run_timestep_fine()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		// reset simulation time
		simVars->timecontrol.current_simulation_time = timeframe_start;
		simVars->timecontrol.max_simulation_time = timeframe_end;
		simVars->timecontrol.current_timestep_nr = 0;
		simVars->timecontrol.current_timestepSize = dt_fine;

		//std::cout << simVars->disc.timestepping_method << std::endl;
		///std::cout << SL_tsm.size() << std::endl;
		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(SL_tsm.begin(), SL_tsm.end(), simVars->disc.timestepping_method) != SL_tsm.end())
		{
			set_previous_solution("fine");
		}

		*(parareal_data_fine) = *(parareal_data_start);

		int nb_timesteps = 0;
		while (nb_timesteps != nb_timesteps_fine)
		{
			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*(parareal_data_fine_previous_timestep) = *(parareal_data_fine);

			runTimestep(parareal_data_fine, "fine");

			simVars->timecontrol.current_simulation_time += simVars->timecontrol.current_timestepSize;
			SWEET_ASSERT(simVars->timecontrol.current_simulation_time <= timeframe_end + 1e-14);
			nb_timesteps++;
		}
	};


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_GenericData& get_reference_2_data_timestep_fine()
	{
		return *(parareal_data_fine);
	};


	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "run_timestep_coarse()" << std::endl;

		// reset simulation time

		simVars_coarse->timecontrol.current_simulation_time = timeframe_start;
		simVars_coarse->timecontrol.max_simulation_time = timeframe_end;
		simVars_coarse->timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(SL_tsm.begin(), SL_tsm.end(), simVars_coarse->disc.timestepping_method) != SL_tsm.end())
			set_previous_solution("coarse");

		*(parareal_data_coarse) = *(parareal_data_start);

		// interpolate to coarse spatial mesh
		if (simVars->parareal.spatial_coarsening)
			parareal_data_coarse_coarse_mesh->restrict(*parareal_data_coarse);
		else
			*parareal_data_coarse_coarse_mesh = *parareal_data_coarse;

		int nb_timesteps = 0;
		while (nb_timesteps != nb_timesteps_coarse)
		{
			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*(parareal_data_coarse_previous_timestep) = *(parareal_data_coarse_coarse_mesh);

			runTimestep(parareal_data_coarse_coarse_mesh, "coarse");
			simVars_coarse->timecontrol.current_simulation_time += simVars_coarse->timecontrol.current_timestepSize;
			SWEET_ASSERT(simVars_coarse->timecontrol.current_simulation_time <= timeframe_end +  1e-14);
			nb_timesteps++;
		}

		// interpolate to coarse spatial mesh
		if (simVars->parareal.spatial_coarsening)
			parareal_data_coarse->pad_zeros(*parareal_data_coarse_coarse_mesh);
		else
			*parareal_data_coarse = *parareal_data_coarse_coarse_mesh;

	};


	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_GenericData& get_reference_2_data_timestep_coarse()
	{
		return *(parareal_data_coarse);
	};

	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	void compute_difference()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "compute_difference()" << std::endl;

		*(parareal_data_error) = *(parareal_data_fine);
		*(parareal_data_error) -= *(parareal_data_coarse);
	};

	/**
	 * return the difference between fine and coarse solution:
	 */
	Parareal_GenericData& get_reference_2_data_timestep_diff()
	{
		return *(parareal_data_error);
	};

	void sim_set_data_diff(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data_diff()" << std::endl;

		// copy to buffers
		*parareal_data_error = i_pararealData;

		parareal_data_error->set_time(timeframe_end);
	};


	/**
	 * return the penult solution after the coarse propagation:
	 */
	Parareal_GenericData& get_reference_2_data_timestep_coarse_previous_timestep()
	{
		return *(parareal_data_coarse_previous_timestep);
	};

	/**
	 * return the penult solution after the fine propagation:
	 */
	Parareal_GenericData& get_reference_2_data_timestep_fine_previous_timestep()
	{
		return *(parareal_data_fine_previous_timestep);
	};


	/**
	 * Compute the data to be forwarded to the next time step
	 * Y^O := Y^C + Y^E
	 *
	 * Return: If true, the error indicator based on the computed error norm between the
	 * old values and new values
	 */
	double compute_output_data(
			bool i_compute_convergence_test
	)
	{
		double convergence = -1;

		if (!i_compute_convergence_test)
		//if (!i_compute_convergence_test || !output_data_valid)
		{
			*(parareal_data_output) = *(parareal_data_coarse);
						*(parareal_data_output) -= *(parareal_data_error);

			//output_data_valid = true;
			return convergence;
		}

		// compute output data
		Parareal_GenericData* tmp = create_new_data_container("fine");
		*tmp = *(parareal_data_coarse);
				*tmp += *(parareal_data_error);

		// compute difference w.r.t. previous output data
		Parareal_GenericData* tmp2 = create_new_data_container("fine");
		*tmp2 = *(parareal_data_output);
		*tmp2 -= *tmp;
		convergence = tmp2->spectral_reduce_maxAbs();

		// store output data
		*(parareal_data_output) = *tmp;

		delete tmp;
		delete tmp2;

		//output_data_valid = true;
		return convergence;

	};


	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_GenericData& get_reference_2_output_data()
	{
		return *(parareal_data_output);
	};


	void output_data_file(
				int iteration_id,
				int time_slice_id,
				bool output_initial_data = false
	)
	{
		if (output_initial_data)
			PInT_Common::output_data_file(
						parareal_data_start,
						iteration_id,
						time_slice_id - 1,
						timeframe_start
			);
		else
			PInT_Common::output_data_file(
						parareal_data_output,
						iteration_id,
						time_slice_id,
						timeframe_end
			);

	}

	void store_parareal_error(
			int iteration_id,
			int time_slice_id,
			std::string path_ref,
			std::string base_solution,	// "ref" or "fine"
			int i_precision = 16
	)
	{

		Parareal_GenericData* parareal_data_ref;
		if (base_solution == "ref")
			parareal_data_ref = parareal_data_ref_exact;
		else if (base_solution == "fine")
			parareal_data_ref = parareal_data_fine_exact;
		else
			SWEETErrorFatal("Wrong base solution for computing parareal errors.");

		int nvar = N;

		std::cout << "ALSKDJASD" << std::endl;
		exit(1);
		PInT_Common::store_pint_error(
						parareal_data_output,
						parareal_data_ref,
						nvar,
						iteration_id,
						time_slice_id,
						timeframe_end,
						path_ref,
						base_solution,
						"parareal",
						i_precision
					);

	}



	void output_data_console(
			int iteration_id,
			int time_slice_id
	)
	{
	};


	void check_for_nan_parareal()
	{
		if (parareal_data_output->check_for_nan())
			SWEETErrorFatal("Instability detected in parareal!");
	};

	void delete_data_container(
					Parareal_GenericData* i_data
	)
	{
		if (i_data)
		{
			delete i_data;
			i_data = nullptr;
		}
	}

	~Parareal_SimulationInstance()
	{

			delete_data_container(parareal_data_start);
			delete_data_container(parareal_data_fine);
			delete_data_container(parareal_data_coarse);
			delete_data_container(parareal_data_coarse_coarse_mesh);
			delete_data_container(parareal_data_output);
			delete_data_container(parareal_data_error);
			delete_data_container(parareal_data_coarse_previous_timestep);
			delete_data_container(parareal_data_coarse_previous_time_slice);
			delete_data_container(parareal_data_fine_previous_timestep);
			delete_data_container(parareal_data_fine_previous_time_slice);
			delete_data_container(parareal_data_ref_exact);
			delete_data_container(parareal_data_fine_exact);
	#if SWEET_DEBUG
			delete_data_container(parareal_data_fine_exact_debug);
	#endif


		if (simVars_coarse)
		{
			delete simVars_coarse;
			simVars_coarse = nullptr;
		}
	}
};

}}

#endif
