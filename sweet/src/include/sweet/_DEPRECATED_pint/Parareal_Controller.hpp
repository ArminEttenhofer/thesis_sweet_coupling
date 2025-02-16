/*
 * PararealController.hpp
 *
 *  Created on: 11 Apr 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_CONTROLLER_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_CONTROLLER_HPP

// Checking if geometry and model have been correctly defined
#if SWEET_PARAREAL_SCALAR
	#if SWEET_PARAREAL_CART2D || SWEET_PARAREAL_SPHERE2D
		#error "More than one geometry has been defined for parareal"
	#endif
#elif SWEET_PARAREAL_CART2D
	#if SWEET_PARAREAL_SCALAR || SWEET_PARAREAL_SPHERE2D
		#error "More than one geometry has been defined for parareal"
	#endif
	#if (!SWEET_PARAREAL_CART2D_SWE) && (!SWEET_PARAREAL_CART2D_BURGERS)
		#error "No model has been defined for parareal on the cart2d"
	#endif
	#if SWEET_PARAREAL_CART2D_SWE && SWEET_PARAREAL_CART2D_BURGERS
		#error "More than one model has been defined for parareal on the cart2d"
	#endif
#elif SWEET_PARAREAL_SPHERE2D
	#if SWEET_PARAREAL_SCALAR || SWEET_PARAREAL_CART2D
		#error "More than one geometry has been defined for parareal"
	#endif
#else
	#error "No geometry has been defined for parareal"
#endif

#if SWEET_PARAREAL==2
#include <mpi.h>
#endif

#include <sweet/_DEPRECATED_pint/Parareal_ConsolePrefix.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_SimulationInstance.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_SimulationVariables.hpp>

#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>

#if SWEET_PARAREAL_SCALAR
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Scalar.hpp>

#elif SWEET_PARAREAL_CART2D
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Cart2DData_Spectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>

#elif SWEET_PARAREAL_SPHERE2D
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Sphere2DData_Spectral.hpp>
#include <sweet/core/sphere2D/Sphere2DOperators.hpp>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <map>


namespace sweet {
namespace DEPRECATED_pint {


/**
 * This class takes over the control and
 * calls methods offered via PararealSimulation.
 *
 * \param t_SimulationInstance	class which implements the Parareal_SimulationInstance interfaces
 */
////template <class t_SimulationInstance>
template <class t_tsmType, int N>
class Parareal_Controller
{

	/**
	 * Pointers to interfaces of simulationInstances
	 * This helps to clearly separate between the allocation of the simulation classes and the parareal interfaces.
	 */
	std::vector<Parareal_SimulationInstance<t_tsmType, N>*> parareal_simulationInstances = {};

	SimulationVariables* simVars;


	// Operators and DataConfig
#if SWEET_PARAREAL_CART2D
	std::vector<Cart2DDataConfig*> cart2DDataConfig;
	std::vector<Cart2DOperators*> op_cart2d;
#elif SWEET_PARAREAL_SPHERE2D
	std::vector<Sphere2DData_Config*> sphere2DDataConfig;
	std::vector<Sphere2DOperators*> op_sphere2D;
	std::vector<Sphere2DOperators*> op_sphere2d_nodealiasing;
#endif


	t_tsmType* timeSteppersFine = nullptr;
	t_tsmType* timeSteppersCoarse = nullptr;


	std::vector<bool> timeframe_do_output = {};

	/**
	 * Pointer to parareal simulation variables.
	 * These variables are used as a singleton
	 */
	Parareal_SimulationVariables *pVars;

	/**
	 * Class which helps prefixing console output
	 */
	Parareal_ConsolePrefix CONSOLEPREFIX;


// MPI
	int mpi_nprocs = 1;
	int mpi_rank = 0;
	std::vector<int> slices_for_proc = {};
	std::vector<int> proc_for_slices = {};
	std::map<int, int> global_2_local_slice;
	int buffer_size;



public:

#if SWEET_PARAREAL==2

	// every communication is made between mpi_proc = 0 and mpi_proc > 0
	void communicate_solution( Parareal_GenericData* io_data, int source_rank, int dest_rank, int tag)
	{

		bool send = false;
		bool recv = false;

		if (mpi_rank == source_rank)
			send = true;
		else if (mpi_rank == dest_rank)
			recv = true;
		else
			return;


	#if SWEET_PARAREAL_SCALAR
		double* serial_data = sweet::Memory::MemBlockAlloc::alloc<double>(N * sizeof(double));
	#elif SWEET_PARAREAL_CART2D
		std::complex<double>* serial_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(N * cart2DDataConfig[0]->spectral_array_data_number_of_elements * sizeof(std::complex<double>));
	#elif SWEET_PARAREAL_SPHERE2D
		std::complex<double>* serial_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(N * sphere2DDataConfig[0]->spectral_array_data_number_of_elements * sizeof(std::complex<double>));
	#endif


		int size = io_data->size();
		///int size = buffer_size;
		tag = tag + 1000 * source_rank + 100 * dest_rank;

		if (send)
		{

			std::cout << "Proc " << mpi_rank << " is sending solution to proc " << dest_rank << " with tag " << tag << std::endl;

			io_data->serialize(serial_data);
	#if SWEET_PARAREAL_SCALAR
			MPI_Send(serial_data, size, MPI_DOUBLE, dest_rank, tag + 0, MPI_COMM_WORLD);
	#else
			MPI_Send(serial_data, size, MPI_DOUBLE_COMPLEX, dest_rank, tag + 0, MPI_COMM_WORLD);
	#endif

		}
		else if (recv)
		{

			std::cout << "Proc " << mpi_rank << " is receiving solution from proc " << source_rank << " with tag " << tag << std::endl;

			MPI_Status status;
	#if SWEET_PARAREAL_SCALAR
			MPI_Recv(serial_data, size, MPI_DOUBLE, source_rank, tag + 0, MPI_COMM_WORLD, &status);
	#else
			MPI_Recv(serial_data, size, MPI_DOUBLE_COMPLEX, source_rank, tag + 0, MPI_COMM_WORLD, &status);
	#endif
			io_data->deserialize(serial_data);

		}

	#if SWEET_PARAREAL_SCALAR
		sweet::Memory::MemBlockAlloc::free(serial_data, N * sizeof(double));
	#elif SWEET_PARAREAL_CART2D
		sweet::Memory::MemBlockAlloc::free(serial_data, N * cart2DDataConfig[0]->grid_number_elements * sizeof(std::complex<double>));
	#elif SWEET_PARAREAL_SPHERE2D
		sweet::Memory::MemBlockAlloc::free(serial_data, N * sphere2DDataConfig[0]->grid_number_elements * sizeof(std::complex<double>));
	#endif


	}


#endif


#if SWEET_PARAREAL_SCALAR
	// Scalar
	Parareal_Controller(SimulationVariables* i_simVars,
						t_tsmType* i_timeSteppersFine,
						t_tsmType* i_timeSteppersCoarse):
		simVars(i_simVars),
		timeSteppersFine(i_timeSteppersFine),
		timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};

#elif SWEET_PARAREAL_CART2D
	// Cart2D
	Parareal_Controller(SimulationVariables* i_simVars,
						std::vector<Cart2DDataConfig*> i_cart2DDataConfig,
						std::vector<Cart2DOperators*> i_op_cart2d,
						t_tsmType* i_timeSteppersFine,
						t_tsmType* i_timeSteppersCoarse):
		simVars(i_simVars),
		cart2DDataConfig(i_cart2DDataConfig),
		op_cart2d(i_op_cart2d),
		timeSteppersFine(i_timeSteppersFine),
		timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};

#elif SWEET_PARAREAL_SPHERE2D
	// Sphere2D
	Parareal_Controller(SimulationVariables* i_simVars,
						std::vector<Sphere2DData_Config*> i_sphere2DDataConfig,
						std::vector<Sphere2DOperators*> &i_op_sphere2D,
						std::vector<Sphere2DOperators*> &i_op_sphere2d_nodealiasing,
						t_tsmType* i_timeSteppersFine,
						t_tsmType* i_timeSteppersCoarse):
		simVars(i_simVars),
		sphere2DDataConfig(i_sphere2DDataConfig),
		op_sphere2D(i_op_sphere2D),
		op_sphere2d_nodealiasing(i_op_sphere2d_nodealiasing),
		timeSteppersFine(i_timeSteppersFine),
		timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};
#endif

	~Parareal_Controller()
	{
		cleanup();
	}


	inline
	void CONSOLEPREFIX_start(
			const char *i_prefix
	)
	{
//		if (pVars->verbosity > 0)
			CONSOLEPREFIX.start(i_prefix);
	}


	inline
	void CONSOLEPREFIX_start(int i_prefix)
	{
//		if (pVars->verbosity > 0)
			CONSOLEPREFIX.start(i_prefix);
	}

	inline
	void CONSOLEPREFIX_end()
	{
		if (pVars->verbosity > 0)
			CONSOLEPREFIX.end();
	}


	void cleanup()
	{
		for (typename std::vector<Parareal_SimulationInstance<t_tsmType, N>*>::iterator it = parareal_simulationInstances.begin();
															it != parareal_simulationInstances.end();
															it++)
			if (*it)
				delete *it;
	}


	void setup(
			//Parareal_SimulationVariables *i_pararealSimVars
	)
	{

		cleanup();

#if SWEET_PARAREAL==2
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
#endif

		pVars = &simVars->parareal;

		if (!pVars->enabled)
			return;

		if (pVars->coarse_slices <= 0)
		{
			std::cerr << "Invalid number of coarse slices" << std::endl;
			exit(1);
		}

		if (pVars->coarse_slices % mpi_nprocs != 0)
		{
			SWEETErrorFatal("Number of coarse slices must be a multiple integer of the number of MPI processes");
		}

		if (pVars->max_simulation_time <= 0)
		{
			std::cerr << "Invalid simulation time" << std::endl;
			exit(1);
		}

		// allocate raw simulation instances
		////simulationInstances = new t_SimulationInstance[pVars->coarse_slices];

		CONSOLEPREFIX.start("[MAIN] ");
		std::cout << "Resetting simulation instances" << std::endl;


		// Distribute slices for each MPI proc
		// mpi_rank = 0 contains all processes to compute all serial parts
		int slices_per_proc = pVars->coarse_slices / mpi_nprocs;

		if (mpi_rank == 0)
			for (int k = 0; k < pVars->coarse_slices; k++)
				global_2_local_slice.emplace(std::make_pair(k, k)); // slice k is the i-th slice treated by proc mpi_rank

		int i = 0;
		for (int k = mpi_rank * slices_per_proc; k < (mpi_rank + 1) * slices_per_proc; k++)
		{
			slices_for_proc.push_back(k); // slice k is treated by proc mpi_rank
			if (mpi_rank > 0)
				global_2_local_slice.emplace(std::make_pair(k, i)); // slice k is the i-th slice treated by proc mpi_rank
			i++;
		}
		if (mpi_rank > 0)
			for (int k = 0; k < pVars->coarse_slices; k++)
				if (k < mpi_rank * slices_per_proc || k >= (mpi_rank + 1) * slices_per_proc)
					global_2_local_slice.emplace(std::make_pair(k, -1)); // slice k is not treated by proc mpi_rank

		///! Distribute MPI proc for each slice
		proc_for_slices = std::vector<int>(pVars->coarse_slices);
		if (mpi_rank == 0)
			for (int i = 0; i < mpi_nprocs; i++)
				for (int k = i * slices_per_proc; k < (i + 1) * slices_per_proc; k++)
					proc_for_slices[k] = i;

#if SWEET_PARAREAL==2
		MPI_Bcast(&proc_for_slices[0], proc_for_slices.size() , MPI_INT, 0, MPI_COMM_WORLD );
#endif


		// size of coarse time step
		double time_slice_size = pVars->max_simulation_time / pVars->coarse_slices;
		if (pVars->coarse_timestepSize < 0)
			pVars->coarse_timestepSize = time_slice_size;


		// convert to pararealsimulationInstances to get Parareal interfaces
		for (int k = 0; k < pVars->coarse_slices; k++)
		{

			// except for mpi_rank = 0, only create instances for slices treated by mpi_rank
			if (mpi_rank > 0)
				if (k < slices_for_proc[0] || k > slices_for_proc.back())
					continue;

			CONSOLEPREFIX_start(k);

			int local_k = global_2_local_slice.at(k);
			parareal_simulationInstances.push_back(new Parareal_SimulationInstance<t_tsmType, N>);
			std::cout << "mpi_rank " << mpi_rank << " setting up instance " << k << " " << local_k << std::endl;
#if SWEET_PARAREAL_SCALAR
				parareal_simulationInstances[local_k]->setup(simVars,
									   timeSteppersFine,
									   timeSteppersCoarse);

#elif SWEET_PARAREAL_CART2D
				parareal_simulationInstances[local_k]->setup(simVars,
									   cart2DDataConfig,
									   op_cart2d,
									   timeSteppersFine,
									   timeSteppersCoarse);

#elif SWEET_PARAREAL_SPHERE2D
				parareal_simulationInstances[local_k]->setup(simVars,
									   sphere2DDataConfig,
									   op_sphere2D,
									   op_sphere2d_nodealiasing,
									   timeSteppersFine,
									   timeSteppersCoarse);
#endif
		}


		CONSOLEPREFIX_start("[MAIN] ");
		std::cout << "Setup time frames" << std::endl;

		/*
		 * SETUP time frame
		 */


		// if time slices are not homogeneous, this should be called by each parareal_simulationInstance
		//for (int k = 0; k < pVars->coarse_slices; k++)
		//{
		CONSOLEPREFIX_start(0);
		parareal_simulationInstances[0]->sim_check_timesteps(time_slice_size);
		//}

		for (int k = 0; k < pVars->coarse_slices; k++)
		////for (int k = slices_for_proc[0]; k <= slices_for_proc.back(); k++)
		{
			if (mpi_rank > 0)
				if (k < slices_for_proc[0] || k > slices_for_proc.back())
					continue;

			CONSOLEPREFIX_start(k);
			int local_k = global_2_local_slice.at(k);
			parareal_simulationInstances[local_k]->sim_set_timeframe(time_slice_size*k, time_slice_size*(k+1));
			timeframe_do_output.push_back(check_do_output(time_slice_size*(k+1)));
		}

		/*
		 * Setup first simulation instance
		 */
		if (mpi_rank == 0)
		{
			CONSOLEPREFIX_start(0);
			parareal_simulationInstances[0]->sim_setup_initial_data();
			CONSOLEPREFIX_end();

			CONSOLEPREFIX_start("[MAIN] ");
			std::cout << "Finished setup parareal" << std::endl;
			CONSOLEPREFIX_end();
		}

#if SWEET_PARAREAL == 2
		if (mpi_rank == 0)
			buffer_size = parareal_simulationInstances[0]->parareal_data_start->size();
		MPI_Bcast(&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD );

		MPI_Barrier(MPI_COMM_WORLD);
#endif

	}


	bool check_do_output(
				double t
			)
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
			simVars->iodata.output_each_sim_seconds < 0 ||
			std::abs(t) < small ||
			std::abs(t - pVars->max_simulation_time) < small ||
			fmod(t, simVars->iodata.output_each_sim_seconds) == 0
		)
			do_output = true;

		return do_output;
	}


	void run()
	{

		std::cout << "Starting run() with mpi_rank " << mpi_rank << std::endl;


#if SWEET_DEBUG
		// DEBUG: full fine simulation
		if (mpi_rank == 0)
		{
			CONSOLEPREFIX_start(0);
			parareal_simulationInstances[0]->run_timestep_fine();
			*(parareal_simulationInstances[0]->parareal_data_fine_exact_debug) = *(parareal_simulationInstances[0]->parareal_data_fine);

			for (int i = 1; i < pVars->coarse_slices; i++)
			{

				CONSOLEPREFIX_start(i - 1);
				Parareal_GenericData &tmp = parareal_simulationInstances[i-1]->get_reference_2_data_timestep_fine();
				Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_2_data_timestep_fine_previous_timestep(); // SL

				// use coarse time step output data as initial data of next coarse time step
				CONSOLEPREFIX_start(i);
				parareal_simulationInstances[i]->sim_set_data(tmp);
				parareal_simulationInstances[i]->sim_set_data_fine_previous_time_slice(tmp2); // SL

				// run coarse time step
				parareal_simulationInstances[i]->run_timestep_fine();

				*(parareal_simulationInstances[i]->parareal_data_fine_exact_debug) = *(parareal_simulationInstances[i]->parareal_data_fine);

			}

		}
#endif

#if SWEET_PARAREAL == 2
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if (mpi_rank == 0)
		{
			// Store initial solution:
			if (pVars->store_iterations)
				parareal_simulationInstances[0]->output_data_file(
						0,
						0,
						true
					);

			CONSOLEPREFIX_start("[MAIN] ");
			std::cout << "Initial propagation" << std::endl;
		}

		/**
		 * Initial propagation
		 */
		if (mpi_rank == 0)
		{
			CONSOLEPREFIX_start(0);
			parareal_simulationInstances[0]->run_timestep_coarse();
			*(parareal_simulationInstances[0]->parareal_data_output) = *(parareal_simulationInstances[0]->parareal_data_coarse);

			for (int i = 1; i < pVars->coarse_slices; i++)
			{
				CONSOLEPREFIX_start(i - 1);
				Parareal_GenericData &tmp = parareal_simulationInstances[i-1]->get_reference_2_data_timestep_coarse();
				Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_2_data_timestep_coarse_previous_timestep(); // SL

				// use coarse time step output data as initial data of next coarse time step
				CONSOLEPREFIX_start(i);
				parareal_simulationInstances[i]->sim_set_data(tmp);
				parareal_simulationInstances[i]->sim_set_data_coarse_previous_time_slice(tmp2); // SL

				// run coarse time step
				parareal_simulationInstances[i]->run_timestep_coarse();

				*(parareal_simulationInstances[i]->parareal_data_output) = *(parareal_simulationInstances[i]->parareal_data_coarse);
			}


			// Store initial propagation:
			if (pVars->store_iterations)
				for (int i = 0; i < pVars->coarse_slices; i++)
					if (timeframe_do_output[i])
						parareal_simulationInstances[i]->output_data_file(
								0,  // 0-th iteration
								i
							);
			// Store initial error relative to reference solution
			if (pVars->load_ref_csv_files)
				for (int i = 0; i < pVars->coarse_slices; i++)
					if (timeframe_do_output[i])
						parareal_simulationInstances[i]->store_parareal_error(
								0,
								i,
								pVars->path_ref_csv_files,
								"ref");
			// Store initial error relative to fine solution
			if (pVars->load_fine_csv_files)
				for (int i = 0; i < pVars->coarse_slices; i++)
					if (timeframe_do_output[i])
						parareal_simulationInstances[i]->store_parareal_error(
								0,
								i,
								pVars->path_fine_csv_files,
								"fine");
		}

#if SWEET_PARAREAL == 2
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		/**
		 * We run as much Parareal iterations as there are coarse slices
		 */
//		int start_slice = 0;

		int k = 0;
		for (; k < pVars->coarse_slices; k++)
		{


			if (mpi_rank == 0)
			{
				CONSOLEPREFIX_start("[MAIN] ");
				std::cout << "Iteration Nr. " << k << std::endl;
			}
			/*
			 * All the following loops should start with 0.
			 * For debugging reasons, we leave it here at 0
			 */


#if SWEET_PARAREAL == 2
			// Communicate coarse solution from proc 0 to proc mpi_rank
			for (int i = k; i < pVars->coarse_slices; i++)
			{
				int working_rank = proc_for_slices[i];
				int local_slice = global_2_local_slice.at(i);

				///Parareal_GenericData* tmp2 = &parareal_simulationInstances[local_slice]->get_reference_2_data_timestep_coarse();
				Parareal_GenericData* tmp2;
				// no need of communication
				if (working_rank == 0 && mpi_rank == 0)
					continue;
				else if (working_rank == mpi_rank) // recv
				{
					tmp2 = parareal_simulationInstances[local_slice]->create_new_data_container("fine");
					communicate_solution(tmp2, 0, mpi_rank, 10000 + i);
					parareal_simulationInstances[local_slice]->sim_set_data(*tmp2);
					delete tmp2;

					tmp2 = parareal_simulationInstances[local_slice]->create_new_data_container("fine");
					communicate_solution(tmp2, 0, mpi_rank, 10000 + 10 + i);
					parareal_simulationInstances[local_slice]->sim_set_data_coarse(*tmp2);
					delete tmp2;
				}
				else if (mpi_rank == 0) // send
				{
					int local_slice_prev = global_2_local_slice.at(i-1);
					tmp2 = &parareal_simulationInstances[local_slice_prev]->get_reference_2_output_data();
					communicate_solution(tmp2, 0, working_rank, 10000 + i);

					tmp2 = &parareal_simulationInstances[local_slice]->get_reference_2_data_timestep_coarse();
					communicate_solution(tmp2, 0, working_rank, 10000 + 10 + i);
				}

			}
#endif

			// SL: send previous timestep for next slice
			for (int i = k; i < pVars->coarse_slices; i++)
			{
				// no previous timestep
				if (i == 0)
					continue;

				// identify proc responsible for this time slice
				int working_rank = proc_for_slices[i];
				int local_slice = global_2_local_slice.at(i);

				// dummy init
				///Parareal_GenericData* tmp2 = &parareal_simulationInstances[local_slice]->get_reference_2_data_timestep_fine_previous_timestep();
				//Parareal_GenericData* tmp2 = parareal_simulationInstances[0]->create_new_data_container();
				Parareal_GenericData* tmp2;

				// there is a previous timestep in this same proc
				if (local_slice > 0 && working_rank == mpi_rank )
				{
					tmp2 = &parareal_simulationInstances[global_2_local_slice.at(i - 1)]->get_reference_2_data_timestep_fine_previous_timestep();
					parareal_simulationInstances[local_slice]->sim_set_data_fine_previous_time_slice(*tmp2);
				}

#if SWEET_PARAREAL == 2
				else if (local_slice == 0 || mpi_rank == 0) // communicate
				{
					if (local_slice == 0) // recv
					{
						///tmp2 = &parareal_simulationInstances[local_slice]->get_reference_2_data_timestep_fine_previous_timestep();
						tmp2 = parareal_simulationInstances[local_slice]->create_new_data_container("fine");
						communicate_solution(tmp2, 0, mpi_rank, 20000 + i);
						parareal_simulationInstances[local_slice]->sim_set_data_fine_previous_time_slice(*tmp2);
						delete tmp2;
					}
					else // send
					{
						tmp2 = &parareal_simulationInstances[global_2_local_slice.at(i - 1)]->get_reference_2_data_timestep_fine_previous_timestep();
						communicate_solution(tmp2, 0, working_rank, 20000 + i);
					}
				}
#endif
			}



			/**
			 * Fine time stepping (in parallel)
			 */
			for (int i = slices_for_proc[0]; i <= slices_for_proc.back(); i++)
			{

				// solution already converged at time slices i < k
				if (i < k)
					continue;

				// identify proc responsible for this time slice
				int working_rank = proc_for_slices[i];
				int local_slice = global_2_local_slice.at(i);

				if (mpi_rank == working_rank)
				{
					CONSOLEPREFIX_start(i);
					parareal_simulationInstances[local_slice]->run_timestep_fine();
				}
			}

#if SWEET_PARAREAL == 2
			MPI_Barrier(MPI_COMM_WORLD);
#endif

			/**
			 * Compute difference between coarse and fine solution
			 */
			for (int i = k; i < pVars->coarse_slices; i++)
			{

				CONSOLEPREFIX_start(i);
				int working_rank = proc_for_slices[i];
				int local_slice = global_2_local_slice.at(i);

				// continue only if this proc is responsible for this time slice
				// or if it is proc 0 (it will receive the computed differences)
				if (local_slice < 0)
					continue;

				// this proc is responsible for this time slice
				if (mpi_rank == working_rank)
					parareal_simulationInstances[local_slice]->compute_difference();

#if SWEET_PARAREAL==2
				// no need to communicate
				if (mpi_rank == 0 && working_rank == 0)
					continue;
				// communicate from each proc to mpi_rank = 0
				///Parareal_GenericData* tmp2 = &parareal_simulationInstances[local_slice]->get_reference_2_data_timestep_diff(); // dummy init
				Parareal_GenericData* tmp2;
				if (mpi_rank == working_rank) // send
				{
					tmp2 = &parareal_simulationInstances[local_slice]->get_reference_2_data_timestep_diff();
					communicate_solution(tmp2, mpi_rank, 0, 30000 + i);
				}
				else if (mpi_rank == 0) // recv
				{
					tmp2 = parareal_simulationInstances[local_slice]->create_new_data_container("fine");
					communicate_solution(tmp2, working_rank, 0, 30000 + i);
					parareal_simulationInstances[local_slice]->sim_set_data_diff(*tmp2);
					delete tmp2;
				}
#endif
			}

			/**
			 * 1) Coarse time stepping (serial)
			 * 2) Compute output + convergence check (serial)
			 * 3) Forward to next frame (serial)
			 */
			if (mpi_rank == 0)
			{
				double max_convergence = -2;
				for (int i = k; i < pVars->coarse_slices; i++)
				{
					CONSOLEPREFIX_start(i);

					if (i > 0)
					{
						// SL:
						Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_2_data_timestep_coarse_previous_timestep();
						parareal_simulationInstances[i]->sim_set_data_coarse_previous_time_slice(tmp2);
					}

					parareal_simulationInstances[i]->run_timestep_coarse();

					// compute convergence
					double convergence = parareal_simulationInstances[i]->compute_output_data(true);
					std::cout << "                        iteration " << k << ", time slice " << i << ", convergence: " << convergence << std::endl;
					if (max_convergence != -1)
						max_convergence = (convergence==-1)?(convergence):(std::max(max_convergence,convergence));


					if (timeframe_do_output[i])
					{
						if (pVars->store_iterations)
							parareal_simulationInstances[i]->output_data_file(
									k + 1,
									i
								);
						if (pVars->load_ref_csv_files)
							parareal_simulationInstances[i]->store_parareal_error(
									k + 1,
									i,
									pVars->path_ref_csv_files,
									"ref");
						if (pVars->load_fine_csv_files)
							parareal_simulationInstances[i]->store_parareal_error(
									k + 1,
									i,
									pVars->path_fine_csv_files,
									"fine");


						CONSOLEPREFIX.start(i);
						parareal_simulationInstances[i]->output_data_console(
								k + 1,
								i
							);

					}

					// last coarse time step slice?
					if (i == pVars->coarse_slices-1)
					{
						// convergence check activated?
						if (pVars->convergence_error_threshold >= 0)
						{
							if (max_convergence >= 0)
							{
								// convergence given?
								if (max_convergence < pVars->convergence_error_threshold)
								{
									CONSOLEPREFIX_start("[MAIN] ");
									std::cout << "Convergence reached at iteration " << k << " with convergence value " << max_convergence << std::endl;
									goto converged;
								}
							}
						}
					}


					// forward to next time slice if it exists
					if (i < pVars->coarse_slices-1)
					{
						CONSOLEPREFIX_start(i);
						Parareal_GenericData &tmp = parareal_simulationInstances[i]->get_reference_2_output_data();

						CONSOLEPREFIX_start(i+1);
						parareal_simulationInstances[i+1]->sim_set_data(tmp);
					}

					parareal_simulationInstances[i]->check_for_nan_parareal();

#if SWEET_DEBUG
					// fine serial solution should be retrieved
					if (i == k)
					{
						std::cout << "Comparing parareal to fine solution at iteration " << i << " and end of timeframe " << i << ". Solutions should be identical." << std::endl;
						parareal_simulationInstances[i]->compare_2_fine_exact();
						std::cout << "Comparison OK!" << std::endl;
					}
#endif

				}
			}

			if (pVars->max_iter >= 0 && k == pVars->max_iter)
				break;
		}

converged:

		CONSOLEPREFIX_end();
	}
};

}}

#endif
