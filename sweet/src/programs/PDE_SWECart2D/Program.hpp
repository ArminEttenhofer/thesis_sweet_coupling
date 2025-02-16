/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_PROGRAMPDESWECART2D_HPP
#define PROGRAMS_PDE_SWECART2D_PROGRAMPDESWECART2D_HPP

#define WITH_TIME_OLD    1

// This is just for the editor to show code as used within precompiler #if ... directives
#include <programs/PDE_SWECart2D/Benchmarks/Shack.hpp>
#include <programs/PDE_SWECart2D/BenchmarksCombined.hpp>
#include <programs/PDE_SWECart2D/NormalModes.hpp>
#include <programs/PDE_SWECart2D/TimeSteppers.hpp>
#include <programs/PDE_SWECart2D/Diagnostics.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/GridMapping.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

// New Time Tree
#include "TimeTree/TimeTreeIR.hpp"

#if SWEET_GUI
#include <sweet/GUI/VisSweet.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

#if SWEET_PARAREAL
#include <sweet/_DEPRECATED_pint/Parareal.hpp>
#endif

#if SWEET_XBRAID
#include <sweet/XBraid/XBraid_sweet_lib.hpp>
#endif

// File writer
#include "FileOutput.hpp"

namespace PDE_SWECart2D {

    class Program
#if SWEET_GUI
        :	public sweet::GUI::SimulationGUICallbacks
#endif
    {
    public:
        sweet::Error::Base error;

        /*
         * Just a class to store simulation data all together
         */
        class DataConfigOps {
        public:
            sweet::Error::Base error;

            sweet::Data::Cart2D::Config cart2DDataConfig;
            sweet::Data::Cart2D::Operators ops;
            sweet::Data::Cart2DComplex::Operators opsComplex;

            DataContainer::Simulation prog;
            DataContainer::Simulation progTmp;

            sweet::Data::Cart2D::DataSpectral t0_prog_h_pert;
            sweet::Data::Cart2D::DataSpectral t0_prog_u;
            sweet::Data::Cart2D::DataSpectral t0_prog_v;

            // Mapping between grids
            sweet::Data::Cart2D::GridMapping gridMapping;

#if SWEET_GUI
            // Data to visualize is stored to this variable
            sweet::Data::Cart2D::DataGrid vis_cart2d_data;

            // Which primitive to use for rendering
            int vis_render_type_of_primitive_id = 1;

            // Which primitive to use for rendering
            int vis_data_id = 0;
#endif

            bool setup(
                    sweet::Data::Cart2D::Shack *i_shackCart2DDataOps,
                    bool i_setup_spectral_transforms = true        // for reset()
            ) {
                /*
                 * Setup Sphere2D Data Config & Operators
                 */
                if (i_setup_spectral_transforms) {
                    cart2DDataConfig.setupAuto(i_shackCart2DDataOps);
                    ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);
                }

                ops.setup(&cart2DDataConfig, i_shackCart2DDataOps);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

                opsComplex.setup(&cart2DDataConfig, i_shackCart2DDataOps);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(opsComplex);

                prog.setup(&cart2DDataConfig);
                progTmp.setup(&cart2DDataConfig);

                t0_prog_h_pert.setup(cart2DDataConfig);
                t0_prog_u.setup(cart2DDataConfig);
                t0_prog_v.setup(cart2DDataConfig);

#if SWEET_GUI
                sweet::Data::Cart2D::Shack shackCart2DDataOps;
#if 0
                // WARNING: We need to use sphere2DDataConfig, since i_shackSphere2DDataOps is not initialized with a reset()
                shackCart2DDataOps.space_res_physical[0] = i_shackSphere2DDataOps->space_res_physical[0];
                shackCart2DDataOps.space_res_physical[1] = i_shackSphere2DDataOps->space_res_physical[1];
#else
                shackCart2DDataOps.space_res_physical[0] = cart2DDataConfig.grid_data_size[0];
                shackCart2DDataOps.space_res_physical[1] = cart2DDataConfig.grid_data_size[1];

#endif
                shackCart2DDataOps.reuse_spectral_transformation_plans = i_shackCart2DDataOps->reuse_spectral_transformation_plans;

                cart2DDataConfig.setupAuto(shackCart2DDataOps);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);
#endif

                if (i_shackCart2DDataOps->space_grid_use_c_staggering)
                    gridMapping.setup(i_shackCart2DDataOps, &cart2DDataConfig);

                return true;
            }

            void clear(bool i_clear_spectral_transforms = true) {
                progTmp.clear();
                prog.clear();

                t0_prog_h_pert.clear();
                t0_prog_u.clear();
                t0_prog_v.clear();

                ops.clear();

                if (i_clear_spectral_transforms)
                    cart2DDataConfig.clear();
            }
        };

        // Simulation data
        DataConfigOps dataConfigOps;

        //////*
        ///// * Just a class to store simulation data all together
        ///// */
        /////class SimDataAndOps
        /////{
        /////public:
        /////	sweet::Error::Base error;

        /////	sweet::Data::Cart2D::Config cart2DDataConfig;
        /////	sweet::Data::Cart2D::Operators ops;

        /////	sweet::Data::Cart2D::DataSpectral prog_h_pert;
        /////	sweet::Data::Cart2D::DataSpectral prog_u;
        /////	sweet::Data::Cart2D::DataSpectral prog_v;

        /////	// TODO: Get rid of me right here
        /////	// Initial values for comparison with analytical solution
        /////	sweet::Data::Cart2D::DataSpectral t0_prog_h_pert;
        /////	sweet::Data::Cart2D::DataSpectral t0_prog_u;
        /////	sweet::Data::Cart2D::DataSpectral t0_prog_v;

        /////	// Forcings
        /////	sweet::Data::Cart2D::DataSpectral force_prog_h_pert;
        /////	sweet::Data::Cart2D::DataSpectral force_prog_u;
        /////	sweet::Data::Cart2D::DataSpectral force_prog_v;

        /////	// Mapping between grids
        /////	sweet::Data::Cart2D::GridMapping gridMapping;
        /////
        /////	// Diagnostics measures
        /////	int last_timestep_nr_update_diagnostics = -1;
        /////
        /////	bool setup(sweet::Data::Cart2D::Shack *i_shackCart2DDataOps)
        /////	{
        /////		/*
        /////		 * Setup Cart2D Data Config & Operators
        /////		 */
        /////		cart2DDataConfig.setupAuto(*i_shackCart2DDataOps);
        /////		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);

        /////		ops.setup(cart2DDataConfig, *i_shackCart2DDataOps);
        /////		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

        /////		prog_h_pert.setup(cart2DDataConfig);
        /////		prog_u.setup(cart2DDataConfig);
        /////		prog_v.setup(cart2DDataConfig);

        /////		t0_prog_h_pert.setup(cart2DDataConfig);
        /////		t0_prog_u.setup(cart2DDataConfig);
        /////		t0_prog_v.setup(cart2DDataConfig);

        /////		force_prog_h_pert.setup(cart2DDataConfig);
        /////		force_prog_u.setup(cart2DDataConfig);
        /////		force_prog_v.setup(cart2DDataConfig);

        /////		last_timestep_nr_update_diagnostics = -1;
        /////
        /////		if (i_shackCart2DDataOps->space_grid_use_c_staggering)
        /////			gridMapping.setup(i_shackCart2DDataOps, &cart2DDataConfig);

        /////		return true;
        /////	}

        /////	void clear()
        /////	{
        /////		prog_h_pert.clear();
        /////		prog_u.clear();
        /////		prog_v.clear();
        /////
        /////		t0_prog_h_pert.clear();
        /////		t0_prog_u.clear();
        /////		t0_prog_v.clear();
        /////
        /////		force_prog_h_pert.clear();
        /////		force_prog_u.clear();
        /////		force_prog_v.clear();

        /////		ops.clear();
        /////		cart2DDataConfig.clear();
        /////	}
        /////};

        /////// Simulation data
        /////SimDataAndOps dataAndOps;


        // time integrators
#if WITH_TIME_OLD
        PDE_SWECart2D::TimeSteppers pdeSWECart2DTimeSteppers;
        ///PDESWESphere2D_TimeSteppers timeSteppers;
#endif
        TimeTree::TimeTree timeSteppersNewTS;
        bool useNewTimeSteppers;

#if SWEET_PARAREAL
        // Implementation of different time steppers
        PDE_SWECart2D::TimeSteppers timeSteppersCoarse;
#endif

        // Handler to all benchmarks
        PDE_SWECart2D::Benchmarks::BenchmarksCombined cart2dBenchmarksCombined;

        FileOutput fileOutput;

        /*
         * Shack directory and shacks to work with
         */
        sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
        sweet::Data::Cart2D::Shack *shackCart2DDataOps;
        sweet::IO::Shack *shackIOData;
        sweet::TimeTree::Shack *shackTimestepControl;
        PDE_SWECart2D::Shack *shackPDESWECart2D;
        PDE_SWECart2D::TimeDiscretization::Shack *shackTimeDisc;
        PDE_SWECart2D::Benchmarks::Shack *shackPDESWECart2DBenchmarks;
        sweet::Parallelization::Shack *shackParallelization;
        ///PDE_SWECart2D::Shack_Diagnostics *shackPDESWECart2DDiagnostics;

#if SWEET_GUI
        // Data to visualize is stored to this variable
        sweet::Data::Cart2D::DataGrid vis_cart2d_data;

        // Which primitive to use for rendering
        int vis_render_type_of_primitive_id = 0;

        // Which primitive to use for rendering
        int vis_data_id = 0;

#endif

        class BenchmarkErrors {
        public:
            //Max difference to initial conditions
            double t0_diff_error_max_abs_h_pert;
            double t0_diff_error_max_abs_u;
            double t0_diff_error_max_abs_v;

            // Error measures L2 norm
            double analytical_error_rms_h;
            double analytical_error_rms_u;
            double analytical_error_rms_v;

            // Error measures max norm
            double analytical_error_maxabs_h;
            double analytical_error_maxabs_u;
            double analytical_error_maxabs_v;

            void setup() {
                t0_diff_error_max_abs_h_pert = -1;
                t0_diff_error_max_abs_u = -1;
                t0_diff_error_max_abs_v = -1;

                analytical_error_rms_h = -1;
                analytical_error_rms_u = -1;
                analytical_error_rms_v = -1;

                analytical_error_maxabs_h = -1;
                analytical_error_maxabs_u = -1;
                analytical_error_maxabs_v = -1;
            }
        };

        BenchmarkErrors benchmarkErrors;

        class NormalModesData {
        public:
            sweet::Error::Base error;

            // Diagnostic information about the projection to
            //    the linear normal wave mode eigenspace (see SWE_bench_normalmodes->hpp)

            sweet::Data::Cart2D::DataSpectral geo;    //Coefficients multiplying geostrophic mode
            sweet::Data::Cart2D::DataSpectral igwest; //Coefficients multiplying west gravity mode
            sweet::Data::Cart2D::DataSpectral igeast; //Coefficients multiplying east gravity mode
            double norm_spec;

            PDE_SWECart2D::NormalModes::NormalModes pdeSWECart2DNormalModes;

        public:
            bool shackRegistration(sweet::Shacks::Dictionary &io_dict) {
                pdeSWECart2DNormalModes.shackRegistration(io_dict);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWECart2DNormalModes);
                return true;
            }

        public:
            bool setup(
                    sweet::Data::Cart2D::Config *cart2DDataConfig
            ) {
                geo.setup(cart2DDataConfig);
                igwest.setup(cart2DDataConfig);
                igeast.setup(cart2DDataConfig);

                return true;
            }

        public:
            bool clear() {
                geo.clear();
                igwest.clear();
                igeast.clear();
                return true;
            }
        };

        NormalModesData *normalmodes;

        //! Diagnostic measures at initial stage, Initialize with 0
        double diagnostics_energy_start = 0;
        double diagnostics_mass_start = 0;
        double diagnostics_potential_enstrophy_start = 0;


        bool compute_error_difference_2_initial_condition = false;
        bool compute_error_2_analytical_solution = false;
        bool compute_normal_modes = false;

        Diagnostics diagnostics;

    public:
        Program(
                int i_argc,
                char *const *const i_argv
        ) :
                useNewTimeSteppers(true),
                shackProgArgDict(i_argc, i_argv),
                shackCart2DDataOps(nullptr),
                shackIOData(nullptr),
                shackTimestepControl(nullptr),
                shackPDESWECart2D(nullptr),
                shackTimeDisc(nullptr),
                shackPDESWECart2DBenchmarks(nullptr),
                shackParallelization(nullptr)
        ///shackPDESWECart2DDiagnostics(nullptr)
        {
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
        }


        bool setup_1_registration() {
            /*
             * Setup argument parsing
             */
            shackProgArgDict.setup();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            /*
             * SHACK: Register classes which we require
             */
            shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
            shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
            shackPDESWECart2D = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::Shack>();
            shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
            shackTimeDisc = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::TimeDiscretization::Shack>();
            shackPDESWECart2DBenchmarks = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack>();
            shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
            ///shackPDESWECart2DDiagnostics = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::Shack_Diagnostics>();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            /*
             * SHACK: Register other things before parsing program arguments
             */
            cart2dBenchmarksCombined.shackRegistration(&shackProgArgDict);
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2dBenchmarksCombined);

            shackProgArgDict.processProgramArgumentsForShack(shackTimeDisc);

            if (shackTimeDisc->timestepping_method.find('(') != std::string::npos)
                useNewTimeSteppers = true;
            else
                useNewTimeSteppers = false;

            if (useNewTimeSteppers) {
                /*
                 * Setup NEW time steppers
                 */
                timeSteppersNewTS.setup_1_registerAllTimesteppers();
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);

                timeSteppersNewTS.setup_2_shackRegistration(&shackProgArgDict);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);

            } else {
#if WITH_TIME_OLD

                /*
                 * Setup legacy time steppers
                 */
                pdeSWECart2DTimeSteppers.shackRegistration(shackProgArgDict);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWECart2DTimeSteppers);
                ////pdeSWECart2DTimeSteppers.setup_1_registerAllTimesteppers();
                ////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWECart2DTimeSteppers);

                ////pdeSWECart2DTimeSteppers.setup_2_shackRegistration(&shackProgArgDict);
                ////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWECart2DTimeSteppers);
#endif
            }


            if (shackPDESWECart2D->normal_mode_analysis_generation) {
                normalmodes = new NormalModesData;
                normalmodes->shackRegistration(shackProgArgDict);
            }

            /////if (!useNewTimeSteppers)
            /////{
            /////	/*
            /////	 * Close shack registration & getting shacks
            /////	 */
            /////	shackProgArgDict.closeRegistration();

            /////	shackProgArgDict.closeGet();
            /////}

            return true;
        }

        void clear_1_shackRegistration() {
            if (shackPDESWECart2D->normal_mode_analysis_generation) {
                delete normalmodes;
                normalmodes = nullptr;
            }

            shackCart2DDataOps = nullptr;
            shackTimestepControl = nullptr;
            shackIOData = nullptr;
            shackTimeDisc = nullptr;
            shackParallelization = nullptr;

            cart2dBenchmarksCombined.clear();
#if WITH_TIME_OLD
            pdeSWECart2DTimeSteppers.clear();
#endif
            shackProgArgDict.clear();
        }

        bool setup_2_processArguments() {
            ///shackProgArgDict.setup();

            /*
             * SHACK: Process arguments
             */
            shackProgArgDict.processProgramArguments();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            shackProgArgDict.processHelpArguments();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            /*
             * Do some validation of program arguments
             */
            shackTimestepControl->validateTimestepSize();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

            return true;
        }

        void clear_2_process_arguments() {
            shackProgArgDict.clear();
        }


        bool setup_3_main() {
            /*
             * Setup Cart2D Data Config & Operators
             */
            dataConfigOps.setup(shackCart2DDataOps);
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

            /*
             * Now we're ready to setup the time steppers
             */
            if (useNewTimeSteppers) {
                timeSteppersNewTS.setup_3_timestepper(
                        shackTimeDisc->timestepping_method,
                        &shackProgArgDict,
                        &dataConfigOps.ops,
                        &dataConfigOps.opsComplex,
                        dataConfigOps.prog
                );
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);

                timeSteppersNewTS.timeIntegrator->setTimeStepSize(shackTimestepControl->currentTimestepSize);

                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);
            } else {
#if WITH_TIME_OLD
                /*
                 * After we setup the cart2d, we can setup the time steppers and their buffers
                 */
                pdeSWECart2DTimeSteppers.setup(&dataConfigOps.ops, &shackProgArgDict);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWECart2DTimeSteppers);

                ////pdeSWECart2DTimeSteppers.setup_3_timestepper(
                ////		shackTimeDisc->timestepping_method,
                ////		&shackProgArgDict,
                ////		&dataConfigOps.ops
                ////	);
                ////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);
#endif
            }


#if 0
            pdeSWECart2DTimeSteppers.timestepper->runTimestep(
                    dataConfigOps.prog.h_pert, dataConfigOps.prog.u, dataConfigOps.prog.v,
                    shackTimestepControl->currentTimestepSize,
                    shackTimestepControl->currentSimulationTime
                );
#endif

#if SWEET_GUI
            vis_cart2d_data.setup(dataConfigOps.cart2DDataConfig);
#endif

            std::cout << "Printing shack information:" << std::endl;
            shackProgArgDict.printShackData();

            cart2dBenchmarksCombined.setupInitialConditions(
                    dataConfigOps.prog.h_pert,
                    dataConfigOps.prog.u,
                    dataConfigOps.prog.v,
                    &dataConfigOps.ops,
                    &dataConfigOps.cart2DDataConfig
            );
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2dBenchmarksCombined);

            dataConfigOps.t0_prog_h_pert = dataConfigOps.prog.h_pert;

            /*
             * Finish registration & getting class interfaces so that nobody can do some
             * strange things with this anymore
             */
            shackProgArgDict.closeRegistration();
            shackProgArgDict.closeGet();

            benchmarkErrors.setup();

            /*
             * Now we should check that all program arguments have really been parsed
             */
            // shackProgArgDict.checkAllArgumentsProcessed();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            if (shackPDESWECart2D->normal_mode_analysis_generation) {
                normalmodes->setup(&dataConfigOps.cart2DDataConfig);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*normalmodes);
            }

            if (shackPDESWECart2D->compute_errors) {
                //Compute difference to initial condition (makes more sense in steady state cases, but useful in others too)
                compute_error_difference_2_initial_condition = true;

#if WITH_TIME_OLD
                //Compute difference to analytical solution (makes more sense in linear cases, but might be useful in others too)
                compute_error_2_analytical_solution = pdeSWECart2DTimeSteppers.linear_only;
#endif
            } else {
                compute_error_difference_2_initial_condition = false;
                compute_error_2_analytical_solution = false;
            }


            fileOutput.setup(shackIOData, shackTimestepControl, shackCart2DDataOps, shackPDESWECart2D,
                             &dataConfigOps.cart2DDataConfig, &dataConfigOps.ops, &dataConfigOps.opsComplex);

            if (shackPDESWECart2DBenchmarks->benchmark_name == "normalmodes")
                compute_normal_modes = true;

            if (compute_normal_modes) {
                update_normal_modes();
                update_diagnostics();
            }

            diagnostics_energy_start = diagnostics.total_energy;
            diagnostics_mass_start = diagnostics.total_mass;
            diagnostics_potential_enstrophy_start = diagnostics.total_potential_enstrophy;

            if (useNewTimeSteppers) {
                shackProgArgDict.closeRegistration();
                shackProgArgDict.closeGet();
            }

            return true;
        }


        void clear_3_main() {
            if (shackPDESWECart2D->normal_mode_analysis_generation) {
                normalmodes->clear();
                delete normalmodes;
                normalmodes = nullptr;
            }

            fileOutput.clear();

#if SWEET_GUI
            vis_cart2d_data.clear();
#endif

            if (useNewTimeSteppers) {
                timeSteppersNewTS.clear();
            } else {
#if WITH_TIME_OLD
                pdeSWECart2DTimeSteppers.clear();
#endif
            }


            dataConfigOps.clear();
        }


        bool setup() {
            if (!setup_1_registration())
                return false;

            if (!setup_2_processArguments())
                return false;

            if (!setup_3_main())
                return false;

            std::cout << "SETUP FINISHED" << std::endl;
            return true;
        }

        void clear() {
            clear_3_main();
            clear_2_process_arguments();
            clear_1_shackRegistration();
        }

        bool reset() {
            // keep pause state if in GUI mode
            bool run_simulation_timesteps = shackTimestepControl->runSimulationTimesteps;

            clear();

            if (!setup()) {
                error.print();
                return false;
            }

            shackTimestepControl->runSimulationTimesteps = run_simulation_timesteps;

            return !error.exists();
        }

        // Update diagnostic variables related to normal modes
        void update_normal_modes() {
            if (!compute_normal_modes)
                return;

#if SWEET_USE_CART2D_SPECTRAL_SPACE
            // Setup cart2dDiagnostics for normal mode projection
            normalmodes->pdeSWECart2DNormalModes.convert_allspectralmodes_2_normalmodes(
                    dataConfigOps.prog.h_pert, dataConfigOps.prog.u, dataConfigOps.prog.v,
                    normalmodes->geo, normalmodes->igwest, normalmodes->igeast
            );

            if (shackTimestepControl->currentTimestepNr == 0) {
                //save the reference normalization parameter
                std::cout << normalmodes->geo.spectral_reduce_rms() << std::endl;
                std::cout << normalmodes->igwest.spectral_reduce_rms() << std::endl;
                std::cout << normalmodes->igeast.spectral_reduce_rms() << std::endl;

                normalmodes->norm_spec = normalmodes->geo.spectral_reduce_sum_sqr_quad() +
                                         normalmodes->igwest.spectral_reduce_sum_sqr_quad() +
                                         normalmodes->igeast.spectral_reduce_sum_sqr_quad();

                normalmodes->norm_spec = std::sqrt(normalmodes->norm_spec);

                if (normalmodes->norm_spec < 10e-14) {
                    normalmodes->norm_spec = 1.0;
                    return;
                }

            }
            normalmodes->geo = normalmodes->geo / normalmodes->norm_spec;
            normalmodes->igwest = normalmodes->igwest / normalmodes->norm_spec;
            normalmodes->igeast = normalmodes->igeast / normalmodes->norm_spec;
#endif
        }


        //Update diagnostic variables related to normal modes
        void dump_normal_modes() {
            if (!compute_normal_modes)
                return;

#if SWEET_USE_CART2D_SPECTRAL_SPACE
            PDE_SWECart2D::Benchmarks::normal_modes n;
            n.dump_all_normal_modes(normalmodes->geo, normalmodes->igwest, normalmodes->igeast);
#endif

            return;
        }

        //Calculate the model diagnostics
        void update_diagnostics() {
            // assure, that the cart2dDiagnostics are only updated for new time steps
            if (diagnostics.last_update_timestep_nr == shackTimestepControl->currentTimestepNr)
                return;

            diagnostics.last_update_timestep_nr = shackTimestepControl->currentTimestepNr;


            if (shackCart2DDataOps->space_grid_use_c_staggering) {
                diagnostics.update_staggered_huv_2_mass_energy_enstrophy(
                        dataConfigOps.ops,
                        shackCart2DDataOps,
                        shackPDESWECart2D,
                        dataConfigOps.prog.h_pert,
                        dataConfigOps.prog.u,
                        dataConfigOps.prog.v
                );
            } else {
                diagnostics.update_nonstaggered_huv_2_mass_energy_enstrophy(
                        dataConfigOps.ops,
                        shackCart2DDataOps,
                        shackPDESWECart2D,
                        dataConfigOps.prog.h_pert,
                        dataConfigOps.prog.u,
                        dataConfigOps.prog.v
                );
            }
        }

        void normal_mode_analysis() {
            normalmodes->pdeSWECart2DNormalModes.normal_mode_analysis(
                    dataConfigOps.prog.h_pert,
                    dataConfigOps.prog.u,
                    dataConfigOps.prog.v,
                    3,
                    &shackProgArgDict,
                    this,
                    &Program::runTimestep
            );

            std::cout << "\n Done normal mode analysis in separate class" << std::endl;
        }

        double precice_max_dt = -1;

        /**
         * Execute a single simulation time step
         */
        bool runTimestep() {
            shackTimestepControl->timestepHelperStart();

            ////std::cout << "AAA " << shackTimestepControl->currentSimulationTime << " " << dataConfigOps.prog.h_pert.spectral_reduce_max_abs() << " " << dataConfigOps.prog.h_pert.toGrid().grid_reduce_max_abs() << std::endl;

            if (shackTimestepControl->timestepHelperStart()) {
                // update time step size if it's changed!
                if (useNewTimeSteppers)
                    timeSteppersNewTS.timeIntegrator->setTimeStepSize(shackTimestepControl->currentTimestepSize);
            }

            // Check the timestep requirements from precice
            if (precice_max_dt >= 0) {
                // Adjust the timestep depending on the coupling
                double minDt = 10e-14;
                if (precice_max_dt - shackTimestepControl->currentTimestepSize < minDt) {
                    // The next time step would be too small or surpass the end
                    shackTimestepControl->currentTimestepSize = precice_max_dt;

                    if (useNewTimeSteppers)
                        timeSteppersNewTS.timeIntegrator->setTimeStepSize(shackTimestepControl->currentTimestepSize);
                }
            }

            if (useNewTimeSteppers) {
                timeSteppersNewTS.runIntegration(
                        dataConfigOps.prog,
                        dataConfigOps.progTmp,
                        shackTimestepControl->currentSimulationTime
                );
                dataConfigOps.prog.swap(dataConfigOps.progTmp);
            } else {
#if WITH_TIME_OLD
                pdeSWECart2DTimeSteppers.timestepper->runTimestep(
                        dataConfigOps.prog.h_pert, dataConfigOps.prog.u, dataConfigOps.prog.v,
                        shackTimestepControl->currentTimestepSize,
                        shackTimestepControl->currentSimulationTime
                );
#endif
            }

            shackTimestepControl->timestepHelperEnd();

#if !SWEET_PARAREAL
            timestep_do_output();
#endif

            if (shackPDESWECart2D->compute_diagnostics)
                update_diagnostics();

            return true;
        }


//////	/**
//////	 * Write file to data and return string of file name
//////	 */
//////	std::string write_file(
//////			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
//////			const char* i_name	//!< name of output variable
//////		)
//////	{
//////		char buffer[1024];
//////
//////		// TODO: convert spectral datato physical
//////
//////		const char* filename_template = shackIOData->outputFileName.c_str();
//////		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);
//////		i_cart2DData.toGrid().file_grid_saveData_ascii(buffer);
//////		return buffer;
//////	}
//////
//////	/**
//////	 * Write current time step info to file
//////	 */
//////
//////	std::string write_output_file(
//////			std::stringstream &buffer
//////		)
//////	{
//////		const char* filename_template = "output_diag_evol.txt";
//////		std::ofstream file(filename_template, std::ofstream::out | std::ofstream::app);
//////		file << std::setprecision(12);
//////  		file << buffer.str() << std::endl;
//////
//////		return buffer.str();
//////	}
//////
//////
//////	/**
//////	 * Write spectrum info to data and return string of file name
//////	 */
//////#if SWEET_USE_CART2D_SPECTRAL_SPACE
//////	std::string write_file_spec(
//////			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
//////			const char* i_name	//!< name of output variable
//////		)
//////	{
//////		char buffer[1024];
//////
//////		const char* filename_template = shackIOData->outputFileName.c_str();
//////		sprintf(buffer, filename_template, i_name, shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScale);
//////		i_cart2DData.file_spectral_abs_saveData_ascii(buffer);
//////		//i_cart2DData.file_spectral_saveData_ascii(buffer);
//////		return buffer;
//////	}
//////#endif

        std::string output_filenames;

    public:
        bool timestep_do_output(
                std::ostream &o_ostream = std::cout
        ) {
            if (shackPDESWECart2D->normal_mode_analysis_generation > 0)
                return false;

            if (!shackIOData->checkDoOutput(shackTimestepControl->currentSimulationTime))
                return false;

            if (shackParallelization->isMPIRoot) {
                fileOutput.write_file_output(
                        dataConfigOps.ops,
                        dataConfigOps.prog.h_pert,
                        dataConfigOps.prog.u,
                        dataConfigOps.prog.v
                );

                if (shackIOData->verbosity > 0) {
                    double progHMin = dataConfigOps.prog.h_pert.toGrid().grid_reduce_min();
                    double progHMax = dataConfigOps.prog.h_pert.toGrid().grid_reduce_max();

                    if (shackParallelization->isMPIRoot) {
                        std::cout << "prog_h min/max:\t" << progHMin << ", " << progHMax << std::endl;
                    }
                }
            }


/////////////		/*
/////////////		 * File output
/////////////		 *
/////////////		 * We write everything in non-staggered output
/////////////		 */
/////////////		// For output, variables need to be on unstaggered A-grid
/////////////		sweet::Data::Cart2D::DataGrid t_h(dataConfigOps.cart2DDataConfig);
/////////////		sweet::Data::Cart2D::DataGrid t_u(dataConfigOps.cart2DDataConfig);
/////////////		sweet::Data::Cart2D::DataGrid t_v(dataConfigOps.cart2DDataConfig);
/////////////
/////////////		if (shackCart2DDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
/////////////		{
/////////////			t_h = dataConfigOps.prog.h_pert.toGrid();
/////////////			dataConfigOps.gridMapping.mapCtoA_u(dataConfigOps.prog.u.toGrid(), t_u);
/////////////			dataConfigOps.gridMapping.mapCtoA_v(dataConfigOps.prog.v.toGrid(), t_v);
/////////////		}
/////////////		else
/////////////		{
/////////////			t_h = dataConfigOps.prog.h_pert.toGrid();
/////////////			t_u = dataConfigOps.prog.u.toGrid();
/////////////			t_v = dataConfigOps.prog.v.toGrid();
/////////////		}
/////////////
/////////////
/////////////		// Dump  data in csv, if output filename is not empty
/////////////		if (shackIOData->outputFileName.size() > 0)
/////////////		{
/////////////			output_filenames = "";
/////////////
/////////////			output_filenames = write_file(t_h, "prog_h_pert");
/////////////			output_filenames += ";" + write_file(t_u, "prog_u");
/////////////			output_filenames += ";" + write_file(t_v, "prog_v");
/////////////
/////////////			output_filenames += ";" + write_file(dataConfigOps.ops.ke(t_u,t_v),"diag_ke");
/////////////
/////////////#if SWEET_USE_CART2D_SPECTRAL_SPACE
/////////////			output_filenames += ";" + write_file_spec(dataConfigOps.ops.ke(t_u,t_v),"diag_ke_spec");
/////////////
/////////////			output_filenames += ";" + write_file_spec(t_h, "prog_h_pert_spec");
/////////////			output_filenames += ";" + write_file_spec(t_u, "prog_u_spec");
/////////////			output_filenames += ";" + write_file_spec(t_v, "prog_v_spec");
/////////////
/////////////			output_filenames += ";" + write_file_spec(dataConfigOps.ops.ke(t_u,t_v).toGrid(), "diag_ke_spec");
/////////////#endif
/////////////
/////////////			output_filenames += ";" + write_file(dataConfigOps.ops.vort(t_u, t_v), "diag_vort");
/////////////			output_filenames += ";" + write_file(dataConfigOps.ops.div(t_u, t_v), "diag_div");
/////////////
/////////////#if SWEET_USE_CART2D_SPECTRAL_SPACE
/////////////			if (compute_normal_modes){
/////////////				output_filenames += ";" + write_file_spec(normalmodes->geo, "nm_geo");
/////////////				output_filenames += ";" + write_file_spec(normalmodes->igwest, "nm_igwest");
/////////////				output_filenames += ";" + write_file_spec(normalmodes->igeast, "nm_igeast");
/////////////			}
/////////////#endif
/////////////		}

            update_normal_modes();
#if SWEET_USE_CART2D_SPECTRAL_SPACE
            if (shackIOData->outputFileName.size() > 0) {
                std::string output_filenames = "";
                if (compute_normal_modes) {
                    output_filenames += ";" + fileOutput.write_file_csv_spec_evol(normalmodes->geo, "nm_geo");
                    output_filenames += ";" + fileOutput.write_file_csv_spec_evol(normalmodes->igwest, "nm_igwest");
                    output_filenames += ";" + fileOutput.write_file_csv_spec_evol(normalmodes->igeast, "nm_igeast");
                }
            }
#endif


            if (shackIOData->verbosity > 0) {
                update_diagnostics();
                computeErrors();

                std::stringstream header;
                std::stringstream rows;

                rows << std::setprecision(16);

                // Prefix
                if (shackTimestepControl->currentTimestepNr == 0)
                    header << "DATA";
                rows << "DATA";

                // Time
                if (shackTimestepControl->currentTimestepNr == 0)
                    header << "\tT";
                rows << "\t" << shackTimestepControl->currentSimulationTime;

#if 1
                // Mass, Energy, Enstrophy
                header << "\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";
                rows << "\t" << diagnostics.total_mass;
                rows << "\t" << diagnostics.total_energy;
                rows << "\t" << diagnostics.total_potential_enstrophy;

                // Mass, Energy, Enstrophy
                header << "\tTOTAL_MASS_REL_ERROR\tTOTAL_ENERGY_REL_ERROR\tPOT_ENSTROPHY_REL_ERROR";
                rows << "\t" << std::abs((diagnostics.total_mass - diagnostics_mass_start) / diagnostics_mass_start);
                rows << "\t"
                     << std::abs((diagnostics.total_energy - diagnostics_energy_start) / diagnostics_energy_start);
                rows << "\t" << std::abs(
                        (diagnostics.total_potential_enstrophy - diagnostics_potential_enstrophy_start) /
                        diagnostics_potential_enstrophy_start);
#endif

#if 1
                if (compute_error_difference_2_initial_condition) {
                    // Difference to initial condition
                    if (shackTimestepControl->currentTimestepNr == 0)
                        header << "\tDIFF_MAXABS_H0\tDIFF_MAXABS_U0\tDIFF_MAXABS_V0";

                    rows << "\t" << benchmarkErrors.t0_diff_error_max_abs_h_pert << "\t"
                         << benchmarkErrors.t0_diff_error_max_abs_u << "\t" << benchmarkErrors.t0_diff_error_max_abs_v;
                }
#endif

#if 1
                if (compute_error_2_analytical_solution) {
                    if (shackTimestepControl->currentTimestepNr == 0)
                        header << "\tREF_DIFF_MAX_H\tREF_DIFF_MAX_U\tREF_DIFF_MAX_V";

                    rows << "\t" << benchmarkErrors.analytical_error_maxabs_h << "\t"
                         << benchmarkErrors.analytical_error_maxabs_u << "\t"
                         << benchmarkErrors.analytical_error_maxabs_v;
                }
#endif

#if 1
                // Normal mode stuff
                if (compute_normal_modes) {
                    // normal modes energy
                    if (shackTimestepControl->currentTimestepNr == 0)
                        header << "\tNM_GEO_RMS\tNM_IGWEST_RMS\tNM_IGEAST_RMS";

                    rows << "\t" << normalmodes->geo.spectral_reduce_rms();
                    rows << "\t" << normalmodes->igwest.spectral_reduce_rms();
                    rows << "\t" << normalmodes->igeast.spectral_reduce_rms();

                    //Dump to file all normal mode evolution
                    dump_normal_modes();
                }
#endif

                // screen output
                if (shackTimestepControl->currentTimestepNr == 0)
                    o_ostream << header.str() << std::endl;

                o_ostream << rows.str() << std::endl;

#if 1
                // output to file
                if (shackTimestepControl->currentTimestepNr == 0)
                    fileOutput.write_output_file(header);
                fileOutput.write_output_file(rows);
#endif

#if 1
                if (diagnostics_mass_start > 0.00001 &&
                    std::abs((diagnostics.total_mass - diagnostics_mass_start) / diagnostics_mass_start) > 10000000.0) {
                    std::cerr << "\n DIAGNOSTICS MASS DIFF TOO LARGE:\t"
                              << std::abs((diagnostics.total_mass - diagnostics_mass_start) / diagnostics_mass_start)
                              << std::endl;
                }
#endif

            }

            shackIOData->advanceNextOutput(shackTimestepControl->currentSimulationTime,
                                           shackTimestepControl->maxSimulationTime);
            return true;
        }


    public:
        void computeErrors() {
#if WITH_TIME_OLD
            if (compute_error_difference_2_initial_condition || compute_error_2_analytical_solution) {
                /**
                 * Compute difference to initial condition
                 */
                if (compute_error_difference_2_initial_condition) {
                    benchmarkErrors.t0_diff_error_max_abs_h_pert = (dataConfigOps.prog.h_pert -
                                                                    dataConfigOps.t0_prog_h_pert).toGrid().grid_reduce_max_abs();
                    benchmarkErrors.t0_diff_error_max_abs_u = (dataConfigOps.prog.u -
                                                               dataConfigOps.t0_prog_u).toGrid().grid_reduce_max_abs();
                    benchmarkErrors.t0_diff_error_max_abs_v = (dataConfigOps.prog.v -
                                                               dataConfigOps.t0_prog_v).toGrid().grid_reduce_max_abs();
                }

                // Calculate linear exact solution, if compute error requests
                if (compute_error_2_analytical_solution) {
                    // Analytical solution at specific time on A-grid
                    sweet::Data::Cart2D::DataSpectral ts_h_pert = dataConfigOps.t0_prog_h_pert;
                    sweet::Data::Cart2D::DataSpectral ts_u = dataConfigOps.t0_prog_u;
                    sweet::Data::Cart2D::DataSpectral ts_v = dataConfigOps.t0_prog_v;

                    // Run exact solution for linear case
                    if (pdeSWECart2DTimeSteppers.l_direct == nullptr) {
                        std::cerr << "Direct solution not available" << std::endl;
                        exit(1);
                    }

                    pdeSWECart2DTimeSteppers.l_direct->runTimestep(
                            ts_h_pert, ts_u, ts_v,
                            shackTimestepControl->currentSimulationTime,    // time step size
                            0                // initial condition given at time 0
                    );

                    benchmarkErrors.analytical_error_rms_h = (ts_h_pert -
                                                              dataConfigOps.prog.h_pert).toGrid().grid_reduce_rms();
                    benchmarkErrors.analytical_error_rms_u = (ts_u - dataConfigOps.prog.u).toGrid().grid_reduce_rms();
                    benchmarkErrors.analytical_error_rms_v = (ts_v - dataConfigOps.prog.v).toGrid().grid_reduce_rms();

                    benchmarkErrors.analytical_error_maxabs_h = (ts_h_pert -
                                                                 dataConfigOps.prog.h_pert).toGrid().grid_reduce_max_abs();
                    benchmarkErrors.analytical_error_maxabs_u = (ts_u -
                                                                 dataConfigOps.prog.u).toGrid().grid_reduce_max_abs();
                    benchmarkErrors.analytical_error_maxabs_v = (ts_v -
                                                                 dataConfigOps.prog.v).toGrid().grid_reduce_max_abs();
                }
            }
#endif
        }

    public:
        void printErrors() {

            if (1) {
                update_diagnostics();

                std::cout << "DIAGNOSTICS ENERGY DIFF:\t"
                          << std::abs((diagnostics.total_energy - diagnostics_energy_start) / diagnostics_energy_start)
                          << std::endl;
                std::cout << "DIAGNOSTICS MASS DIFF:\t"
                          << std::abs((diagnostics.total_mass - diagnostics_mass_start) / diagnostics_mass_start)
                          << std::endl;
                std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs(
                        (diagnostics.total_potential_enstrophy - diagnostics_potential_enstrophy_start) /
                        diagnostics_potential_enstrophy_start) << std::endl;

            }

            if (shackPDESWECart2D->compute_errors) {
                std::cout << "BENCHMARK DIFF H(t=0):\t" << benchmarkErrors.t0_diff_error_max_abs_h_pert << std::endl;
                std::cout << "BENCHMARK DIFF U(t=0):\t" << benchmarkErrors.t0_diff_error_max_abs_u << std::endl;
                std::cout << "BENCHMARK DIFF V(t=0):\t" << benchmarkErrors.t0_diff_error_max_abs_v << std::endl;

                std::cout << "[MULE] error_end_linf_h_pert: " << benchmarkErrors.t0_diff_error_max_abs_h_pert
                          << std::endl;
                std::cout << "[MULE] error_end_linf_u: " << benchmarkErrors.t0_diff_error_max_abs_u << std::endl;
                std::cout << "[MULE] error_end_linf_v: " << benchmarkErrors.t0_diff_error_max_abs_v << std::endl;
                std::cout << std::endl;

                std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << benchmarkErrors.analytical_error_rms_h << std::endl;
                std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << benchmarkErrors.analytical_error_rms_u << std::endl;
                std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << benchmarkErrors.analytical_error_rms_v << std::endl;

                std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << benchmarkErrors.analytical_error_maxabs_h
                          << std::endl;
                std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << benchmarkErrors.analytical_error_maxabs_u
                          << std::endl;
                std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << benchmarkErrors.analytical_error_maxabs_v
                          << std::endl;
            }
        }


    public:
        bool should_quit() {
            if (shackTimestepControl->maxTimestepsNr != -1)
                if (shackTimestepControl->maxTimestepsNr <= shackTimestepControl->currentTimestepNr)
                    return true;

            if (shackTimestepControl->maxSimulationTime != -1)
                if (shackTimestepControl->maxSimulationTime <= shackTimestepControl->currentSimulationTime +
                                                               shackTimestepControl->maxSimulationTime *
                                                               1e-10)    // care about roundoff errors with 1e-10
                    return true;

            return false;
        }


#if SWEET_GUI

        /**
         * postprocessing of frame: do time stepping
         */
        void vis_post_frame_processing(
                int i_num_iterations
        )
        {
            if (shackTimestepControl->runSimulationTimesteps)
                for (int i = 0; i < i_num_iterations && !should_quit(); i++)
                    runTimestep();
        }


        struct VisStuff
        {
            const sweet::Data::Cart2D::DataSpectral* data;
            const char *description;
        };

        /**
         * Arrays for online visualisation and their textual description
         */
        VisStuff vis_arrays[3] =
        {
                {&dataConfigOps.prog.h_pert,	"h"},
                {&dataConfigOps.prog.u,	"u"},
                {&dataConfigOps.prog.v,	"v"},
        };


        void vis_getDataArray(
                const sweet::Data::Cart2D::DataGrid **o_dataArray,
                double *o_aspect_ratio,
                int *o_render_primitive,
                void **o_bogus_data,
                double *o_viz_min,
                double *o_viz_max,
                bool *viz_reset
        )
        {
            if (vis_data_id < 0)
            {
                // Analytical solution at specific time on A-grid
                sweet::Data::Cart2D::DataSpectral ts_h_pert = dataConfigOps.t0_prog_h_pert;
                sweet::Data::Cart2D::DataSpectral ts_u = dataConfigOps.t0_prog_u;
                sweet::Data::Cart2D::DataSpectral ts_v = dataConfigOps.t0_prog_v;

                if (vis_data_id == -1 || vis_data_id == -2 )
                {
                    // Missing to setup REXIFunctions, so invalid phi function set

                    if (pdeSWECart2DTimeSteppers.l_direct)
                    {
                        // Run exact solution for linear case
                        pdeSWECart2DTimeSteppers.l_direct->runTimestep(
                                ts_h_pert, ts_u, ts_v,
                                shackTimestepControl->currentSimulationTime,
                                0			// initial condition given at time 0
                        );
                    }
                    else
                    {
                        SWEETErrorFatal("Direct solution not available");
                    }
                }

                sweet::Data::Cart2D::DataSpectral vis_tmp;
                switch(vis_data_id)
                {
                case -1:
                    vis_tmp = ts_h_pert+shackPDESWECart2D->h0;			// Exact linear solution
                    break;

                case -2:
                    vis_tmp = ts_h_pert-dataConfigOps.prog.h_pert;			// difference to exact linear solution
                    break;

                case -3:
                    vis_tmp = dataConfigOps.t0_prog_h_pert-dataConfigOps.prog.h_pert;	// difference to initial condition
                    break;

                case -4:
                    vis_tmp = dataConfigOps.ops.diff_c_x(dataConfigOps.prog.v) - dataConfigOps.ops.diff_c_y(dataConfigOps.prog.u);	// relative vorticity
                    break;
                case -5:
                    vis_tmp = dataConfigOps.prog.h_pert;			//Perturbation of depth
                    break;
                case -6:
                    vis_tmp = normalmodes->geo ;	// geostrophic mode
                    break;
                case -7:
                    vis_tmp = normalmodes->igwest;	// inertia grav mode west
                    break;
                case -8:
                    vis_tmp = normalmodes->igeast;	// inertia grav mode east
                    break;
                }

                vis_cart2d_data = vis_tmp.toGrid();
                *o_dataArray = &vis_cart2d_data;
                *o_aspect_ratio = shackCart2DDataOps->cart2d_domain_size[1] / shackCart2DDataOps->cart2d_domain_size[0];
                return;
            }

            int id = vis_data_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

            sweet::Data::Cart2D::DataSpectral vis_tmp;
            if (id == 0)
            {
                vis_tmp = *vis_arrays[id].data+shackPDESWECart2D->h0;
            }
            else
            {
                vis_tmp = *(vis_arrays[id].data);
            }

            vis_cart2d_data = vis_tmp.toGrid();
            *o_dataArray = &vis_cart2d_data;

            *o_aspect_ratio = shackCart2DDataOps->cart2d_domain_size[1] / shackCart2DDataOps->cart2d_domain_size[0];
        }



        /**
         * return status string for window title
         */
        const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
        {
            // first, update cart2dDiagnostics if required
            update_normal_modes();
            update_diagnostics();

            const char* description = "";
            if (vis_data_id >= 0)
            {
                int id = vis_data_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
                description = vis_arrays[id].description;
            }
            else
            {
                switch (vis_data_id)
                {
                case -1:
                    description = "Direct solution for h (linear only)";
                    break;

                case -2:
                    description = "Diff in h to exact linear solution";
                    break;

                case -3:
                    description = "Diff in h to initial condition";
                    break;
                case -4:
                    description = "Relative vorticity";
                    break;
                case -5:
                    description = "Depth perturbation";
                    break;
                case -6:
                    description = "Geostrophic wave";
                    break;
                case -7:
                    description = "Inertia-gravity west wave";
                    break;
                case -8:
                    description = "Inertia-gravity east wave";
                    break;
                }
            }

            static char title_string[2048];

            //sprintf(title_string, "Time (days): %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
            sprintf(title_string, "Time (days): %f (%.2f d), k: %i, dt: %.3e, Vis: %s, MaxVal: %.6e, MinVal: %.6e ",
                    shackTimestepControl->currentSimulationTime,
                    shackTimestepControl->currentSimulationTime/(60.0*60.0*24.0),
                    shackTimestepControl->currentTimestepNr,
                    shackTimestepControl->currentTimestepSize,
                    description,
                    //shackPDESWECart2DDiagnostics->total_mass,
                    //shackPDESWECart2DDiagnostics->total_energy,
                    //shackPDESWECart2DDiagnostics->total_potential_enstrophy,
                    vis_cart2d_data.grid_reduce_max(),
                    vis_cart2d_data.grid_reduce_min()
                );

            return title_string;
        }



        void vis_pause()
        {
            shackTimestepControl->runSimulationTimesteps = !shackTimestepControl->runSimulationTimesteps;
        }



        void vis_keypress(int i_key)
        {
            switch(i_key)
            {
            case 'v':
                vis_data_id++;
                break;

            case 'V':
                vis_data_id--;
                break;
            }
        }
#endif


        bool instability_detected() {
            return !(dataConfigOps.prog.h_pert.toGrid().grid_reduce_boolean_all_finite() &&
                     dataConfigOps.prog.u.toGrid().grid_reduce_boolean_all_finite() &&
                     dataConfigOps.prog.v.toGrid().grid_reduce_boolean_all_finite()
            );
        }

    };

}

#endif
