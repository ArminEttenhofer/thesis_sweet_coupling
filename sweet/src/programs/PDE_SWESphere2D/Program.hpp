/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#define WITH_TIME_OLD    1


#ifndef PROGRAMS_PDE_SWESPHERE2D_PROGRAM_HPP
#define PROGRAMS_PDE_SWESPHERE2D_PROGRAM_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#include "Benchmarks/Shack.hpp"
#include "Shack.hpp"


#if SWEET_GUI
#include <sweet/GUI/VisSweet.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#endif

#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>

#include <sweet/Tools/StopwatchBox.hpp>

// Benchmarks
#include "BenchmarksCombined.hpp"

#if WITH_TIME_OLD
// Old time steppers
#include "TimeOld/PDESWESphere2D_TimeSteppers.hpp"

#endif

// New Time Tree
#include "TimeTree/TimeTree.hpp"

#include "DataContainer/Simulation.hpp"
#include "DataContainer/Topography.hpp"

// Diagnostics
#include "Diagnostics.hpp"

// Normal mode analysis
#include "NormalModeAnalysis.hpp"

// File writer
#include "FileOutput.hpp"


namespace PDE_SWESphere2D {

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

            sweet::Data::Sphere2D::Config sphere2DDataConfig;
            sweet::Data::Sphere2D::Operators ops;
            sweet::Data::Sphere2DComplex::Operators opsComplex;

            DataContainer::Simulation prog;
            DataContainer::Simulation progTmp;

            DataContainer::Topography topography;

            sweet::Data::Sphere2D::DataSpectral t0_prog_phi_pert;
            sweet::Data::Sphere2D::DataSpectral t0_prog_div;
            sweet::Data::Sphere2D::DataSpectral t0_prog_vrt;


#if SWEET_GUI
            sweet::Data::Cart2D::Config cart2DDataConfig;

            // Data to visualize is stored to this variable
            sweet::Data::Cart2D::DataGrid vis_cart2d_data;

            // Which primitive to use for rendering
            int vis_render_type_of_primitive_id = 1;

            // Which primitive to use for rendering
            int vis_data_id = 0;
#endif


            bool setup(
                    sweet::Data::Sphere2D::Shack *i_shackSphere2DDataOps,
                    bool i_setup_spectral_transforms = true        // for reset()
            ) {
                /*
                 * Setup Sphere2D Data Config & Operators
                 */
                if (i_setup_spectral_transforms) {
                    sphere2DDataConfig.setupAuto(i_shackSphere2DDataOps);
                    ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DDataConfig);
                }

                ops.setup(&sphere2DDataConfig, i_shackSphere2DDataOps);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

                opsComplex.setup(&sphere2DDataConfig, i_shackSphere2DDataOps);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(opsComplex);

                prog.setup(&sphere2DDataConfig);
                progTmp.setup(&sphere2DDataConfig);

                topography.setup(&sphere2DDataConfig);

#if SWEET_GUI
                sweet::Data::Cart2D::Shack shackCart2DDataOps;
#if 0
                // WARNING: We need to use sphere2DDataConfig, since i_shackSphere2DDataOps is not initialized with a reset()
                shackCart2DDataOps.space_res_physical[0] = i_shackSphere2DDataOps->space_res_physical[0];
                shackCart2DDataOps.space_res_physical[1] = i_shackSphere2DDataOps->space_res_physical[1];
#else
                shackCart2DDataOps.space_res_physical[0] = sphere2DDataConfig.grid_num_lon;
                shackCart2DDataOps.space_res_physical[1] = sphere2DDataConfig.grid_num_lat;

#endif
                shackCart2DDataOps.reuse_spectral_transformation_plans = i_shackSphere2DDataOps->reuse_spectral_transformation_plans;

                cart2DDataConfig.setupAuto(shackCart2DDataOps);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);
#endif

                return true;
            }

            void clear(bool i_clear_spectral_transforms = true) {
                progTmp.clear();
                prog.clear();

                topography.clear();

                t0_prog_phi_pert.clear();
                t0_prog_div.clear();
                t0_prog_vrt.clear();

                ops.clear();

                if (i_clear_spectral_transforms)
                    sphere2DDataConfig.clear();
            }
        };

        // Simulation data
        DataConfigOps dataConfigOps;

        // time integrators
#if WITH_TIME_OLD
        PDESWESphere2D_TimeSteppers timeSteppers;
#endif
        TimeTree::TimeTree timeSteppersNewTS;
        bool useNewTimeSteppers;

        // Handler to all benchmarks
        Benchmarks::BenchmarksCombined sphere2DBenchmarks;

        FileOutput fileOutput;

        /*
         * Shack directory and shacks to work with
         */
        sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
        sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
        sweet::IO::Shack *shackIOData;
        sweet::TimeTree::Shack *shackTimestepControl;
        sweet::Parallelization::Shack *shackParallelization;
        ShackTimeDiscretization *shackTimeDisc;
        Benchmarks::Shack *shackBenchmarks;
        PDE_SWESphere2D::Shack *shackPDESWESphere2D;

        double simulationTimeLastOutputSimtime;

        Diagnostics diagnostics;


    public:
        Program(
                int i_argc,
                char *const *const i_argv
        ) :
                useNewTimeSteppers(true),
                shackProgArgDict(i_argc, i_argv),
                shackSphere2DDataOps(nullptr),
                shackIOData(nullptr),
                shackTimestepControl(nullptr),
                shackParallelization(nullptr),
                shackTimeDisc(nullptr),
                shackBenchmarks(nullptr),
                shackPDESWESphere2D(nullptr) {
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
        }


        bool setup_1_shackRegistration() {
            /*
             * Setup argument parsing
             */
            shackProgArgDict.setup();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            /*
             * SHACK: Register classes which we require
             */
            shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
            shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
            shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
            shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
            shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackTimeDiscretization>();
            shackBenchmarks = shackProgArgDict.getAutoRegistration<Benchmarks::Shack>();
            shackPDESWESphere2D = shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Shack>();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            /*
             * SHACK: Register other things before parsing program arguments
             */

            /*
             * Setup benchmarks
             */
            sphere2DBenchmarks.setup_1_registerAllBenchmark();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);

            sphere2DBenchmarks.setup_2_shackRegistration(&shackProgArgDict);
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);

            shackProgArgDict.processProgramArgumentsForShack(shackTimeDisc);


            if (
                    shackTimeDisc->timestepping_method.find('(') != std::string::npos ||
                    shackTimeDisc->timestepping_method == "help" ||
                    shackTimeDisc->timestepping_method == "helpall"
                    )
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
                timeSteppers.setup_1_registerAllTimesteppers();
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

                timeSteppers.setup_2_shackRegistration(&shackProgArgDict);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);
#endif
            }

            /*
             * Process HELP arguments
             */
            shackProgArgDict.processHelpArguments();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            if (!useNewTimeSteppers) {
                /*
                 * Close shack registration & getting shacks
                 */
                shackProgArgDict.closeRegistration();

                shackProgArgDict.closeGet();
            }

            return true;
        }

        void clear_1_shackRegistration() {
            shackSphere2DDataOps = nullptr;
            shackTimestepControl = nullptr;
            shackParallelization = nullptr;
            shackIOData = nullptr;
            shackTimeDisc = nullptr;

            sphere2DBenchmarks.clear();
#if WITH_TIME_OLD
            timeSteppers.clear();
#endif
            shackProgArgDict.clear();
        }


        bool setup_2_processArguments() {
            /*
             * SHACK: Process arguments
             */
            shackProgArgDict.processProgramArguments();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);


            /*
             * BENCHMARK: Detect particular benchmark to use
             */
            sphere2DBenchmarks.setup_3_benchmarkDetection();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);

            /*
             * Setup benchmark itself
             */
            sphere2DBenchmarks.setup_4_benchmarkSetup_1_withoutOps();


            return true;
        }

        void clear_2_processArguments() {
            shackProgArgDict.clear();
        }


        bool setup_3_dataAndOps(bool i_setup_spectral_transforms = true) {
            /*
             * Setup the data fields
             */
            dataConfigOps.setup(shackSphere2DDataOps, i_setup_spectral_transforms);
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

            /*
             * Setup benchmark itself
             */
            sphere2DBenchmarks.setup_5_benchmarkSetup_2_withOps(&dataConfigOps.ops);

            /*
             * Now we're ready to setup the time steppers
             */
            if (useNewTimeSteppers) {
                bool retval = timeSteppersNewTS.setup_3_timestepper(
                        shackTimeDisc->timestepping_method,
                        &shackProgArgDict,
                        &dataConfigOps.ops,
                        &dataConfigOps.opsComplex,
                        dataConfigOps.prog
                );

                if (!retval) {
                    int helpVerbosity = 0;
                    if (shackTimeDisc->timestepping_method == "helpall") {
                        helpVerbosity = 1;
                        error.clear();

                        timeSteppersNewTS.outputHelp(std::cout, "", helpVerbosity);
                    } else {
                        ERROR_FORWARD(timeSteppersNewTS);
                    }

                    std::cout << "Finishing now..." << std::endl;
                    return false;
                }

                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);

                if (shackIOData->verbosity > 2) {
                    if (shackParallelization->isMPIRoot)
                        timeSteppersNewTS.print();
                }

                /*
                 * Everything is setup now, hence we also set the time step size
                 */
                timeSteppersNewTS.timeIntegrator->setTimeStepSize(shackTimestepControl->currentTimestepSize);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);


            } else {
#if WITH_TIME_OLD
                timeSteppers.setup_3_timestepper(
                        shackTimeDisc->timestepping_method,
                        &shackProgArgDict,
                        &dataConfigOps.ops
                );
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);
#endif
            }

            if (useNewTimeSteppers) {
                shackProgArgDict.closeRegistration();
                shackProgArgDict.closeGet();
            }


            /*
             * Load topography
             */
            sphere2DBenchmarks.benchmark->setup_topography();

            /*
             * Load initial state of benchmark
             */
            sphere2DBenchmarks.benchmark->getInitialState(
                    dataConfigOps.prog.phi_pert,
                    dataConfigOps.prog.vrt,
                    dataConfigOps.prog.div
            );

            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);

            /*
             * Do some validation of program arguments
             */
            shackTimestepControl->validateTimestepSize();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);


            /*
             * Backup data at t=0
             */
            dataConfigOps.t0_prog_phi_pert = dataConfigOps.prog.phi_pert;
            dataConfigOps.t0_prog_vrt = dataConfigOps.prog.vrt;
            dataConfigOps.t0_prog_div = dataConfigOps.prog.div;

            /*
             * Finish registration & getting class interfaces so that nobody can do some
             * strange things with this anymore
             */
            shackProgArgDict.closeRegistration();
            shackProgArgDict.closeGet();

            /*
             * Now we should check that all program arguments have really been parsed
             */
            shackProgArgDict.checkAllArgumentsProcessed();
            ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

            if (shackPDESWESphere2D->compute_diagnostics)
                diagnostics.setup(&(dataConfigOps.ops), shackPDESWESphere2D, 0);

            fileOutput.setup(shackIOData, shackTimestepControl, shackPDESWESphere2D);

            simulationTimeLastOutputSimtime = -1;

            return true;
        }

        void clear_3_data(bool i_clear_spectral_transforms = true) {
            fileOutput.clear();

#if SWEET_GUI
            dataConfigOps.vis_cart2d_data.clear();
#endif

            if (useNewTimeSteppers) {
                timeSteppersNewTS.clear();
            } else {
#if WITH_TIME_OLD
                timeSteppers.clear();
#endif
            }

            dataConfigOps.clear(i_clear_spectral_transforms);
        }

        bool setup(
                bool i_setup_spectral_transforms = true
        ) {
            sweet::Tools::StopwatchBox::getInstance().main_setup.start();

            if (!setup_1_shackRegistration())
                return false;

            if (!setup_2_processArguments())
                return false;

            if (!setup_3_dataAndOps(i_setup_spectral_transforms)) {
                return false;
            }

            if (shackParallelization->isMPIRoot) {
                std::cout << "Printing shack information:" << std::endl;

                shackProgArgDict.printShackData();

                std::cout << "SETUP FINISHED" << std::endl;
            }

            sweet::Tools::StopwatchBox::getInstance().main_setup.stop();

            /*
             * Output data for the first time step as well if output of datafiels is requested
             */
            if (shackIOData->outputEachSimTime >= 0)
                _timestepDoOutput();

            //SWEET_ASSERT(dataConfigOps.cart2DDataConfig.)
            return true;
        }

        void clear(bool i_clear_spectral_transforms = true) {
            clear_3_data(i_clear_spectral_transforms);
            clear_2_processArguments();
            clear_1_shackRegistration();
        }

        bool reset() {
            // keep pausing simulation
            bool run_simulation_timesteps = shackTimestepControl->runSimulationTimesteps;

            clear(false);

            if (!setup(false)) {
                error.print();
                return false;
            }

            shackTimestepControl->runSimulationTimesteps = run_simulation_timesteps;

            return !error.exists();
        }

        void printSimulationErrors() {
            sweet::Data::Sphere2D::DataSpectral diff = dataConfigOps.t0_prog_phi_pert - dataConfigOps.prog.phi_pert;

            double lmax_error = diff.toGrid().grid_reduce_max_abs();
            double rms_error = diff.toGrid().grid_reduce_rms();

            if (shackParallelization->isMPIRoot) {
                std::cout << "Error compared to initial condition" << std::endl;
                std::cout << "Lmax error: " << lmax_error << std::endl;
                std::cout << "RMS error: " << rms_error << std::endl;
            }
        }

        ~Program() {
            clear();
        }

        double precice_max_dt = -1;

        bool runTimestep() {
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
                timeSteppers.timestepper->runTimestep(
                        dataConfigOps.prog.phi_pert, dataConfigOps.prog.vrt, dataConfigOps.prog.div,
                        shackTimestepControl->currentTimestepSize,
                        shackTimestepControl->currentSimulationTime
                );
#endif
            }

            shackTimestepControl->timestepHelperEnd();

            if (shackIOData->verbosity > 2)
                if (shackParallelization->isMPIRoot) {
                    double output_time = shackTimestepControl->currentSimulationTime *
                                         shackTimestepControl->outputSimulationTimeMultiplier;
                    std::cout << shackTimestepControl->currentTimestepNr << ": " << output_time << std::endl;
                }

            if (shackPDESWESphere2D->compute_diagnostics)
                update_diagnostics();

            return true;
        }


        bool should_quit() {
            return shackTimestepControl->isFinalTimestepReached();
        }


        void update_diagnostics() {
            if (diagnostics.last_update_timestep_nr == shackTimestepControl->currentTimestepNr)
                return;

            SWEET_ASSERT(shackPDESWESphere2D->compute_diagnostics);

            diagnostics.update_phi_vrt_div_2_mass_energy_enstrophy(
                    &dataConfigOps.ops,
                    dataConfigOps.prog.phi_pert,
                    dataConfigOps.prog.vrt,
                    dataConfigOps.prog.div,

                    shackSphere2DDataOps->sphere2d_radius,
                    shackPDESWESphere2D->gravitation
            );
        }


        void _timestepDoOutput() {
            if (shackPDESWESphere2D->compute_diagnostics) {
                update_diagnostics();

                if (shackParallelization->isMPIRoot) {
                    // Print header
                    if (shackTimestepControl->currentTimestepNr == 0)
                        diagnostics.printTabularHeader();

                    diagnostics.printTabularRow(shackTimestepControl->currentSimulationTime);
                }
            }

            if (shackPDESWESphere2D->compute_errors) {
                /*
                 * Check for stationary solutions
                 */
                if (
                        shackBenchmarks->benchmark_name != "williamson2" &&
                        shackBenchmarks->benchmark_name != "williamson2_linear" &&
                        shackBenchmarks->benchmark_name != "galewsky_nobump" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_linear" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_topography" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_1" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_2" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_4" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_8" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_16" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_32" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_64" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_128" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_256" &&
                        shackBenchmarks->benchmark_name != "geostrophic_balance_512"
                        ) {

                    if (shackParallelization->isMPIRoot) {
                        std::cout << "Benchmark name: " << shackBenchmarks->benchmark_name << std::endl;
                    }
                    SWEETErrorFatal("Analytical solution not available for this benchmark");
                }

                sweet::Data::Sphere2D::DataSpectral anal_solution_phi_pert(dataConfigOps.sphere2DDataConfig);
                sweet::Data::Sphere2D::DataSpectral anal_solution_vrt(dataConfigOps.sphere2DDataConfig);
                sweet::Data::Sphere2D::DataSpectral anal_solution_div(dataConfigOps.sphere2DDataConfig);

                sphere2DBenchmarks.benchmark->getInitialState(anal_solution_phi_pert, anal_solution_vrt,
                                                              anal_solution_div);

                /*
                 * Compute difference
                 */
                sweet::Data::Sphere2D::DataSpectral diff_phi = dataConfigOps.prog.phi_pert - anal_solution_phi_pert;
                sweet::Data::Sphere2D::DataSpectral diff_vrt = dataConfigOps.prog.vrt - anal_solution_vrt;
                sweet::Data::Sphere2D::DataSpectral diff_div = dataConfigOps.prog.div - anal_solution_div;

                double error_phi = diff_phi.toGrid().grid_reduce_max_abs();
                double error_vrt = diff_vrt.toGrid().grid_reduce_max_abs();
                double error_div = diff_div.toGrid().grid_reduce_max_abs();

                if (shackParallelization->isMPIRoot) {
                    int nTimeSteps = shackTimestepControl->currentTimestepNr;
                    std::cout << "[MULE] errors." << std::setw(8) << std::setfill('0') << nTimeSteps << ": ";

                    std::cout << "simtime=" << shackTimestepControl->currentSimulationTime;
                    std::cout << "\terror_linf_phi=" << error_phi;
                    std::cout << "\terror_linf_vrt=" << error_vrt;
                    std::cout << "\terror_linf_div=" << error_div;
                    std::cout << std::endl;
                }
            }

            if (shackParallelization->isMPIRoot) {
                fileOutput.write_file_output(
                        dataConfigOps.ops,
                        dataConfigOps.prog.phi_pert,
                        dataConfigOps.prog.div,
                        dataConfigOps.prog.vrt
                );

                if (shackIOData->verbosity > 0) {
                    double progPhiMin = dataConfigOps.prog.phi_pert.toGrid().grid_reduce_min();
                    double progPhiMax = dataConfigOps.prog.phi_pert.toGrid().grid_reduce_max();

                    if (shackParallelization->isMPIRoot) {
                        std::cout << "prog_phi min/max:\t" << progPhiMin << ", " << progPhiMax << std::endl;
                    }
                }
            }

            if (shackIOData->outputEachSimTime > 0)
                while (shackIOData->_outputNextSimTime <= shackTimestepControl->currentSimulationTime)
                    shackIOData->_outputNextSimTime += shackIOData->outputEachSimTime;
        }


    public:
        bool timestepHandleOutput() {
            if (shackIOData->outputEachSimTime < 0)
                return false;

            if (shackTimestepControl->currentSimulationTime == simulationTimeLastOutputSimtime)
                return false;

            simulationTimeLastOutputSimtime = shackTimestepControl->currentSimulationTime;

            if (shackTimestepControl->currentSimulationTime <
                shackTimestepControl->maxSimulationTime - shackIOData->outputEachSimTime * 1e-10) {
                if (shackIOData->_outputNextSimTime > shackTimestepControl->currentSimulationTime)
                    return false;
            }

            if (shackParallelization->isMPIRoot)
                if (shackIOData->verbosity > 0)
                    std::cout << std::endl;

            _timestepDoOutput();

            return true;
        }


        bool detect_instability() {
            if (dataConfigOps.prog.phi_pert.spectral_is_first_nan_or_inf()) {
                if (shackParallelization->isMPIRoot) {
                    std::cout << "Infinity value detected" << std::endl;
                    std::cerr << "Infinity value detected" << std::endl;
                }
                return true;
            }

            return false;
        }

        void normalmode_analysis() {
            if (useNewTimeSteppers) {
                SWEETErrorFatal("Not supported, yet");
            } else {
                NormalModeAnalysis::normal_mode_analysis(
                        dataConfigOps.prog.phi_pert,
                        dataConfigOps.prog.vrt,
                        dataConfigOps.prog.div,

                        shackIOData,
                        shackTimestepControl,

                        shackPDESWESphere2D->normal_mode_analysis_generation,

                        this,
                        &Program::runTimestep
                );
            }
        }

        void output_timings() {
            if (shackParallelization->isMPIRoot) {
                std::cout << std::endl;
                sweet::Tools::StopwatchBox::getInstance().output();

                std::cout << "***************************************************" << std::endl;
                std::cout << "* Other timing information (direct)" << std::endl;
                std::cout << "***************************************************" << std::endl;
                std::cout << "[MULE] shackTimestepControl.currentTimestepNr: "
                          << shackTimestepControl->currentTimestepNr << std::endl;
                std::cout << "[MULE] shackTimestepControl.currentTimestepSize: "
                          << shackTimestepControl->currentTimestepSize << std::endl;
                std::cout << std::endl;
                std::cout << "***************************************************" << std::endl;
                std::cout << "* Other timing information (derived)" << std::endl;
                std::cout << "***************************************************" << std::endl;
                std::cout << "[MULE] simulation_benchmark_timings.time_per_time_step (secs/ts): "
                          << sweet::Tools::StopwatchBox::getInstance().main_timestepping() /
                             (double) shackTimestepControl->currentTimestepNr << std::endl;
            }
        }


#if SWEET_GUI
        /*!
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


        int max_vis_types = 9;


        void vis_getDataArray(
                const sweet::Data::Cart2D::DataGrid **o_dataArray,
                double *o_aspect_ratio,
                int *o_vis_render_type_of_primitive_id,
                void **o_bogus_data,
                double *o_viz_min,
                double *o_viz_max,
                bool *viz_reset
        )
        {
            // request rendering of sphere2D
            *o_vis_render_type_of_primitive_id = dataConfigOps.vis_render_type_of_primitive_id;
            *o_bogus_data = &dataConfigOps.sphere2DDataConfig;

            int id = dataConfigOps.vis_data_id % max_vis_types;

            switch (id)
            {
                default:

                case 0:
                    dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(sweet::Data::Sphere2D::DataSpectral(dataConfigOps.prog.phi_pert), dataConfigOps.cart2DDataConfig);
                    break;

                case 1:
                    dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(sweet::Data::Sphere2D::DataSpectral(dataConfigOps.prog.vrt), dataConfigOps.cart2DDataConfig);
                    break;

                case 2:
                    dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(sweet::Data::Sphere2D::DataSpectral(dataConfigOps.prog.div), dataConfigOps.cart2DDataConfig);
                    break;

                case 3:
                    dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(shackPDESWESphere2D->h0 + sweet::Data::Sphere2D::DataSpectral(dataConfigOps.prog.phi_pert)/shackPDESWESphere2D->gravitation, dataConfigOps.cart2DDataConfig);
                    break;

                case 4:
                {
                    sweet::Data::Sphere2D::DataGrid u(dataConfigOps.sphere2DDataConfig);
                    sweet::Data::Sphere2D::DataGrid v(dataConfigOps.sphere2DDataConfig);

                    // Don't use Robert, since we're not interested in the Robert formulation here
                    dataConfigOps.ops.vrtdiv_2_uv(dataConfigOps.prog.vrt, dataConfigOps.prog.div, u, v);
                    dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(u, dataConfigOps.cart2DDataConfig);
                    break;
                }

                case 5:
                {
                    sweet::Data::Sphere2D::DataGrid u(dataConfigOps.prog.vrt.sphere2DDataConfig);
                    sweet::Data::Sphere2D::DataGrid v(dataConfigOps.prog.vrt.sphere2DDataConfig);

                    // Don't use Robert, since we're not interested in the Robert formulation here
                    dataConfigOps.ops.vrtdiv_2_uv(dataConfigOps.prog.vrt, dataConfigOps.prog.div, u, v);
                    dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(v, dataConfigOps.cart2DDataConfig);
                    break;
                }

                case 6:
                case 7:
                case 8:
                {
                    sweet::Data::Sphere2D::DataSpectral anal_solution_phi_pert(dataConfigOps.sphere2DDataConfig);
                    sweet::Data::Sphere2D::DataSpectral anal_solution_vrt(dataConfigOps.sphere2DDataConfig);
                    sweet::Data::Sphere2D::DataSpectral anal_solution_div(dataConfigOps.sphere2DDataConfig);

                    sphere2DBenchmarks.benchmark->getInitialState(
                            anal_solution_phi_pert,
                            anal_solution_vrt,
                            anal_solution_div
                        );

                    switch (id)
                    {
                    case 6:
                        dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(dataConfigOps.prog.phi_pert - anal_solution_phi_pert, dataConfigOps.cart2DDataConfig);
                        break;

                    case 7:
                        dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(dataConfigOps.prog.vrt - anal_solution_vrt, dataConfigOps.cart2DDataConfig);
                        break;

                    case 8:
                        dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(dataConfigOps.prog.div - anal_solution_div, dataConfigOps.cart2DDataConfig);
                        break;
                    }
                }
            }

            double viz_min = dataConfigOps.vis_cart2d_data.grid_reduce_min();
            double viz_max = dataConfigOps.vis_cart2d_data.grid_reduce_max();

            viz_max = std::max(std::abs(viz_max), std::abs(viz_min));
            viz_min = -viz_max;

            *o_viz_min = viz_min;
            *o_viz_max = viz_max;


            *o_dataArray = &dataConfigOps.vis_cart2d_data;
            *o_aspect_ratio = 0.5;
        }



        /*!
         * return status string for window title
         */
        const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
        {
            std::ostringstream ss;
            std::string sep = ",  ";

            o_replace_commas_with_newline = true;

            if (shackParallelization->useMPI) {
                ss << "Rank=" << shackParallelization->mpiRank;
                ss << ",";
            }

            const char *fields_array[] = {
                    "phi_pert",
                    "vrt",
                    "div",
                    "h",
                    "u",
                    "v",
                    "phi diff t0",
                    "u diff t0",
                    "v diff t0",
            };

            int days = ((int)(shackTimestepControl->currentSimulationTime/(60.0*60.0*24.0)));
            int hours = ((int)(shackTimestepControl->currentSimulationTime/(60.0*60.0)))%24;
            int minutes = ((int)(shackTimestepControl->currentSimulationTime/(60.0)))%60;
            int seconds = ((int)shackTimestepControl->currentSimulationTime) % 60;

            ss << "time="
                    << days
                    << "d-"
                    << (hours < 10 ? "0" : "") << hours
                    << "h:"
                    << (minutes < 10 ? "0" : "") << minutes
                    << "m:"
                    << (seconds < 10 ? "0" : "") << seconds
                    << "s"
                    << sep;


            int id = (dataConfigOps.vis_data_id % max_vis_types + max_vis_types) % max_vis_types;
            ss << "id=" << dataConfigOps.vis_data_id << sep;
            ss << "field=" << fields_array[id] << sep;
            ss << "max=" << dataConfigOps.vis_cart2d_data.grid_reduce_max() << sep;
            ss << "min=" << dataConfigOps.vis_cart2d_data.grid_reduce_min() << sep;

            ss << "time.step.nr="	<< shackTimestepControl->currentTimestepNr << sep;
            ss << "time.step.size="	<< shackTimestepControl->currentTimestepSize << sep;

            if (shackPDESWESphere2D->compute_diagnostics)
            {
                update_diagnostics();

                ss << "diag.total_mass="	<< diagnostics.total_mass << sep;
                ss << "diag.total_energy="	<< diagnostics.total_energy << sep;
                ss << "diag.total_potential_enstrophy="	<< diagnostics.total_potential_enstrophy;
                ss << ",  ";
            }
            ss << "FIN";

            return ss.str();
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
                dataConfigOps.vis_data_id++;
                break;

            case 'V':
                dataConfigOps.vis_data_id--;
                break;

            case 'b':
                dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id + 1) % 2;
                break;

            case 'B':
                dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id - 1) % 2;
                break;
            }
        }
#endif
    };

}

#endif
