/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/TimeOld/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/TimeTree/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/TimeTree/TimeStepper
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/Benchmarks/
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */


#include <programs/PDE_SWECart2D/Program.hpp>
#include <precice/precice.hpp>
#include <eigen3/Eigen/Core>
#include "PDE_SWECart2D/CouplingUtility.hpp"
#include <iostream>
#include <sys/stat.h>
#include <errno.h>

#if SWEET_MPI
int mpi_rank;
#endif

bool isMPIRoot() {
#if SWEET_MPI
    return mpi_rank == 0;
#else
    return true;
#endif
}


int main_mpi(int i_argc, char *i_argv[]) {
    sweet::Tools::StopwatchBox::getInstance().main.start();

    PDE_SWECart2D::Program simulation(i_argc, i_argv);
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

    simulation.setup();

    Coupling_Shack coupling{};
    simulation.shackProgArgDict.processProgramArgumentsForShack(&coupling);

    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

#if SWEET_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

#if SWEET_GUI
    if (simulation.shackIOData->guiEnabled)
    {
        sweet::GUI::VisSweet visSweet(simulation);
    }
    else
#endif
    std::vector<double> timestamps{};
    std::vector<double> point_a{};
    int a_x = 100;
    int b_x = 150;
    int c_x = 190;
    std::vector<double> point_b{};
    std::vector<double> point_c{};

    {
#if SWEET_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);
        if (simulation.shackPDESWECart2D->normal_mode_analysis_generation > 0) {
            simulation.normal_mode_analysis();
        } else {
            simulation.timestep_do_output();

            constexpr int FOAM_originX = 120;
            constexpr int FOAM_originY = 70;
            constexpr int FOAM_sizeX = 60;
            constexpr int FOAM_sizeY = 60;
            constexpr int FOAM_grid_size_z = 5;
            constexpr int OFFSET = 2;


            const double dz = simulation.shackPDESWECart2D->h0 / FOAM_grid_size_z;

            // Initialize coupling
            std::string config_name{"precice-config.xml"};
            if (coupling.couple_alpha) config_name = "precice-config-alpha.xml";
            precice::Participant precice{"SWEET", config_name, 0, 1};
            const auto grid_size = simulation.shackCart2DDataOps->space_res_physical;

            // Init preCICE mesh
            std::vector<double> coordinates{};
            std::vector<int> vertexIDs_write{};
            const auto size_x = simulation.shackCart2DDataOps->cart2d_domain_size[0];
            const auto size_y = simulation.shackCart2DDataOps->cart2d_domain_size[1];
            const double dx = size_x / grid_size[0];
            const double dy = size_y / grid_size[1];

            // North
            for (int i = 0; i < FOAM_sizeX; i++) {
                for (int z = 0; z < FOAM_grid_size_z; z++) {
                    vertexIDs_write.emplace_back(10000 * 0 + 100 * i + z);
                    coordinates.emplace_back(FOAM_originX * dx + i * dx + dx / 2);
                    coordinates.emplace_back(FOAM_originY * dy + dy / 2 - dy);
                    coordinates.emplace_back(z * dz + dz / 2);
                }
            }

            // East
            for (int i = 0; i < FOAM_sizeY; i++) {
                for (int z = 0; z < FOAM_grid_size_z; z++) {
                    vertexIDs_write.emplace_back(10000 * 1 + 100 * i + z);
                    coordinates.emplace_back(FOAM_originX * dx + FOAM_sizeX * dx + dx / 2);
                    coordinates.emplace_back(FOAM_originY * dy + i * dy + dy / 2);
                    coordinates.emplace_back(z * dz + dz / 2);
                }
            }

            // South
            for (int i = 0; i < FOAM_sizeX; i++) {
                for (int z = 0; z < FOAM_grid_size_z; z++) {
                    vertexIDs_write.emplace_back(10000 * 2 + 100 * i + z);
                    coordinates.emplace_back(FOAM_originX * dx + i * dx + dx / 2);
                    coordinates.emplace_back(FOAM_originY * dy + FOAM_sizeY * dy + dy / 2);
                    coordinates.emplace_back(z * dz + dz / 2);
                }
            }

            // West
            for (int i = 0; i < FOAM_sizeY; i++) {
                for (int z = 0; z < FOAM_grid_size_z; z++) {
                    vertexIDs_write.emplace_back(10000 * 3 + 100 * i + z);
                    coordinates.emplace_back(FOAM_originX * dx + dx / 2 - dx);
                    coordinates.emplace_back(FOAM_originY * dy + i * dy + dy / 2);
                    coordinates.emplace_back(z * dz + dz / 2);
                }
            }

            const double h0 = simulation.shackPDESWECart2D->h0;

            if (coupling.use_precice) {
                precice.setMeshVertices("SWE_Sides", coordinates, vertexIDs_write);

                if (coupling.couple_alpha) {
                    std::vector<double> heights(vertexIDs_write.size());
                    std::fill(heights.begin(), heights.end(), h0);
                    if (precice.requiresInitialData()) {
                        precice.writeData("SWE_Sides", "Alpha_SWE_FOAM", vertexIDs_write, heights);
                    }
                }
            }

            coordinates.clear();
            std::vector<int> vertexIDs_read{};
            for (int i = 0; i < FOAM_sizeX; ++i) {
                for (int j = 0; j < FOAM_sizeY; ++j) {
                    vertexIDs_read.emplace_back(i * FOAM_sizeY + j);
                    coordinates.emplace_back(dx * (0.5 + FOAM_originX + i));
                    coordinates.emplace_back(dy * (0.5 + FOAM_originY + j));
                    coordinates.emplace_back(simulation.shackPDESWECart2D->h0);
                }
            }

            if (coupling.use_precice) {
                precice.setMeshVertices("SWE_Top", coordinates, vertexIDs_read);
                precice.initialize();
            }

            sweet::Tools::StopwatchBox::getInstance().main_timestepping.start();

            {
                auto h_pert_temp = simulation.dataConfigOps.prog.h_pert.toGrid();
                timestamps.push_back(simulation.shackTimestepControl->currentSimulationTime);
                point_a.push_back(h_pert_temp.grid_get(100, a_x));
                point_b.push_back(h_pert_temp.grid_get(100, b_x));
                point_c.push_back(h_pert_temp.grid_get(100, c_x));
            }
            while (!simulation.should_quit() && (!coupling.use_precice || precice.isCouplingOngoing())) {
                if (coupling.use_precice) simulation.precice_max_dt = precice.getMaxTimeStepSize();
                simulation.runTimestep();

                // Instability
                if (simulation.shackPDESWECart2D->instability_checks) {
                    if (simulation.instability_detected())
                        SWEETErrorFatal("INSTABILITY DETECTED");
                }

                auto h = simulation.dataConfigOps.prog.h_pert.toGrid();

                if (coupling.use_precice) {
                    auto u = simulation.dataConfigOps.prog.u.toGrid();
                    auto v = simulation.dataConfigOps.prog.v.toGrid();

                    std::vector<double> velocities{};
                    std::vector<double> heights{};

                    // North
                    for (int i = 0; i < FOAM_sizeX; i++) {
                        for (int z = 0; z < FOAM_grid_size_z; z++) {
                            velocities.emplace_back(u.grid_get(FOAM_originY - 1, i + FOAM_originX));
                            velocities.emplace_back(v.grid_get(FOAM_originY - 1, i + FOAM_originX));
                            //                            velocities.emplace_back(0);
                            //                            velocities.emplace_back(0);
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_get(FOAM_originY - 1, i + FOAM_originX) + h0);
                                //                                heights.emplace_back(h0);
                            }
                        }
                    }

                    // East
                    for (int i = 0; i < FOAM_sizeY; i++) {
                        for (int z = 0; z < FOAM_grid_size_z; z++) {
                            velocities.emplace_back(u.grid_get(FOAM_originY + i, FOAM_originX + FOAM_sizeX));
                            velocities.emplace_back(v.grid_get(FOAM_originY + i, FOAM_originX + FOAM_sizeX));
                            //                            velocities.emplace_back(0);
                            //                            velocities.emplace_back(0);
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_get(FOAM_originY + i, FOAM_originX + FOAM_sizeX) + h0);
                                //                                heights.emplace_back(h0);
                                //Test
                            }
                        }
                    }

                    // South
                    for (int i = 0; i < FOAM_sizeX; i++) {
                        for (int z = 0; z < FOAM_grid_size_z; z++) {
                            velocities.emplace_back(u.grid_get(FOAM_originY + FOAM_sizeY, i + FOAM_originX));
                            velocities.emplace_back(v.grid_get(FOAM_originY + FOAM_sizeY, i + FOAM_originX));
                            //                            velocities.emplace_back(0);
                            //                            velocities.emplace_back(0);
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_get(FOAM_originY + FOAM_sizeY, i + FOAM_originX) + h0);
                                //                                heights.emplace_back(h0);
                            }
                        }
                    }

                    // West
                    for (int i = 0; i < FOAM_sizeY; i++) {
                        for (int z = 0; z < FOAM_grid_size_z; z++) {
                            velocities.emplace_back(u.grid_get(FOAM_originY + i, FOAM_originX - 1));
                            velocities.emplace_back(v.grid_get(FOAM_originY + i, FOAM_originX - 1));
                            //                            velocities.emplace_back(0);
                            //                            velocities.emplace_back(0);
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_get(FOAM_originY + i, FOAM_originX - 1) + h0);
                                //                                heights.emplace_back(h0);
                            }
                        }
                    }

                    precice.writeData("SWE_Sides", "Velocity_SWE_FOAM", vertexIDs_write, velocities);
                    if (coupling.couple_alpha) {
                        precice.writeData("SWE_Sides", "Alpha_SWE_FOAM", vertexIDs_write, heights);
                    }

                    velocities = std::vector<double>(3 * vertexIDs_read.size());
                    heights = std::vector<double>(vertexIDs_read.size());
                    precice.readData("SWE_Top", "Velocity_FOAM_SWE", vertexIDs_read, 0, velocities);
                    precice.readData("SWE_Top", "Alpha_FOAM_SWE", vertexIDs_read, 0, heights);

                    for (int i = OFFSET; i < FOAM_sizeX - OFFSET; ++i) {
                        for (int j = OFFSET; j < FOAM_sizeY - OFFSET; ++j) {
                            u.grid_setValue(j + FOAM_originY, i + FOAM_originX,
                                            velocities.at(3 * (i * FOAM_sizeY + j)));
                            v.grid_setValue(j + FOAM_originY, i + FOAM_originX,
                                            velocities.at(1 + 3 * (i * FOAM_sizeY + j)));
                            h.grid_setValue(j + FOAM_originY, i + FOAM_originX, heights.at(i * FOAM_sizeY + j) - h0);
                        }
                    }

                    simulation.dataConfigOps.prog.u.loadCart2DDataGrid(u);
                    simulation.dataConfigOps.prog.v.loadCart2DDataGrid(v);
                    simulation.dataConfigOps.prog.h_pert.loadCart2DDataGrid(h);

                    precice.advance(simulation.shackTimestepControl->currentTimestepSize);
                }

                timestamps.push_back(simulation.shackTimestepControl->currentSimulationTime);
                point_a.push_back(h.grid_get(100, a_x));
                point_b.push_back(h.grid_get(100, b_x));
                point_c.push_back(h.grid_get(100, c_x));
            }

            if (coupling.use_precice) precice.finalize();
            sweet::Tools::StopwatchBox::getInstance().main_timestepping.stop();
        }
    }

    if (isMPIRoot()) {
        ////if (simulation.shackIOData->outputFileName.size() > 0)
        ////	std::cout << "[MULE] reference_filenames: " << simulation.output_filenames << std::endl;


        const char *directory = "output/points";
        if (mkdir(directory, 0755) == 0 || errno == EEXIST) {
            std::string filename = "output/points/points.csv";
            std::ofstream file(filename);

            file << std::fixed << std::setprecision(2);

            file << "Time,Point_A,Point_B,Point_C\n";
            for (int i = 0; i < timestamps.size(); i++) {
                file << timestamps[i] << ", " << point_a[i] << ", " << point_b[i] << ", " << point_c[i] << "\n";
            }

            file.close();
        }

        if (simulation.fileOutput.output_reference_filenames.size() > 0)
            std::cout << "[MULE] reference_filenames: " << simulation.fileOutput.output_reference_filenames <<
                      std::endl;


        // End of run output results
        std::cout << "***************************************************" << std::endl;
        std::cout << "Number of time steps: " << simulation.shackTimestepControl->currentTimestepNr << std::endl;
        std::cout << "Time per time step: " << sweet::Tools::StopwatchBox::getInstance().main_timestepping() / (double)
                simulation.shackTimestepControl->currentTimestepNr << " sec/ts" << std::endl;
        std::cout << "Last time step size: " << simulation.shackTimestepControl->currentTimestepSize << std::endl;

        simulation.computeErrors();
        simulation.printErrors();

        std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;
    }

    sweet::Tools::StopwatchBox::getInstance().main.stop();

#if SWEET_MPI
    if (mpi_rank == 0)
#endif
    {
        std::cout << std::endl;
        sweet::Tools::StopwatchBox::getInstance().output();
    }

    if (isMPIRoot()) {
        std::cout << "FIN" << std::endl;
    }
    return 0;
}


int main(int i_argc, char *i_argv[]) {
#if SWEET_MPI
#if SWEET_THREADING
    int provided;
    MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE)
        SWEETErrorFatal("MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye...");
#else
    MPI_Init(&i_argc, &i_argv);
#endif
#endif

    int retval = main_mpi(i_argc, i_argv);

#if SWEET_MPI
    MPI_Finalize();
#endif

    return retval;
}
