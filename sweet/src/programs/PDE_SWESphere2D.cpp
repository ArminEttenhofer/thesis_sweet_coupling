/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeOld
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeTree
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeTree/TimeStepper
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/Benchmarks
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <precice/precice.hpp>
#include <eigen3/Eigen/Core>
#include "PDE_SWECart2D/CouplingUtility.hpp"
#include <iostream>
#include <sys/stat.h>
#include <errno.h>


#if SWEET_MPI
#	include <mpi.h>
#endif

#include "PDE_SWESphere2D/Program.hpp"


#if SWEET_MPI
int mpi_comm_rank;
int mpi_comm_size;
#endif

bool isMPIRoot() {
#if SWEET_MPI
    return mpi_comm_rank == 0;
#else
    return true;
#endif
}


int main_mpi(int i_argc, char *i_argv[]) {
#if SWEET_MPI
#if SWEET_THREADING
    int provided;
    MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE)
        SWEETErrorFatal("MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye...");
#else
    MPI_Init(&i_argc, &i_argv);
#endif

MPI_Comm_rank(MPI_COMM_WORLD, &mpi_comm_rank);
MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
#endif

    PDE_SWESphere2D::Program simulation(i_argc, i_argv);
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

    if (!simulation.setup()) {
        SWEET_ASSERT(simulation.error.exists());
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);
    }

    Coupling_Shack coupling{};
    simulation.shackProgArgDict.processProgramArgumentsForShack(&coupling);

    bool PRECICE = true;
    if (PRECICE) {
        coupling.use_precice = true;
        coupling.couple_alpha = true;
    }


    {
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

#if SWEET_GUI
        if (simulation.shackIOData->guiEnabled)
        {
            sweet::GUI::VisSweet visSweet(simulation);
        }
        else
#endif
        {
            std::vector<double> timestamps{};
            std::vector<double> point_a{};
            int a_x = 0;
            int b_x = 25000;
            int c_x = 45000;
            std::vector<double> point_b{};
            std::vector<double> point_c{};

            simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
            ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

            if (simulation.shackPDESWESphere2D->normal_mode_analysis_generation > 0) {
                simulation.normalmode_analysis();
            } else {
                simulation.timestepHandleOutput();

                constexpr int FOAM_centerX = 25000;
                constexpr int FOAM_centerY = 0;
                constexpr int FOAM_sizeX = 15000;
                constexpr int FOAM_sizeY = 15000;
                constexpr int OFFSET = 4;

                // Initialize coupling
                std::string config_name{"../../precice-config.xml"};
                if (coupling.couple_alpha) config_name = "../../precice-config-alpha.xml";
                precice::Participant precice{"SWEET", config_name, 0, 1};

                // Init preCICE mesh
                std::vector<double> coordinates{};
                std::vector<int> vertexIDs_write{};
                int long_size = simulation.shackSphere2DDataOps->space_res_physical[0];
                int lat_size = simulation.shackSphere2DDataOps->space_res_physical[1];
                double radius = simulation.shackSphere2DDataOps->sphere2d_radius;
                double circumference = 2 * M_PI * radius;

                double dx = circumference / long_size;
                double dy = circumference / 2 / lat_size;

                a_x = std::floor(0.5 * circumference / dx);
                b_x = std::floor((0.5 * circumference + b_x) / dx);
                c_x = std::floor((0.5 * circumference + c_x) / dx);

                int latSouth = (int) std::ceil((FOAM_centerY - FOAM_sizeY / 2. + 0.25 * circumference) / dy) - 1;
                int latNorth = (int) std::ceil((FOAM_centerY + FOAM_sizeY / 2. + 0.25 * circumference) / dy) - 1;
                int longWest = (int) std::ceil((FOAM_centerX - FOAM_sizeX / 2. + 0.5 * circumference) / dx) - 1;
                int longEast = (int) std::ceil((FOAM_centerX + FOAM_sizeX / 2. + 0.5 * circumference) / dx) - 1;

                std::cout << "South: " << latSouth << " North: " << latNorth << " West: " << longWest << " East: "
                          << longEast << "\n";

                // North
                for (int i = longWest; i <= longEast; i++) {
                    vertexIDs_write.emplace_back(10000 * 0 + 100 * i);
                    coordinates.emplace_back(i * dx + dx / 2);
                    coordinates.emplace_back(latNorth * dy + dy / 2 + 0.25 * circumference);
                    coordinates.emplace_back(0);
                }

                // East
                for (int i = latSouth; i <= latNorth; i++) {
                    vertexIDs_write.emplace_back(10000 * 1 + 100 * i);
                    coordinates.emplace_back(longEast * dx + dx / 2);
                    coordinates.emplace_back(i * dy + dy / 2 + 0.25 * circumference);
                    coordinates.emplace_back(0);
                }

                // South
                for (int i = longWest; i <= longEast; i++) {
                    vertexIDs_write.emplace_back(10000 * 2 + 100 * i);
                    coordinates.emplace_back(i * dx + dx / 2);
                    coordinates.emplace_back(latSouth * dy + dy / 2 + 0.25 * circumference);
                    coordinates.emplace_back(0);
                }

                // West
                for (int i = latSouth; i <= latNorth; i++) {
                    vertexIDs_write.emplace_back(10000 * 3 + 100 * i);
                    coordinates.emplace_back(longWest * dx + dx / 2);
                    coordinates.emplace_back(i * dy + dy / 2 + 0.25 * circumference);
                    coordinates.emplace_back(0);
                }

                const double h0 = simulation.shackPDESWESphere2D->h0;

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
                for (int i = longWest + OFFSET; i <= longEast - OFFSET; ++i) {
                    for (int j = latSouth + OFFSET; j <= latNorth - OFFSET; ++j) {
                        vertexIDs_read.emplace_back(i * FOAM_sizeY + j);
                        coordinates.emplace_back(dx * (0.5 + i));
                        coordinates.emplace_back(dy * (0.5 + j) + 0.25 * circumference);
                        coordinates.emplace_back(simulation.shackPDESWESphere2D->h0);
                    }
                }

                if (coupling.use_precice) {
                    precice.setMeshVertices("SWE_Top", coordinates, vertexIDs_read);
                    precice.initialize();
                }

                {
                    auto h_pert_temp = simulation.dataConfigOps.prog.phi_pert.toGrid() /
                                       simulation.shackPDESWESphere2D->gravitation;
                    timestamps.push_back(simulation.shackTimestepControl->currentSimulationTime);
                    point_a.push_back(h_pert_temp.grid_getValue(a_x, lat_size / 2));
                    point_b.push_back(h_pert_temp.grid_getValue(b_x, lat_size / 2));
                    point_c.push_back(h_pert_temp.grid_getValue(c_x, lat_size / 2));
                }

                sweet::Tools::StopwatchBox::getInstance().main_timestepping.start();

                while (!simulation.should_quit() && (!coupling.use_precice || precice.isCouplingOngoing())) {
                    if (coupling.use_precice) simulation.precice_max_dt = precice.getMaxTimeStepSize();
                    simulation.runTimestep();

                    if (simulation.shackPDESWESphere2D->instability_checks) {
                        if (isMPIRoot()) {
                            if (simulation.detect_instability()) {
                                std::cerr << "INSTABILITY DETECTED" << std::endl;
                                exit(1);
                            }
                        }
                    }

                    auto h = simulation.dataConfigOps.prog.phi_pert.toGrid() /
                             simulation.shackPDESWESphere2D->gravitation;

                    if (coupling.use_precice) {
                        auto &vrt = simulation.dataConfigOps.prog.vrt;
                        auto &div = simulation.dataConfigOps.prog.div;

                        sweet::Data::Sphere2D::DataGrid u{simulation.dataConfigOps.sphere2DDataConfig};
                        sweet::Data::Sphere2D::DataGrid v{simulation.dataConfigOps.sphere2DDataConfig};
                        simulation.dataConfigOps.ops.vrtdiv_2_uv(vrt, div, u, v);

                        std::vector<double> velocities{};
                        std::vector<double> heights{};

                        // North
                        for (int i = longWest; i <= longEast; i++) {
                            velocities.emplace_back(u.grid_getValue(i, latNorth));
                            velocities.emplace_back(v.grid_getValue(i, latNorth));
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_getValue(i, latNorth) + h0);
                            }
                        }

                        // East
                        for (int i = latSouth; i <= latNorth; i++) {
                            velocities.emplace_back(u.grid_getValue(longEast, i));
                            velocities.emplace_back(v.grid_getValue(longEast, i));
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_getValue(longEast, i) + h0);
                            }
                        }

                        // South
                        for (int i = longWest; i <= longEast; i++) {
                            velocities.emplace_back(u.grid_getValue(i, latSouth));
                            velocities.emplace_back(v.grid_getValue(i, latSouth));
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_getValue(i, latSouth) + h0);
                            }
                        }

                        // West
                        for (int i = latSouth; i <= latNorth; i++) {
                            velocities.emplace_back(u.grid_getValue(longWest, i));
                            velocities.emplace_back(v.grid_getValue(longWest, i));
                            velocities.emplace_back(0);
                            if (coupling.couple_alpha) {
                                heights.emplace_back(h.grid_getValue(longWest, i) + h0);
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

                        int buf_i = 0;
                        for (int i = longWest + OFFSET; i <= longEast - OFFSET; ++i) {
                            for (int j = latSouth + OFFSET; j <= latNorth - OFFSET; ++j) {
                                u.grid_setValue(i, j, velocities.at(3 * buf_i));
                                v.grid_setValue(i, j, velocities.at(1 + 3 * buf_i));
                                h.grid_setValue(i, j, heights.at(buf_i) - h0);
                                buf_i++;
                                if (i == 50 && j == 120) std::cout << heights.at(i * (latNorth - latSouth + 1) + j) << "=============================\n";
                            }
                        }

                        simulation.dataConfigOps.ops.uv_2_vrtdiv(u, v, vrt, div);
                        simulation.dataConfigOps.prog.phi_pert.loadSphere2DDataGrid(
                                h * simulation.shackPDESWESphere2D->gravitation);

                        precice.advance(simulation.shackTimestepControl->currentTimestepSize);
                    }

                    timestamps.push_back(simulation.shackTimestepControl->currentSimulationTime);
                    point_a.push_back(h.grid_getValue(a_x, lat_size / 2));
                    point_b.push_back(h.grid_getValue(b_x, lat_size / 2));
                    point_c.push_back(h.grid_getValue(c_x, lat_size / 2));

                    simulation.timestepHandleOutput();
                }

                if (isMPIRoot())
                    std::cout << "TIMESTEPPING FINISHED" << std::endl;

                sweet::Tools::StopwatchBox::getInstance().main_timestepping.stop();
            }
            if (isMPIRoot()) {
                const char *directory = "points";
                if (mkdir(directory, 0755) == 0 || errno == EEXIST) {
                    std::string filename = "points/points.csv";
                    std::ofstream file(filename);

                    file << std::fixed << std::setprecision(2);

                    file << "Time,Point_A,Point_B,Point_C\n";
                    for (int i = 0; i < timestamps.size(); i++) {
                        file << timestamps[i] << ", " << point_a[i] << ", " << point_b[i] << ", " << point_c[i] << "\n";
                    }

                    file.close();
                }

                if (!simulation.fileOutput.output_reference_filenames.empty())
                    std::cout << "[MULE] reference_filenames: " << simulation.fileOutput.output_reference_filenames
                              << std::endl;
            }
        }

        // End of run output results
        simulation.output_timings();
    }

    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

    if (isMPIRoot())
        std::cout << "FIN" << std::endl;

#if SWEET_MPI
    MPI_Finalize();
#endif

    return 0;
}


int main(int i_argc, char *i_argv[]) {
    int retval = main_mpi(i_argc, i_argv);

    //MemBlockAlloc::shutdown();

    return retval;
}

