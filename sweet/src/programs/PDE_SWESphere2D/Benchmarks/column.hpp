#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_COLUMN_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_COLUMN_HPP

#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/LibMath/GaussQuadrature.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
    namespace Benchmarks {

        class column : public BaseInterface {

        private:
            // Uses Haversine formula
            double great_circle_distance(double lat1, double lon1, double lat2, double lon2, double r) {
                auto dlat = lat1 - lat2;
                auto dlon = lon1 - lon2;

                double res = 1 - std::cos(dlat) + std::cos(lat1) * std::cos(lat2) * (1 - std::cos(dlon));
                res = res / 2;
                res = std::sqrt(res);
                res = 2 * r * std::asin(res);

                return res;
            }

        public:
            column() {
            }

            bool shackRegistration(
                    sweet::Shacks::Dictionary *io_shackDict
            ) {
                BaseInterface::shackRegistration(io_shackDict);
                return true;
            }

            std::string benchmark_name;

            bool implements_benchmark(
                    const std::string &i_benchmark_name
            ) {
                benchmark_name = i_benchmark_name;

                return i_benchmark_name == "column";
            }

            void setup_1_shackData() {
                if (shackParallelization->isMPIRoot) {
                    std::cout << "!!! WARNING !!!" << std::endl;
                    std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
                    std::cout << "!!! WARNING !!!" << std::endl;
                }

                // Setup world parameters
                shackPDESWESphere2D->sphere2d_rotating_coriolis_omega = 0;
                shackPDESWESphere2D->gravitation = 9.80616;
                shackSphere2DDataOps->sphere2d_radius = 1e5 / (2 * M_PI); // Match size to Cart2D (100km circumference)
//                shackPDESWESphere2D->viscosity = 1e-6;
                shackPDESWESphere2D->h0 = 1000;
            }

            void setup_2_withOps(
                    sweet::Data::Sphere2D::Operators *io_ops
            ) {
                ops = io_ops;
            }


            void clear() {}

            std::string getHelp() {
                std::ostringstream stream;
                stream << "  'column': 100m column of water" << std::endl;
                return stream.str();
            }

            void getInitialState(
                    sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
                    sweet::Data::Sphere2D::DataSpectral &o_vrt,
                    sweet::Data::Sphere2D::DataSpectral &o_div
            ) {
                const sweet::Data::Sphere2D::Config *sphere2DDataConfig = o_phi_pert.sphere2DDataConfig;

                constexpr double col_lat = 0;
                constexpr double col_lon = M_PI;
                constexpr double col_r = 7500;
                constexpr double col_h = 200;


                /*
                 * Setup V=0
                 */
                o_vrt.spectral_setZero();
                o_div.spectral_setZero();

                auto phi_phys = o_phi_pert.toGrid();
                phi_phys.grid_setZero();
                phi_phys.grid_update_lambda(
                        [&](double lon, double lat, double &o_data) {
                            double dist = great_circle_distance(col_lat, col_lon, lat, lon,
                                                                shackSphere2DDataOps->sphere2d_radius);

                            double actual_height = 0.5 * (std::cos((dist/col_r) * M_PI) + 1) * col_h;

                            if (dist < col_r) {
                                o_data = actual_height * shackPDESWESphere2D->gravitation;
                            } else {
                                o_data = 0;
                            }
                        }
                );

                o_phi_pert.loadSphere2DDataGrid(phi_phys);
            }

            void setup_topography() override {
                // initialize the topography
                shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
            }
        };
    }
}

#endif
