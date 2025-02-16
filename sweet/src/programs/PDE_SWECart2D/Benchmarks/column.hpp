//
// Created by armin on 11/13/24.
//

#ifndef SWEET_PRECICE_COUPLING_COLUMN_HPP
#define SWEET_PRECICE_COUPLING_COLUMN_HPP

#include <programs/PDE_SWECart2D/Benchmarks/BaseInterface.hpp>
#include <stdlib.h>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <cmath>
#include <sweet/LibMath/GaussQuadrature.hpp>
#include <sweet/Shacks/Dictionary.hpp>


namespace PDE_SWECart2D {
    namespace Benchmarks {

        class column : public PDE_SWECart2D::Benchmarks::BaseInterface {
            double f;
            double g;
            double sx;
            double sy;
            const double height = 200;
            const double pos_x = 0.5;
            const double pos_y = 0.5;
            const double size = 0.075;

        public:
            void setup_depth(
                    sweet::Data::Cart2D::DataSpectral &o_depth,
                    bool i_with_bump = true
            ) {

                sweet::Data::Cart2D::DataGrid depth_phys(o_depth.cart2DDataConfig);

                for (int j = 0; j < shackCart2DDataOps->space_res_physical[1]; j++) {
                    for (int i = 0; i < shackCart2DDataOps->space_res_physical[0]; i++) {

                        double x =
                                (i + 0.5) / shackCart2DDataOps->space_res_physical[0]; //*shackDict.sim.domain_size[0];
                        double y =
                                (j + 0.5) / shackCart2DDataOps->space_res_physical[1]; //*shackDict.sim.domain_size[1];

                        double x_dist = std::abs(x - pos_x);
                        double y_dist = std::abs(y - pos_y);

                        // Create a circular column of constant height
                        double dist = std::sqrt(x_dist * x_dist + y_dist * y_dist);
                        double actual_height = 0.5 * (std::cos((dist/size) * M_PI) + 1) * height;
                        if (dist <= size) {
                            depth_phys.grid_setValue(j, i, depth_phys.grid_get(j, i) + actual_height);
                        }
                    }
                }

                o_depth.loadCart2DDataGrid(depth_phys);
            }

            void setup_velocity(
                    sweet::Data::Cart2D::DataSpectral &o_u,
                    sweet::Data::Cart2D::DataSpectral &o_v
            ) {

                o_v.spectral_setZero();
                o_u.spectral_setZero();
            }

        public:
            bool setupBenchmark(
                    sweet::Data::Cart2D::DataSpectral &o_h,
                    sweet::Data::Cart2D::DataSpectral &o_u,
                    sweet::Data::Cart2D::DataSpectral &o_v
            ) {
                f = shackPDESWECart2D->cart2d_rotating_f0;
                g = shackPDESWECart2D->gravitation;
                sx = shackCart2DDataOps->cart2d_domain_size[0];
                sy = shackCart2DDataOps->cart2d_domain_size[1];


                std::cout << "Generating column initial conditions.";

                setup_velocity(o_u, o_v);
                setup_depth(o_h);

                std::cout << "   Done! " << std::endl;

                return true;
            }


        };


    }
}

#endif //SWEET_PRECICE_COUPLING_COLUMN_HPP
