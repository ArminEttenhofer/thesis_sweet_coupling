//
// Created by armin on 11/14/24.
//

#ifndef SWEET_PRECICE_COUPLING_IMPULSE_HPP
#define SWEET_PRECICE_COUPLING_IMPULSE_HPP

//
// Created by armin on 11/13/24.
//

#include <programs/PDE_SWECart2D/Benchmarks/BaseInterface.hpp>
#include <stdlib.h>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <cmath>
#include <sweet/LibMath/GaussQuadrature.hpp>
#include <sweet/Shacks/Dictionary.hpp>


namespace PDE_SWECart2D {
    namespace Benchmarks {

        class impulse : public PDE_SWECart2D::Benchmarks::BaseInterface {
            double f;
            double g;
            double sx;
            double sy;
            const double width = 0.1;
            const double pos_x = 0.2;
            const double pos_y = 0.5;
            const double speed = 30;
            const double length = 0.4;

        public:
            void setup_depth(
                    sweet::Data::Cart2D::DataSpectral &o_depth
            ) {
                sweet::Data::Cart2D::DataGrid depth_phys(o_depth.cart2DDataConfig);
                o_depth.loadCart2DDataGrid(depth_phys);
            }

            void setup_velocity(
                    sweet::Data::Cart2D::DataSpectral &o_u,
                    sweet::Data::Cart2D::DataSpectral &o_v
            ) {

                sweet::Data::Cart2D::DataGrid u_phys(o_u.cart2DDataConfig);

                for (int j = 0; j < shackCart2DDataOps->space_res_physical[1]; j++) {
                    for (int i = 0; i < shackCart2DDataOps->space_res_physical[0]; i++) {

                        double x =
                                (i + 0.5) / shackCart2DDataOps->space_res_physical[0]; //*shackDict.sim.domain_size[0];
                        double y =
                                (j + 0.5) / shackCart2DDataOps->space_res_physical[1]; //*shackDict.sim.domain_size[1];

                        double x_dist = x - pos_x;
                        double y_dist = std::abs(y - pos_y);

                        // Create a circular column of constant height
                        if (y_dist < width / 2 && x_dist >= 0 && x_dist <= length) {
                            u_phys.grid_setValue(j, i, u_phys.grid_get(j, i) + speed);
                        }
                    }
                }

                o_u.loadCart2DDataGrid(u_phys);
                o_v.spectral_setZero();
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


                std::cout << "Generating impulse initial conditions.";

                setup_velocity(o_u, o_v);
                setup_depth(o_h);

                std::cout << "   Done! " << std::endl;

                return true;
            }


        };


    }
}

#endif //SWEET_PRECICE_COUPLING_IMPULSE_HPP
