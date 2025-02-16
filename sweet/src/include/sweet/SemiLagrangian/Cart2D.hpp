/*
 * SemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef INCLUDE_SWEET_SEMILAGRANGIAN_CART2D_HPP
#define INCLUDE_SWEET_SEMILAGRANGIAN_CART2D_HPP

#include <sweet/Data/Cart2D/Convert/DataGrid_2_Vector.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>
#include <sweet/Data/Cart2D/Staggering.hpp>
#include <sweet/Data/Vector/Convert/Vector_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Tools/StopwatchBox.hpp>


namespace sweet {
namespace SemiLagrangian {

class Cart2D
{
	sweet::Data::Cart2D::DataSampler sample2D;
	const sweet::Data::Cart2D::Config *cart2DDataConfig;

public:
	Cart2D()	:
		cart2DDataConfig(nullptr)
	{
	}


	void setup(
		double i_domain_size[2],
		const sweet::Data::Cart2D::Config *i_cart2DDataConfig
	)
	{
		cart2DDataConfig = i_cart2DDataConfig;
		sample2D.setup(i_domain_size, cart2DDataConfig);
	}


	/**
	 * Stable extrapolation Two-Time-Level Scheme, Mariano Hortal,
	 *     Development and testing of a new two-time-level semi-lagrangian scheme (settls) in the ECMWF forecast model.
	 * Quaterly Journal of the Royal Meterological Society
	 *
	 * r_d = r_a - dt/2 * (2 * v_n(r_d) - v_{n-1}(r_d) + v_n(r_a))
	 *
	 * v^{iter} := (dt*v_n - dt*0.5*v_{n-1})
	 * r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
	 */
	void semi_lag_departure_points_settls(
			const sweet::Data::Cart2D::DataGrid &i_u_prev,	// Velocities at time t-1
			const sweet::Data::Cart2D::DataGrid &i_v_prev,
			const sweet::Data::Cart2D::DataGrid &i_u, 		// Velocities at time t
			const sweet::Data::Cart2D::DataGrid &i_v,

			const sweet::Data::Vector::Vector<double> &i_posx_a,	// Position of arrival points x / y
			const sweet::Data::Vector::Vector<double> &i_posy_a,

			double i_dt,				//!< time step size
			sweet::Data::Vector::Vector<double> &o_posx_d, 	//!< Position of departure points x / y
			sweet::Data::Vector::Vector<double> &o_posy_d,

			double i_domain_size[2],	//!< domain size

			const sweet::Data::Cart2D::Staggering *i_staggering,	//!< staggering, if any (ux, uy, vx, vy)

			int i_timestepping_order,

			int max_iters,
			double convergence_tolerance
	)
	{
#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().main_timestepping_semi_lagrangian.start();
#endif

		sweet::Data::Cart2D::Staggering s;
		if (i_staggering == nullptr)
		{
			s.setup_a_staggering();
			i_staggering = &(const sweet::Data::Cart2D::Staggering&)s;
		}

		std::size_t num_points = i_posx_a.numberOfElements;

		if (i_timestepping_order == 1)
		{
			o_posx_d = i_posx_a - i_dt*sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(i_u, false);
			o_posy_d = i_posy_a - i_dt*sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(i_v, false);
		}
		else if (i_timestepping_order == 2)
		{
			sweet::Data::Vector::Vector<double> u_prev = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(i_u_prev, false);
			sweet::Data::Vector::Vector<double> v_prev = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(i_v_prev, false);

			sweet::Data::Vector::Vector<double> u = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(i_u, false);
			sweet::Data::Vector::Vector<double> v = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(i_v, false);

			double dt = i_dt;

			sweet::Data::Cart2D::DataGrid u_extrap = sweet::Data::Vector::Convert::Vector_2_Cart2D_DataGrid::convert(2.0*u - u_prev, cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid v_extrap = sweet::Data::Vector::Convert::Vector_2_Cart2D_DataGrid::convert(2.0*v - v_prev, cart2DDataConfig);

			//Departure point tmp
			sweet::Data::Vector::Vector<double> rx_d_new(num_points);
			sweet::Data::Vector::Vector<double> ry_d_new(num_points);

			//Previous departure point
			sweet::Data::Vector::Vector<double> rx_d_prev = i_posx_a;
			sweet::Data::Vector::Vector<double> ry_d_prev = i_posy_a;

			// initialize departure points with arrival points
			o_posx_d = i_posx_a;
			o_posy_d = i_posy_a;

			int iters = 0;
			double diff = 999;
			for (; iters < max_iters; iters++)
			{
				// r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
				rx_d_new = i_posx_a - dt*0.5 *
					(u + sample2D.bilinear_scalar(
						u_extrap,
						o_posx_d, o_posy_d, i_staggering->u[0], i_staggering->u[1]
				));
				ry_d_new = i_posy_a - dt*0.5 *
					(v + sample2D.bilinear_scalar(
						v_extrap,
						o_posx_d, o_posy_d, i_staggering->v[0], i_staggering->v[1]
				));

				diff = (rx_d_new - rx_d_prev).reduce_maxAbs()/i_domain_size[0] + (ry_d_new - ry_d_prev).reduce_maxAbs()/i_domain_size[1];
				rx_d_prev = rx_d_new;
				ry_d_prev = ry_d_new;

				SWEET_THREADING_SPACE_PARALLEL_FOR
				for (std::size_t i = 0; i < num_points; i++)
				{
					o_posx_d.data[i] = sample2D.wrapPeriodic(rx_d_new.data[i], sample2D.domain_size[0]);
					o_posy_d.data[i] = sample2D.wrapPeriodic(ry_d_new.data[i], sample2D.domain_size[1]);
				}

				if (diff < convergence_tolerance)
				   break;
			}


			if (convergence_tolerance > 0)
			{
				if (diff > convergence_tolerance)
				{
					std::cout << "WARNING: Over convergence tolerance" << std::endl;
					std::cout << "+ Iterations: " << iters << std::endl;
					std::cout << "+ maxAbs: " << diff << std::endl;
					std::cout << "+ Convergence tolerance: " << convergence_tolerance << std::endl;
				}
			}
		}
		else
		{
			SWEETErrorFatal("This time integration order is not implemented");
		}


#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().main_timestepping_semi_lagrangian.stop();
#endif

	}
};

}}

#endif
