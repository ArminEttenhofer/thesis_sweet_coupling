/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_ERK_N_ERK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_ERK_N_ERK_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"



class PDESWESphere2DTS_l_erk_n_erk	: public PDESWESphere2DTS_BaseInterface
{
public:
	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk_linear;
	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk_nonlinear;

	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order,	//!< order of RK time stepping method
			int i_order2
	);

	PDESWESphere2DTS_l_erk_n_erk();

	virtual ~PDESWESphere2DTS_l_erk_n_erk();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

public:
	void euler_timestep_update_linear(
			const sweet::Data::Sphere2D::DataSpectral &i_h,
			const sweet::Data::Sphere2D::DataSpectral &i_u,
			const sweet::Data::Sphere2D::DataSpectral &i_v,

			sweet::Data::Sphere2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_v_t,	//!< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			const sweet::Data::Sphere2D::DataSpectral &i_h,
			const sweet::Data::Sphere2D::DataSpectral &i_u,
			const sweet::Data::Sphere2D::DataSpectral &i_v,

			sweet::Data::Sphere2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_v_t,	//!< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt,
			double i_simulation_timestamp
	);

public:
	void euler_timestep_update(
			const sweet::Data::Sphere2D::DataSpectral &i_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_vort,
			const sweet::Data::Sphere2D::DataSpectral &i_div,

			sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);

};

#endif
