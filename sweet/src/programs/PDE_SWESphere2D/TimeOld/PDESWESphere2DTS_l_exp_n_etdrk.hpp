/*
 * PDESWESphere2DTS_l_phi0_n_edt.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_N_ETDRK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_N_ETDRK_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_erk_n_erk.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"


class PDESWESphere2DTS_l_exp_n_etdrk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool shackRegistration(sweet::Shacks::Dictionary *io_shackDict) override;

	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			const std::string &i_exp_method,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestepSize,
			bool i_use_rexi_sphere2d_solver_preallocation
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;

	std::string getIDString() override;

private:
	PDESWESphere2DTS_l_erk_n_erk ts_l_erk_n_erk;

	PDESWESphere2DTS_l_exp ts_phi0_exp;
	PDESWESphere2DTS_l_exp ts_phi1_exp;
	PDESWESphere2DTS_l_exp ts_phi2_exp;

	PDESWESphere2DTS_l_exp ts_ups0_exp;
	PDESWESphere2DTS_l_exp ts_ups1_exp;
	PDESWESphere2DTS_l_exp ts_ups2_exp;
	PDESWESphere2DTS_l_exp ts_ups3_exp;



private:
	void euler_timestep_update_linear(
			const sweet::Data::Sphere2D::DataSpectral &i_h,
			const sweet::Data::Sphere2D::DataSpectral &i_u,
			const sweet::Data::Sphere2D::DataSpectral &i_v,

			sweet::Data::Sphere2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_v_t,	//!< time updates

			double i_max_timestamp
	);



private:
	void euler_timestep_update_nonlinear(
			const sweet::Data::Sphere2D::DataSpectral &i_h,
			const sweet::Data::Sphere2D::DataSpectral &i_u,
			const sweet::Data::Sphere2D::DataSpectral &i_v,

			sweet::Data::Sphere2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_v_t,	//!< time updates

			double i_max_timestamp
	);


public:
	PDESWESphere2DTS_l_exp_n_etdrk();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphere2DTS_l_exp_n_etdrk();
};

#endif
