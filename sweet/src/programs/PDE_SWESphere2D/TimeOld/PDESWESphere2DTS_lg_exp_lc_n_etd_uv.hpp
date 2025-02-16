/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_N_ETD_UV_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_N_ETD_UV_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"
#include "PDESWESphere2DTS_ln_erk_split_uv.hpp"


class PDESWESphere2DTS_lg_exp_lc_n_etd_uv	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			int i_timestepping_order,
			double i_timestepSize,
			bool i_with_na,
			bool i_with_nr
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;
	void printHelp() override;


	PDESWESphere2DTS_ln_erk_split_uv ts_ln_erk_split_uv;

	sweet::Data::Sphere2D::DataSpectral NU_phi_prev, NU_vrt_prev, NU_div_prev;
	sweet::Data::Sphere2D::DataSpectral NU_phi_prev_2, NU_vrt_prev_2, NU_div_prev_2;

	bool with_na;
	bool with_nr;

	PDESWESphere2DTS_l_exp ts_phi0_exp;
	PDESWESphere2DTS_l_exp ts_phi1_exp;
	PDESWESphere2DTS_l_exp ts_phi2_exp;
	PDESWESphere2DTS_l_exp ts_phi3_exp;

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		ts_phi0_exp.shackRegistration(io_shackDict);
		ts_phi1_exp.shackRegistration(io_shackDict);
		ts_phi2_exp.shackRegistration(io_shackDict);
		ts_phi3_exp.shackRegistration(io_shackDict);
		return true;
	}


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
	PDESWESphere2DTS_lg_exp_lc_n_etd_uv();


	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphere2DTS_lg_exp_lc_n_etd_uv();
};

#endif
