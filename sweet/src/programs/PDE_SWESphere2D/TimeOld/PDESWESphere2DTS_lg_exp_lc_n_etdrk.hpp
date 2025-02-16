/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_N_ETDRK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_N_ETDRK_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"
#include "PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"


class PDESWESphere2DTS_lg_exp_lc_n_etdrk	: public PDESWESphere2DTS_BaseInterface
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
		int i_timestepping_order2,
		double i_timestepSize
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;


private:
	PDESWESphere2DTS_lg_erk_lc_n_erk ts_lg_erk_lc_n_erk;

	PDESWESphere2DTS_l_exp ts_phi0_rexi;
	PDESWESphere2DTS_l_exp ts_phi1_rexi;
	PDESWESphere2DTS_l_exp ts_phi2_rexi;

	PDESWESphere2DTS_l_exp ts_ups0_rexi;
	PDESWESphere2DTS_l_exp ts_ups1_rexi;
	PDESWESphere2DTS_l_exp ts_ups2_rexi;
	PDESWESphere2DTS_l_exp ts_ups3_rexi;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		ts_lg_erk_lc_n_erk.shackRegistration(io_shackDict);

		ts_phi0_rexi.shackRegistration(io_shackDict);
		ts_phi1_rexi.shackRegistration(io_shackDict);
		ts_phi2_rexi.shackRegistration(io_shackDict);

		ts_ups0_rexi.shackRegistration(io_shackDict);
		ts_ups1_rexi.shackRegistration(io_shackDict);
		ts_ups2_rexi.shackRegistration(io_shackDict);
		ts_ups3_rexi.shackRegistration(io_shackDict);
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
	PDESWESphere2DTS_lg_exp_lc_n_etdrk();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphere2DTS_lg_exp_lc_n_etdrk();
};

#endif
