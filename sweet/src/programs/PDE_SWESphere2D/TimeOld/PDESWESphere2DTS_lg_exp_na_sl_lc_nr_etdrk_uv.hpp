/*
 * PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv
 *
 * Created on: 24 Mar 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_NA_SL_LC_NR_ETDRK_UV_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_NA_SL_LC_NR_ETDRK_UV_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"
#include "PDESWESphere2DTS_ln_erk_split_uv.hpp"


class PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;
	void printHelp() override;

private:
	PDESWESphere2DTS_ln_erk_split_uv ts_ln_erk_split_uv;


private:
	enum NLRemainderTreatment_enum{
		NL_REMAINDER_IGNORE,
		NL_REMAINDER_NONLINEAR,
	};

	NLRemainderTreatment_enum nonlinear_remainder_treatment;

public:
	sweet::Data::Sphere2D::DataSpectral U_phi_prev, U_vrt_prev, U_div_prev;

	PDESWESphere2DTS_l_exp ts_phi0_exp;
	PDESWESphere2DTS_l_exp ts_phi2_exp;

	PDESWESphere2DTS_l_exp ts_psi1_exp;
	PDESWESphere2DTS_l_exp ts_psi2_exp;

	sweet::SemiLagrangian::Sphere2D semiLagrangian;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		ts_ln_erk_split_uv.shackRegistration(io_shackDict);

		ts_phi0_exp.shackRegistration(io_shackDict);
		ts_phi2_exp.shackRegistration(io_shackDict);

		ts_psi1_exp.shackRegistration(io_shackDict);
		ts_psi2_exp.shackRegistration(io_shackDict);
		return true;
	}



public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestepSize,

			NLRemainderTreatment_enum i_nonlinear_remainder_treatment
	);


public:
	PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE2D) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE2D)
	void set_previous_solution(
				sweet::Data::Sphere2D::DataSpectral &i_phi_prev,
				sweet::Data::Sphere2D::DataSpectral &i_vrt_prev,
				sweet::Data::Sphere2D::DataSpectral &i_div_prev
	) override
	{
		U_phi_prev = i_phi_prev;
		U_vrt_prev = i_vrt_prev;
		U_div_prev = i_div_prev;
	}
#endif

	virtual ~PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv();
};

#endif
