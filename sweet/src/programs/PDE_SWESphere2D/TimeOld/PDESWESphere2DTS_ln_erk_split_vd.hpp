/*
 * PDESWESphere2DTS_split_lg_lc_na_nr_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LN_ERK_SPLIT_VD_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LN_ERK_SPLIT_VD_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"


class PDESWESphere2DTS_ln_erk_split_vd	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order,	//!< order of RK time stepping method
			bool i_lg,
			bool i_lc,
			bool i_na,
			bool i_nr,
			bool i_antialiasing_for_each_term
	);


public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;

		if (
				i_timestepping_method == "l_na_erk_split_vd"	||
				i_timestepping_method == "l_na_erk_split_aa_vd"	||
				i_timestepping_method == "l_erk_split_vd"	||
				i_timestepping_method == "l_erk_split_aa_vd"	||
				i_timestepping_method == "ln_erk_split_vd"		||
				i_timestepping_method == "ln_erk_split_aa_vd"
		)
			return true;

		return false;
	}

	std::string getIDString() override
	{
		return "ln_erk_split_vd";
	}

public:
	PDESWESphere2DTS_ln_erk_split_vd();

	virtual ~PDESWESphere2DTS_ln_erk_split_vd();

private:
	bool use_lg = false;
	bool use_lc = false;
	bool use_na = false;
	bool use_nr = false;

	bool anti_aliasing_for_each_term = false;


	// Sampler
	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk;


public:
	void euler_timestep_update_lg(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_pert_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_lc(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_lc_spectral_only(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_na(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_nr(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_set_tendencies(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_set_tendencies_na_only(
			const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_U_div,

			sweet::Data::Sphere2D::DataSpectral &o_U_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_vrt_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_U_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_U_phi,
			sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_U_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	void run_timestep_na(
			sweet::Data::Sphere2D::DataSpectral &io_U_phi,
			sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_U_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

};

#endif
