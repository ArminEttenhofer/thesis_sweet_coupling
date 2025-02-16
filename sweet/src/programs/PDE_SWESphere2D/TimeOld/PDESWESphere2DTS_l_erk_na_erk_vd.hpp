/*
 * PDESWESphere2DTS_l_erk_na_erk_vd.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_ERK_NA_ERK_VD_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_ERK_NA_ERK_VD_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_ln_erk_split_vd.hpp"



class PDESWESphere2DTS_l_erk_na_erk_vd	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	PDESWESphere2DTS_ln_erk_split_vd l_erk_split_vd;
	PDESWESphere2DTS_ln_erk_split_vd na_erk_split_vd;

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		l_erk_split_vd.shackRegistration(io_shackDict);
		na_erk_split_vd.shackRegistration(io_shackDict);
		return true;
	}

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order,	//!< order of RK time stepping method
			int i_order2
	);

public:
	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	PDESWESphere2DTS_l_erk_na_erk_vd();


	virtual ~PDESWESphere2DTS_l_erk_na_erk_vd();
};

#endif
