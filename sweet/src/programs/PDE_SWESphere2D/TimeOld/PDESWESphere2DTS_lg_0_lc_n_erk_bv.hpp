/*
 * Author: Pedor Peixoto <ppeixoto@usp.br>
 * based on stuff from:
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *         
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_0_LC_N_ERK_BV_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_0_LC_N_ERK_BV_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"
#include "PDESWESphere2DTS_lg_irk.hpp"



class PDESWESphere2DTS_lg_0_lc_n_erk_bv	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order	//!< order of RK time stepping method
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		if (
			i_timestepping_method == "lg_0_lc_n_erk_bv" 
		)
			return true;

		return false;
	}

	std::string getIDString() override
	{
		std::string s = "lg_0_lc_n_erk_bv";
		return s;
	}

	void printHelp() override;


private:
	int timestepping_order;

	double timestep_size;

	/*
	 * Non-linear time steppers
	 */
	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk;

public:
	PDESWESphere2DTS_lg_0_lc_n_erk_bv();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void euler_timestep_update(
		const sweet::Data::Sphere2D::DataSpectral &i_phi, //prog
		const sweet::Data::Sphere2D::DataSpectral &i_vrt, //prog
		const sweet::Data::Sphere2D::DataSpectral &i_div, //prog

		sweet::Data::Sphere2D::DataSpectral &o_phi_t, //updated with euler
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t, //updated with euler
		sweet::Data::Sphere2D::DataSpectral &o_div_t, //updated with euler

		double i_simulation_timestamp
	);

	virtual ~PDESWESphere2DTS_lg_0_lc_n_erk_bv();
};

#endif
