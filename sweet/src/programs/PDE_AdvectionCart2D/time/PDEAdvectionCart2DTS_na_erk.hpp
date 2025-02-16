/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_TIME_PDEADVECTIONCART2DTS_NA_ERK_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_TIME_PDEADVECTIONCART2DTS_NA_ERK_HPP

#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_BaseInterface.hpp>
#include <programs/PDE_AdvectionCart2D/time/ShackPDEAdvectionCart2DTimeDiscretization.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp"


class PDEAdvectionCart2DTS_na_erk	:
		public PDEAdvectionCart2DTS_BaseInterface
{
public:
	int timestepping_order;

	sweet::DEPRECATED_TimesteppingExplicitRKCart2DData timestepping_rk;


public:
	PDEAdvectionCart2DTS_na_erk();

	~PDEAdvectionCart2DTS_na_erk();

	bool setup(sweet::Data::Cart2D::Operators *io_ops);

public:
	void euler_timestep_update(
			const sweet::Data::Cart2D::DataSpectral &i_phi,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_vort,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_div,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_phi_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_vort_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);

public:
	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_phi,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_vort,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_div,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);
};

#endif
