/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_TIME_PDEADVECTIONCART2DTS_NA_SL_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_TIME_PDEADVECTIONCART2DTS_NA_SL_HPP

#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_BaseInterface.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>

#include <limits>
#include <sweet/Error/Base.hpp>
#include <sweet/SemiLagrangian/Cart2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "../benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp"


class PDEAdvectionCart2DTS_na_sl	:
		public PDEAdvectionCart2DTS_BaseInterface
{
	int timestepping_order;

	sweet::Data::Cart2D::DataSampler sampler2D;
	sweet::SemiLagrangian::Cart2D semiLagrangian;

	sweet::Data::Cart2D::DataSpectral prog_u_prev, prog_v_prev;


	/**
	 * Position in physical space given in longitude/latitude angles
	 */
	sweet::Data::Vector::Vector<double> posx_a, posy_a;

public:
	PDEAdvectionCart2DTS_na_sl();

	~PDEAdvectionCart2DTS_na_sl();

	bool setup(sweet::Data::Cart2D::Operators *io_ops);

private:
	void _setup();

public:
	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_phi,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);
};

#endif
