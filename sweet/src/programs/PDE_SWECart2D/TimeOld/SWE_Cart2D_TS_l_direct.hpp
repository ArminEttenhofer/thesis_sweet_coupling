/*
 * SWE_Cart2D_TS_l_direct.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_DIRECT_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_DIRECT_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <programs/PDE_SWECart2D/TimeOld/Shack.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/GridMapping.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Data/Cart2D/Staggering.hpp>
#include <limits>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/ExpIntegration/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "../Shack.hpp"



class SWE_Cart2D_TS_l_direct	:
		public PDESWECart2DTS_BaseInterface
{
	typedef double T;

	sweet::ExpIntegration::ExpFunction<T> expFunctions;

	sweet::Data::Cart2D::GridMapping cart2DDataGridMapping;

	// Precompute Z = Z(k1, k2, Dt) = Q * e^{\Lambda*Dt} * Q^{-1}
	std::vector<std::vector<std::array<std::array<std::complex<T>, 3>, 3>>> Z;  // Z[k1][k2][0,1,2][0,1,2];
	double dt_precompute_phin = 0.; // store dt for which Z has been precomputed in order to check if it is necessary to recompute it.


public:

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_cgrid(
			sweet::Data::Cart2D::DataSpectral &io_h_pert,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,		//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,		//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid_cart2ddata(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);



	void run_timestep_agrid_cart2ddatacomplex(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	bool setup(
			sweet::Data::Cart2D::Operators *io_ops,
			const std::string &i_function_name
	);

	bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	);

	virtual ~SWE_Cart2D_TS_l_direct() {}
};

#endif
