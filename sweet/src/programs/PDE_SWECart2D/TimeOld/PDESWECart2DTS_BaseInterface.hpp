/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_PDESWECART2DTS_BASEINTERFACE_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_PDESWECART2DTS_BASEINTERFACE_HPP

#include <programs/PDE_SWECart2D/TimeOld/Shack.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/Shack.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>

#include "../Benchmarks/Shack.hpp"
#include "../Shack.hpp"


#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
#endif

class PDESWECart2DTS_BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */

	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::ExpIntegration::Shack *shackExpIntegration;
	PDE_SWECart2D::TimeDiscretization::Shack *shackPDESWETimeDisc;
	PDE_SWECart2D::Benchmarks::Shack *shackPDESWEBenchmark;
	PDE_SWECart2D::Shack *shackPDESWECart2D;
	sweet::SemiLagrangian::Shack *shackSemiLagrangian;

	sweet::Data::Cart2D::Operators *ops;

	PDESWECart2DTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackCart2DDataOps(nullptr),
		shackIOData(nullptr),
		shackExpIntegration(nullptr),
		shackPDESWETimeDisc(nullptr),
		shackPDESWEBenchmark(nullptr),
		shackPDESWECart2D(nullptr),
		shackSemiLagrangian(nullptr),
		ops(nullptr)
	{

	}

	virtual bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackIOData = io_shackDict->getAutoRegistration<sweet::IO::Shack>();
		shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ExpIntegration::Shack>();
		shackPDESWETimeDisc = io_shackDict->getAutoRegistration<PDE_SWECart2D::TimeDiscretization::Shack>();
		shackPDESWEBenchmark = io_shackDict->getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack>();
		shackPDESWECart2D = io_shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
		shackSemiLagrangian = io_shackDict->getAutoRegistration<sweet::SemiLagrangian::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


	virtual bool setup(
		sweet::Data::Cart2D::Operators *io_ops
	)
	{
		ops = io_ops;
		return true;
	}


public:
	virtual void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h_pert,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt,		//!< time step size
			double i_sim_timestamp
	) = 0;

#if (SWEET_PARAREAL && SWEET_PARAREAL_CART2D) || (SWEET_XBRAID && SWEET_XBRAID_CART2D)
	void runTimestep(
			sweet::DEPRECATED_pint::Parareal_GenericData* io_data,

			double i_dt,		//!< time step size
			double i_sim_timestamp
	)
	{
		sweet::Data::Cart2D::DataSpectral h_pert = *(io_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[0]);
		sweet::Data::Cart2D::DataSpectral u = *(io_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[1]);
		sweet::Data::Cart2D::DataSpectral v = *(io_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[2]);

		runTimestep(h_pert, u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[0]) = h_pert;
		*(io_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[1]) = u;
		*(io_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[2]) = v;

	}

	// for parareal SL
	virtual void set_previous_solution(
				sweet::Data::Cart2D::DataSpectral &i_h_prev,
				sweet::Data::Cart2D::DataSpectral &i_u_prev,
				sweet::Data::Cart2D::DataSpectral &i_v_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			sweet::DEPRECATED_pint::Parareal_GenericData* i_data
	)
	{
		sweet::Data::Cart2D::DataSpectral h_prev = *i_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[0];
		sweet::Data::Cart2D::DataSpectral u_prev = *i_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[1];
		sweet::Data::Cart2D::DataSpectral v_prev = *i_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[2];

		set_previous_solution(h_prev, u_prev, v_prev);
	};
#endif

	virtual ~PDESWECart2DTS_BaseInterface() {}
};

#endif
