/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_BASEINTERFACE_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_BASEINTERFACE_HPP


#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/Shack.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "ShackTimeDiscretization.hpp"
#include "../Benchmarks/Shack.hpp"
#include "../Shack.hpp"



#if SWEET_PARAREAL || SWEET_XBRAID
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
#endif


class PDESWESphere2DTS_BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */

	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::ExpIntegration::Shack *shackExpIntegration;
	sweet::SemiLagrangian::Shack *shackTimesteppingSemiLagrangianSphere2DData;
	sweet::Parallelization::Shack *shackParallelization;

	ShackTimeDiscretization *shackPDESWETimeDisc;
	PDE_SWESphere2D::Benchmarks::Shack *shackPDESWEBenchmark;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;

	const sweet::Data::Sphere2D::Operators *ops;

	std::string timestepping_method;
	int timestepping_order;
	int timestepping_order2;

	sweet::Data::Sphere2D::DataGrid fg;

	PDESWESphere2DTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphere2DDataOps(nullptr),
		shackIOData(nullptr),
		shackExpIntegration(nullptr),
		shackTimesteppingSemiLagrangianSphere2DData(nullptr),
		shackParallelization(nullptr),
		shackPDESWETimeDisc(nullptr),
		shackPDESWEBenchmark(nullptr),
		shackPDESWESphere2D(nullptr),
		ops(nullptr),
		timestepping_order(-1),
		timestepping_order2(-1)
	{

	}

	virtual ~PDESWESphere2DTS_BaseInterface()
	{
	}


	bool setupFG()
	{
		SWEET_ASSERT(shackPDESWESphere2D != nullptr);
		if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
			fg = ops->getFG_fSphere2D(shackPDESWESphere2D->sphere2d_fsphere2d_f0);
		else
			fg = ops->getFG_rotatingSphere2D(shackPDESWESphere2D->sphere2d_rotating_coriolis_omega);

		return true;
	}


	virtual bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		shackIOData = io_shackDict->getAutoRegistration<sweet::IO::Shack>();
		shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ExpIntegration::Shack>();
		shackTimesteppingSemiLagrangianSphere2DData = io_shackDict->getAutoRegistration<sweet::SemiLagrangian::Shack>();
		shackParallelization = io_shackDict->getAutoRegistration<sweet::Parallelization::Shack>();

		shackPDESWETimeDisc = io_shackDict->getAutoRegistration<ShackTimeDiscretization>();
		shackPDESWEBenchmark = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();
		shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool shackRegistration(
			PDESWESphere2DTS_BaseInterface* i_baseInterface
	) 
	{
		*this = *i_baseInterface;
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}


	virtual bool setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		timestepping_method = i_timestepping_method;
		ops = io_ops;
		return true;
	}


public:
	virtual void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_h_pert,
			sweet::Data::Sphere2D::DataSpectral &io_u,
			sweet::Data::Sphere2D::DataSpectral &io_v,

			double i_dt,		//!< time step size
			double i_sim_timestamp
	) = 0;


	virtual bool implementsTimesteppingMethod(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string getIDString() = 0;

	virtual void printHelp()
	{
	}


	virtual void printImplementedTimesteppingMethods(
		std::ostream &o_ostream,
		const std::string &i_prefix
	)
	{
		o_ostream << i_prefix << "TODO" << std::endl;
	}


/*
 * Martin@Joao Please let's discuss this.
 */
// needed for parareal (instead of using directly shackDict.disc.timestepping_method)
#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE2D) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE2D)
	void runTimestep(
			sweet::DEPRECATED_pint::Parareal_GenericData* io_data,

			double i_dt,		//!< time step size
			double i_sim_timestamp
	)
	{
		sweet::Data::Sphere2D::DataSpectral h = *(io_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[0]);
		sweet::Data::Sphere2D::DataSpectral u = *(io_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[1]);
		sweet::Data::Sphere2D::DataSpectral v = *(io_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[2]);

		runTimestep(h, u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[0]) = h;
		*(io_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[1]) = u;
		*(io_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[2]) = v;

	}


	// for parareal SL
	virtual void set_previous_solution(
				sweet::Data::Sphere2D::DataSpectral &i_phi_prev,
				sweet::Data::Sphere2D::DataSpectral &i_vrt_prev,
				sweet::Data::Sphere2D::DataSpectral &i_div_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			sweet::DEPRECATED_pint::Parareal_GenericData* i_data
	)
	{
		sweet::Data::Sphere2D::DataSpectral phi_prev = *i_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[0];
		sweet::Data::Sphere2D::DataSpectral vrt_prev = *i_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[1];
		sweet::Data::Sphere2D::DataSpectral div_prev = *i_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[2];

		set_previous_solution(phi_prev, vrt_prev, div_prev);
	};

#endif
};

#endif
