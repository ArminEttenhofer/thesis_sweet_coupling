/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *      
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <iostream>
#include <sweet/Error/Base.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>




class SimulationTestRK
{
public:
	sweet::Error::Base error;

	sweet::DEPRECATED_TimesteppingExplicitRKCart2DData timestepping;

	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::Error::Base error;

		sweet::Data::Cart2D::Config cart2DDataConfig;
		sweet::Data::Cart2D::Operators ops;

		sweet::Data::Cart2D::DataSpectral prog_h;
		sweet::Data::Cart2D::DataSpectral prog_u;
		sweet::Data::Cart2D::DataSpectral prog_v;

		sweet::Data::Cart2D::DataGrid prog_h_phys;
		sweet::Data::Cart2D::DataGrid prog_u_phys;
		sweet::Data::Cart2D::DataGrid prog_v_phys;


		bool setup(sweet::Data::Cart2D::Shack *i_shackCart2DDataOps)
		{
			/*
			 * Setup Cart2D Data Config & Operators
			 */
			if (!cart2DDataConfig.setupAuto(*i_shackCart2DDataOps))
				return error.forwardWithPositiveReturn(cart2DDataConfig.error);

			if (!ops.setup(
					cart2DDataConfig,
					*i_shackCart2DDataOps
				))
				return error.forwardWithPositiveReturn(ops.error);

			prog_h.setup(cart2DDataConfig);
			prog_u.setup(cart2DDataConfig);
			prog_v.setup(cart2DDataConfig);

			prog_h_phys.setup(cart2DDataConfig);
			prog_u_phys.setup(cart2DDataConfig);
			prog_v_phys.setup(cart2DDataConfig);

			return true;
		}

		void clear()
		{
			prog_h_phys.clear();
			prog_u_phys.clear();
			prog_v_phys.clear();

			prog_h.clear();
			prog_u.clear();
			prog_v.clear();

			ops.clear();
			cart2DDataConfig.clear();
		}
	};

	// Simulation data
	Data data;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::IO::Shack *shackIOData;

	int function_order = -1;
	int timestepping_order = -1;

public:
	SimulationTestRK(
			int i_argc,
			char *const * const i_argv,
			int i_function_order,
			int i_timestepping_order
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackCart2DDataOps(nullptr),
		shackTimestepControl(nullptr),
		function_order(i_function_order),
		timestepping_order(i_timestepping_order)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}


	bool setup()
	{
		shackProgArgDict.setup();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.printShackData();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

		data.setup(shackCart2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(data);

		data.prog_h_phys.grid_setValue(test_function(function_order, 0));
		data.prog_u_phys.grid_setValue(0);
		data.prog_v_phys.grid_setValue(0);

		data.prog_h.loadCart2DDataGrid(data.prog_h_phys);
		data.prog_u.loadCart2DDataGrid(data.prog_u_phys);
		data.prog_v.loadCart2DDataGrid(data.prog_v_phys);

		return true;
	}


	void clear()
	{
		data.clear();

		shackCart2DDataOps = nullptr;
		shackProgArgDict.clear();

		shackCart2DDataOps = nullptr;
		shackTimestepControl = nullptr;
	}


	/**
	 * Function of 4th order.
	 */
public:
	double test_function(
			int i_order,
			double z
	)
	{
		switch(i_order)
		{
		case 0:
			return	2.0;

		case 1:
			return	z
					+2.0;

		case 2:
			return	+(5./12.)*z*z
					+z
					+2.0;

		case 3:
			return	-(1./2.)*z*z*z
					+(5./12.)*z*z
					+z
					+2.0;

		case 4:
			return	(1./12.)*z*z*z*z
					-(1./2.)*z*z*z
					+(5./12.)*z*z
					+z
					+2.0;
		}

		SWEETErrorFatal("Not supported (i_order)");
		return 0;
	}



public:
	double test_function_diff(int i_order, double z)
	{
		switch(i_order)
		{
		case 0:
			return	0.0;

		case 1:
			return	1.0;

		case 2:
			return	+(10./12.)*z
					+1.0;

		case 3:
			return	-(3./2.)*z*z
					+(10./12.)*z
					+1.0;

		case 4:
			return	(4./12.)*z*z*z
					-(3./2.)*z*z
					+(10./12.)*z
					+1.0;
		}
		SWEETErrorFatal("Not supported test_function_diff (i_order)");
		return 0;
	}

	void p_run_euler_timestep_update(
			const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

			double i_current_timestamp = -1
	)
	{
		sweet::Data::Cart2D::DataGrid o_h_t_phys(data.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid o_u_t_phys(data.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid o_v_t_phys(data.cart2DDataConfig);

		o_h_t_phys.grid_setValue(test_function_diff(function_order, i_current_timestamp));
		o_u_t_phys.grid_setValue(0);
		o_v_t_phys.grid_setValue(0);

		o_h_t.loadCart2DDataGrid(o_h_t_phys);
		o_u_t.loadCart2DDataGrid(o_u_t_phys);
		o_v_t.loadCart2DDataGrid(o_v_t_phys);

		shackTimestepControl->currentTimestepNr++;
	}


	void runTimestep()
	{
		// either set time step size to 0 for autodetection or to
		// a positive value to use a fixed time step size
		SWEET_ASSERT(shackTimestepControl->currentTimestepSize > 0);

		shackTimestepControl->timestepHelperStart();

		timestepping.runTimestep(
				this,
				&SimulationTestRK::p_run_euler_timestep_update,	//!< pointer to function to compute euler time step updates
				data.prog_h, data.prog_u, data.prog_v,
				shackTimestepControl->currentTimestepSize,
				timestepping_order,
				shackTimestepControl->currentSimulationTime
			);

		shackTimestepControl->timestepHelperEnd();
	}


	bool should_quit()
	{
		return false;
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{
	for (int fun_order = 0; fun_order <= 4; fun_order++)
	{
		for (int timestepping_order = 1; timestepping_order <= 4; timestepping_order++)
		{
			SimulationTestRK simulation(i_argc, i_argv, fun_order, timestepping_order);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

			simulation.setup();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

			while(true)
			{
				if (simulation.shackTimestepControl->isFinalTimestepReached())
					break;

				simulation.runTimestep();


				double value_numerical = simulation.data.prog_h.toGrid().grid_get(0,0);
				double value_exact = simulation.test_function(simulation.function_order, simulation.shackTimestepControl->currentSimulationTime);

				std::cout << "t=" << simulation.shackTimestepControl->currentSimulationTime;
				std::cout << ", ";
				std::cout << "num=" << value_numerical;
				std::cout << ", ";
				std::cout << "exact=" << value_exact;
				std::cout << std::endl;
			}

			double value_numerical = simulation.data.prog_h.toGrid().grid_get(0,0);

			double value_exact = simulation.test_function(simulation.function_order, simulation.shackTimestepControl->currentSimulationTime);

			double error = std::abs(value_numerical - value_exact);

			if (fun_order <= timestepping_order)
			{
				if (error > 1e-8)
				{
					std::cout << "ERROR threshold exceeded!" << std::endl;
					return 1;
				}

				std::cout << "OK" << std::endl;
			}
			else
			{
				std::cout << "OK (errors expected)" << std::endl;

				if (error < 1e-8)
				{
					std::cout << "WARNING: Error expected to be larger, however relatively small error found" << std::endl;
					return 1;
				}

			}
		}
	}


	return 0;
}
