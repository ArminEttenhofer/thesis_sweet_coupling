#ifndef PROGRAMS_ODE_GENERIC_XBRAID_APP_HPP
#define PROGRAMS_ODE_GENERIC_XBRAID_APP_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/App.hpp>
#include <programs/ODE_Generic/XBraid/DataContainer.hpp>
#include <programs/ODE_Generic/XBraid/TimeTree.hpp>

#include <programs/ODE_Generic/DE_Dahlquist/XBraid/DataContainer.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/XBraid/TimeTree.hpp>

namespace ODE_Generic {
namespace XBraid {

class App :
	public sweet::XBraid::App
{

public:
	int N = -1;

public:

	std::string ode = "";

	// General Shacks
	////////sweet::Data::Cart2D::Shack*					shackDataOps;

	// Level-defined shacks
	std::vector<ODE_Generic::Shack*>			shacksModel_levels;
	///std::vector<t_ShackTimeDiscretization*>		shacksTimeDisc_levels;
	////std::vector<sweet::Data::Cart2D::Shack*>		shacksDataOps_levels;

	// General config and operators
	/////sweet::Data::Cart2D::Config*			base_config = nullptr;
	/////sweet::Data::Cart2D::Operators*			base_op = nullptr;
	/////sweet::Data::Cart2DComplex::Operators*		base_op_complex = nullptr;

	// Level-defined config and operators
	/////std::vector<sweet::Data::Cart2D::Config*>		config_levels = {};
	/////std::vector<sweet::Data::Cart2D::Operators*>		op_levels = {};
	/////std::vector<sweet::Data::Cart2DComplex::Operators*>	op_complex_levels = {};

public:
	App( 		MPI_Comm				i_comm_t,
			int					i_rank,
			double					i_tstart,
			double					i_tstop,
			int 					i_ntime,
			std::string				i_ode
			/////sweet::Data::Cart2D::Config*		i_config,
			/////sweet::Data::Cart2D::Operators*		i_op,
			/////sweet::Data::Cart2DComplex::Operators*	i_op_complex
	)
		:
		sweet::XBraid::App(i_comm_t, i_rank, i_tstart, i_tstop, i_ntime)
	{
		ode = i_ode;
		if (ode == "dahlquist")
		{
			N = 1;
		}
		else
			SWEETErrorFatal("Unknown ODE '"+ode+"'");
		////base_config = i_config;
		////base_op = i_op;
		////base_op_complex = i_op_complex;
	}

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{

		sweet::XBraid::App::shackRegistration(io_shackDict);

		std::cout << "GLOBAL REGISTRATION - XBRAID SWECART2D" << std::endl;
		/////shackDataOps = io_shackDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();

		shackRegistrationLevels(&io_shackDict);

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(io_shackDict);

		return true;
	}

public:
	bool shackRegistrationLevels(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{

		sweet::XBraid::App::shackRegistrationLevels(io_shackDict);

		std::cout << "LOCAL REGISTRATION" << std::endl;
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			std::cout << " --> REGISTRATION LEVEL " << level << " - XBRAID SWECART2D " << std::endl;

			ODE_Generic::Shack* shackModel_level = nullptr;
			//////sweet::Data::Cart2D::Shack* shackDataOps_level = nullptr;

			shackModel_level = shacksDict_levels[level]->getAutoRegistration<ODE_Generic::Shack>();
			//////shackDataOps_level = shacksDict_levels[level]->getAutoRegistration<sweet::Data::Cart2D::Shack>();

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shacksDict_levels[level]);

			shacksModel_levels.push_back(shackModel_level);
			/////shacksDataOps_levels.push_back(shackDataOps_level);
		}

		return true;
	}



public:
	bool setupDataConfigOps() override
	{

		// set Config, Operators and spectral resolution for each level
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			if (ode == "dahlquist")
				///data_container.push_back(new ODE_Generic::DE_Dahlquist::XBraid::DataContainer);
				data_container.push_back(createNewVector(level));
			else
				SWEETErrorFatal("Unknown ODE '"+ode+"'");
			ODE_Generic::XBraid::DataContainer* data_container2 = (ODE_Generic::XBraid::DataContainer*) data_container.back();
			data_container2->setup(level);
		}


		//////// Set correct resolution in shacks
		//////for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		//////{
		//////	for (int j = 0; j < 2; j++)
		//////	{
		//////		shacksDataOps_levels[level]->space_res_physical[j] = config_levels[level]->grid_res[j];
		//////		shacksDataOps_levels[level]->space_res_spectral[j] = config_levels[level]->spectral_data_size[j];
		//////	}
		//////}

		// get buffer size
		////! To be updated depending on the tsm
		// Overestimated
		size_buffer = 0;
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
			size_buffer += N * sizeof(std::complex<double>);
			////size_buffer += N * config_levels[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);


		return true;

	}

	braid_Int
	Clone(
			braid_Vector	i_U,
			braid_Vector*	o_V
		) override
	{
		ODE_Generic::XBraid::DataContainer* U = (ODE_Generic::XBraid::DataContainer*) i_U;
		sweet::XBraid::Vector* v = createNewVector(U->level);
		ODE_Generic::XBraid::DataContainer* V = (ODE_Generic::XBraid::DataContainer*) v;

		*V = *U;
		*o_V = (braid_Vector) V;

		return 0;
	}


public:
	int getMaxLevel() override
	{
		return (int)data_container.size() - 1;
	}


public:
	int getBufferSize() override
	{
		////return N * config_levels[0]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
		return N * sizeof(std::complex<double>);
	}

public:
	sweet::XBraid::Vector* createNewVector(int i_level) override
	{
		sweet::XBraid::Vector* U;
		if (ode == "dahlquist")
			U = new ODE_Generic::DE_Dahlquist::XBraid::DataContainer;
		else
			SWEETErrorFatal("Unknown ODE '"+ode+"'");
		ODE_Generic::XBraid::DataContainer* U2 = (ODE_Generic::XBraid::DataContainer*) U;
		U2->setup(i_level);
		return U;
	}

	sweet::XBraid::TimeTree* createConfigureNewTimestepper(int i_level) override
	{

		// Set tsm and tso to instance of shackTimeDiscretization
		///////shacksTimeDisc_levels[level]->timestepping_method = tsms[i_level];

		sweet::XBraid::TimeTree* tsm;
		if (ode == "dahlquist")
			tsm = new ODE_Generic::DE_Dahlquist::XBraid::TimeTree;
		else
			SWEETErrorFatal("Unknown ODE '"+ode+"'");

		ODE_Generic::XBraid::TimeTree* tsm2 = (ODE_Generic::XBraid::TimeTree *) tsm;


		ODE_Generic::XBraid::DataContainer* U = (ODE_Generic::XBraid::DataContainer*) prog_initial_solution;

		tsm2->setup(
				tsms[i_level],
				shacksDict_levels[i_level],
				//////op_levels[i_level],
				//////op_complex_levels[i_level],
				U,
				shacksTimestepControl_levels[i_level]->currentTimestepSize
		);

		return tsm;
	}

};

}}

#endif
