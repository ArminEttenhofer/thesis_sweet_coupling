#ifndef PROGRAMS_PDE_SWESPHERE2D_XBRAID_APP_HPP
#define PROGRAMS_PDE_SWESPHERE2D_XBRAID_APP_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/App.hpp>
#include <programs/PDE_SWESphere2D/XBraid/DataContainer.hpp>
#include <programs/PDE_SWESphere2D/XBraid/TimeTree.hpp>

namespace PDE_SWESphere2D {
namespace XBraid {

class App :
	public sweet::XBraid::App
{

public:
	const int N = 3;

public:

	// General Shacks
	sweet::Data::Sphere2D::Shack*					shackDataOps;

	// Level-defined shacks
	std::vector<PDE_SWESphere2D::Shack*>				shacksModel_levels;
	std::vector<sweet::Data::Sphere2D::Shack*>			shacksDataOps_levels;

	// General config and operators
	sweet::Data::Sphere2D::Config*					base_config = nullptr;
	sweet::Data::Sphere2D::Operators*				base_op = nullptr;
	sweet::Data::Sphere2DComplex::Operators*			base_op_complex = nullptr;

	// Level-defined config and operators
	std::vector<sweet::Data::Sphere2D::Config*>			config_levels = {};
	std::vector<sweet::Data::Sphere2D::Operators*>			op_levels = {};
	std::vector<sweet::Data::Sphere2DComplex::Operators*>		op_complex_levels = {};

public:
	App( 		MPI_Comm				i_comm_t,
			int					i_rank,
			double					i_tstart,
			double					i_tstop,
			int 					i_ntime,
			sweet::Data::Sphere2D::Config*		i_config,
			sweet::Data::Sphere2D::Operators*		i_op,
			sweet::Data::Sphere2DComplex::Operators*	i_op_complex
	)
		:
		sweet::XBraid::App(i_comm_t, i_rank, i_tstart, i_tstop, i_ntime)
	{
		base_config = i_config;
		base_op = i_op;
		base_op_complex = i_op_complex;
	}

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{

		sweet::XBraid::App::shackRegistration(io_shackDict);

		std::cout << "GLOBAL REGISTRATION - XBRAID SWESPHERE2D" << std::endl;
		shackDataOps = io_shackDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();

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
			std::cout << " --> REGISTRATION LEVEL " << level << " - XBRAID SWESPHERE2D " << std::endl;

			PDE_SWESphere2D::Shack* shackModel_level = nullptr;
			sweet::Data::Sphere2D::Shack* shackDataOps_level = nullptr;

			shackModel_level = shacksDict_levels[level]->getAutoRegistration<PDE_SWESphere2D::Shack>();
			shackDataOps_level = shacksDict_levels[level]->getAutoRegistration<sweet::Data::Sphere2D::Shack>();

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shacksDict_levels[level]);

			shacksModel_levels.push_back(shackModel_level);
			shacksDataOps_levels.push_back(shackDataOps_level);
		}

		return true;
	}



public:
	bool setupDataConfigOps() override
	{

		// set Config, Operators, data container and spectral resolution for each level
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			data_container.push_back(new PDE_SWESphere2D::XBraid::DataContainer);
			if (shackXBraid->xbraid_spatial_coarsening)
			{
				int ndim = base_config->ndim;
				int N_physical[ndim];
				int N_spectral[ndim];
				for (int j = 0; j < ndim; j++)
				{
					N_physical[j] = -1;
					// proportional to time step
					if (shackXBraid->xbraid_spatial_coarsening == 1)
						N_spectral[j] = std::max(4, int(shackDataOps->space_res_spectral[j] / std::pow(shackXBraid->xbraid_cfactor, level)));
					else if (shackXBraid->xbraid_spatial_coarsening > 1)
					{
						if (level == 0)
							N_spectral[j] = std::max(4, shackDataOps->space_res_spectral[j]);
						else
							N_spectral[j] = std::max(4, shackXBraid->xbraid_spatial_coarsening);
					}
					else
						SWEETErrorFatal("Invalid parameter xbraid_spatial_coarsening");
				}

				config_levels.push_back(new sweet::Data::Sphere2D::Config);
				config_levels.back()->setupAuto(shacksDataOps_levels[level]);
				op_levels.push_back(new sweet::Data::Sphere2D::Operators(config_levels.back(), shacksDataOps_levels[level]));
				op_complex_levels.push_back(new sweet::Data::Sphere2DComplex::Operators(config_levels.back(), shacksDataOps_levels[level]));

				std::cout << "Spectral resolution at level " << level << " : ";
				for (int j = 0; j < ndim; j++)
					std::cout << N_spectral[j] << " ";
				std::cout << std::endl;
			}
			else
			{
				config_levels.push_back(base_config);
				op_levels.push_back(base_op);
				op_complex_levels.push_back(base_op_complex);
			}
			PDE_SWESphere2D::XBraid::DataContainer* data_container2 = (PDE_SWESphere2D::XBraid::DataContainer*) data_container.back();
			data_container2->setup(config_levels.back(), level);
		}


		// Set correct resolution in shacks
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			shacksDataOps_levels[level]->space_res_physical[0] = config_levels[level]->grid_num_lon;
			shacksDataOps_levels[level]->space_res_physical[1] = config_levels[level]->grid_num_lat;

			shacksDataOps_levels[level]->space_res_spectral[0] = config_levels[level]->spectral_modes_m_max + 1;
			shacksDataOps_levels[level]->space_res_spectral[1] = config_levels[level]->spectral_modes_n_max + 1;
		}

		// get buffer size
		////! To be updated depending on the tsm
		// Overestimated
		size_buffer = 0;
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
			size_buffer += N * config_levels[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);


		return true;

	}

	braid_Int
	Clone(
			braid_Vector	i_U,
			braid_Vector*	o_V
		) override
	{
		PDE_SWESphere2D::XBraid::DataContainer* U = (PDE_SWESphere2D::XBraid::DataContainer*) i_U;
		sweet::XBraid::Vector* v = createNewVector(U->level);
		PDE_SWESphere2D::XBraid::DataContainer* V = (PDE_SWESphere2D::XBraid::DataContainer*) v;

		*V = *U;
		*o_V = (braid_Vector) V;

		return 0;
	}


public:
	int getMaxLevel() override
	{
		return (int)this->data_container.size() - 1;
	}


public:
	int getBufferSize() override
	{
		return N * config_levels[0]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
	}

public:
	sweet::XBraid::Vector* createNewVector(int i_level) override
	{
		sweet::XBraid::Vector* U = new PDE_SWESphere2D::XBraid::DataContainer;
		PDE_SWESphere2D::XBraid::DataContainer* U2 = (PDE_SWESphere2D::XBraid::DataContainer*) U;
		U2->setup(config_levels[i_level], i_level);
		return U;
	}

	sweet::XBraid::TimeTree* createConfigureNewTimestepper(int i_level) override
	{

		// Set tsm and tso to instance of shackTimeDiscretization

		sweet::XBraid::TimeTree* tsm = new PDE_SWESphere2D::XBraid::TimeTree;
		PDE_SWESphere2D::XBraid::TimeTree* tsm2 = (PDE_SWESphere2D::XBraid::TimeTree *) tsm;
		PDE_SWESphere2D::XBraid::DataContainer* U = (PDE_SWESphere2D::XBraid::DataContainer*) prog_initial_solution;

		tsm2->setup(
				tsms[i_level],
				shacksDict_levels[i_level],
				op_levels[i_level],
				op_complex_levels[i_level],
				U,
				shacksTimestepControl_levels[i_level]->currentTimestepSize
		);

		return tsm;
	}

};

}}

#endif
