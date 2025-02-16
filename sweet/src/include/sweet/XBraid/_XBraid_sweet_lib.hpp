/*
 * Xbraid_sweet_lib.hpp
 *
 *  Created on: 10 Jun 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef INCLUDE_SWEET_XBRAID_XBRAID_SWEET_LIB_HPP
#define INCLUDE_SWEET_XBRAID_XBRAID_SWEET_LIB_HPP

#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
#include <sweet/_DEPRECATED_pint/PInT_Common.hpp>

#include <xbraid/braid.hpp>

#if SWEET_GUI
#include<sweet/GUI/VisSweet.hpp>
#endif

#include <algorithm>

#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/XBraid/Shack.hpp>

#include <sweet/Data/GenericContainer/ConfigBase.hpp>

#include "GeometryDependentDefinitions.hpp"


//////////#if SWEET_XBRAID_SCALAR
//////////	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Scalar.hpp>
//////////	#include <sweet/Data/Scalar/Shack.hpp>
//////////	#include <sweet/Data/Scalar/Config.hpp>
//////////	#include <sweet/Data/Scalar/Operators.hpp>
//////////	#include <programs/ODE_Scalar/ODEScalarTimeSteppers.hpp>
//////////	#include <programs/ODE_Scalar/ShackODEScalar.hpp>
//////////	#include <programs/ODE_Scalar/time/ShackODEScalarTimeDiscretization.hpp>
//////////	#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/TimeTreeIR.hpp>
//////////
//////////namespace sweet {
//////////namespace XBraid {
//////////
//////////	typedef ODEScalarTimeSteppers t_tsmType;
//////////	typedef ODE_Generic::DE_Dahlquist::TimeTree t_tsmType_NewTS;
//////////	typedef ODE_Generic::DE_Dahlquist::DataContainer::Simulation t_dataContainer;
//////////	typedef ShackODEScalar t_ShackModel;
//////////	typedef ShackODEScalarTimeDiscretization t_ShackTimeDiscretization;
//////////	typedef ShackODEScalarBenchmarks t_ShackBenchmarks;
//////////	typedef sweet::Data::Scalar::Shack t_ShackDataOps;
//////////	typedef sweet::Data::Scalar::Config	t_Config;
//////////	typedef sweet::Data::Scalar::Operators	t_Operators;
//////////	typedef sweet::Data::Scalar::Operators	t_OperatorsComplex;
//////////	#define N_vec 1
//////////
//////////}}
//////////
//////////#elif SWEET_XBRAID_CART2D
//////////	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Cart2DData_Spectral.hpp>
//////////	#include <sweet/Data/Cart2D/Shack.hpp>
//////////	#if SWEET_XBRAID_CART2D_BURGERS
//////////		#include <programs/pde_burgersCart2D/PDEBurgersCart2D_TimeSteppers.hpp>
//////////		#include <programs/pde_burgersCart2D/ShackPDEBurgersCart2D.hpp>
//////////		#include <programs/pde_burgersCart2D/time/ShackPDEBurgersCart2DTimeDiscretization.hpp>
//////////
//////////namespace sweet {
//////////namespace XBraid {
//////////
//////////		typedef BurgersCart2DTimeSteppers t_tsmType;
//////////		typedef BurgersCart2DTimeSteppersNewTS t_tsmType_NewTS;
//////////		typedef ShackPDEBurgersCart2D t_ShackModel;
//////////		typedef ShackPDEBurgersCart2DTimeDiscretization t_ShackTimeDiscretization;
//////////		typedef ShackPDEBurgersCart2DBenchmarks t_ShackBenchmarks;
//////////}}
//////////
//////////		#define N_vec 2
//////////	#elif SWEET_XBRAID_CART2D_SWE
//////////		#include <programs/PDE_SWECart2D/TimeSteppers.hpp>
//////////		#include <programs/PDE_SWECart2D/BenchmarksCombined.hpp>
//////////		#include <programs/PDE_SWECart2D/Shack.hpp>
//////////		#include <programs/PDE_SWECart2D/TimeOld/Shack.hpp>
//////////
//////////namespace sweet {
//////////namespace XBraid {
//////////
//////////		typedef PDE_SWECart2D::TimeSteppers t_tsmType;
//////////		typedef PDE_SWECart2D::TimeSteppersNewTS t_tsmType_NewTS;
//////////		typedef PDE_SWECart2D::Shack t_ShackModel;
//////////		typedef PDE_SWECart2D::TimeDiscretization::Shack t_ShackTimeDiscretization;
//////////		typedef PDE_SWECart2D::Benchmarks::Shack t_ShackBenchmarks;
//////////		typedef sweet::Data::Cart2D::Config	t_Config;
//////////		typedef sweet::Data::Cart2D::Operators	t_Operators;
//////////}}
//////////
//////////		#define N_vec 3
//////////	#endif
//////////#elif SWEET_XBRAID_SPHERE2D
//////////	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Sphere2DData_Spectral.hpp>
//////////	#include <sweet/Data/Sphere2D/Shack.hpp>
//////////	#include <programs/PDE_SWESphere2D/TimeOld/PDESWESphere2D_TimeSteppers.hpp>
//////////	#include <programs/PDE_SWESphere2D/TimeTree/TimeTree.hpp>
//////////	#include <programs/PDE_SWESphere2D/BenchmarksCombined.hpp>
//////////	#include <programs/PDE_SWESphere2D/Shack.hpp>
//////////	#include <programs/PDE_SWESphere2D/TimeOld/ShackTimeDiscretization.hpp>
//////////	#include <programs/PDE_SWESphere2D/DataContainer.hpp>
//////////
//////////namespace sweet {
//////////namespace XBraid {
//////////	typedef PDESWESphere_TimeSteppers t_tsmType;
//////////	typedef PDE_SWESphere::TimeTree::TimeTree t_tsmType_NewTS;
//////////	typedef PDE_SWESphere::Shack t_ShackModel;
//////////	typedef ShackTimeDiscretization t_ShackTimeDiscretization;
//////////	typedef PDE_SWESphere::Benchmarks::Shack t_ShackBenchmarks;
//////////	typedef PDE_SWESphere::DataContainer t_dataContainer;
//////////	typedef sweet::Data::Sphere2D::Config	t_Config;
//////////	typedef sweet::Data::Sphere2D::Operators	t_Operators;
//////////}}
//////////	#define N_vec 3
//////////#endif
//////////
//////////namespace sweet {
//////////namespace XBraid {
//////////
//////////#if SWEET_XBRAID_CART2D
//////////	// Grid Mapping (staggered grid)
//////////	sweet::Data::Cart2D::GridMapping gridMapping;
//////////#endif
//////////
//////////}}



namespace sweet {
namespace XBraid {


class sweet_BraidVector;
class sweet_BraidApp;



////////* --------------------------------------------------------------------
/////// * XBraid vector 
/////// * Stores the state of the simulation for a given time step
/////// * Define BraidVector, can contain anything, and be named anything
/////// * --> Put all time-dependent information here
/////// * -------------------------------------------------------------------- */
///////class sweet_BraidVector
///////{
///////public:
///////	sweet::DEPRECATED_pint::Parareal_GenericData*	data = nullptr;
///////
///////	t_dataContainer*		data_DE = nullptr;
///////	t_dataContainer*		data_DE_tmp = nullptr;
///////
///////	t_Config*		config = nullptr;
///////
//////////////#if SWEET_XBRAID_SCALAR
//////////////	ODE_Scalar::DataContainer*		data_DE = nullptr;
//////////////	ODE_Scalar::DataContainer*		data_DE_tmp = nullptr;
//////////////#elif SWEET_XBRAID_CART2D
//////////////	sweet::Data::Cart2D::PlaneData_Config* cart2DDataConfig;
//////////////	#if SWEET_XBRAID_PLANE_SWE
//////////////		PDE_SWECart2D::DataContainer*		data_DE = nullptr;
//////////////		PDE_SWECart2D::DataContainer*		data_DE_tmp = nullptr;
//////////////	#endif
//////////////#elif SWEET_XBRAID_SPHERE2D
//////////////	sweet::Data::Sphere2D::Config* sphere2DDataConfig;
//////////////	PDE_SWESphere2D::DataContainer*		data_DE = nullptr;
//////////////	PDE_SWESphere2D::DataContainer*		data_DE_tmp = nullptr;
//////////////#endif
///////
///////	int level;
///////
///////
///////	sweet_BraidVector(
///////				t_Config* i_config,
////////////#if SWEET_XBRAID_CART2D
////////////				sweet::Data::Cart2D::Config* i_cart2DDataConfig,
////////////#elif SWEET_XBRAID_SPHERE2D
////////////				sweet::Data::Sphere2D::Config* i_sphere2DDataConfig,
////////////#endif
///////				int i_level
///////	)
///////		:
///////		config(i_config),
////////////#if SWEET_XBRAID_CART2D
////////////		cart2DDataConfig(i_cart2DDataConfig),
////////////#elif SWEET_XBRAID_SPHERE2D
////////////		sphere2DDataConfig(i_sphere2DDataConfig),
////////////#endif
///////		level(i_level)
///////	{
///////		allocate_data();
///////	}
///////
///////	virtual ~sweet_BraidVector()
///////	{
///////		if (data)
///////		{
///////			delete data;
///////			data = nullptr;
///////		}
///////
///////		if (data_DE)
///////		{
///////			data_DE->clear();
///////			///delete data_DE;
///////			data_DE = nullptr;
///////		}
///////
///////		if (data_DE_tmp)
///////		{
///////			data_DE_tmp->clear();
///////			///delete data_DE_tmp;
///////			data_DE_tmp = nullptr;
///////		}
///////
///////
///////	}
///////
///////	sweet_BraidVector(const sweet_BraidVector &i_vector)
///////	{
///////		config = i_vector.config;
////////////#if SWEET_XBRAID_CART2D
////////////		cart2DDataConfig = i_vector.cart2DDataConfig;
////////////#elif SWEET_XBRAID_SPHERE2D
////////////		sphere2DDataConfig = i_vector.sphere2DDataConfig;
////////////#endif
///////		*data = *i_vector.data;
///////		///*data_DE = *i_vector.data_DE;
///////		data_DE->op_setVector(*i_vector.data_DE);
///////		level = i_vector.level;
///////	};
///////
///////	sweet_BraidVector& operator=(const sweet_BraidVector &i_vector)
///////	{
///////		config = i_vector.config;
/////////////#if SWEET_XBRAID_CART2D
/////////////		cart2DDataConfig = i_vector.cart2DDataConfig;
/////////////#elif SWEET_XBRAID_SPHERE2D
/////////////		sphere2DDataConfig = i_vector.sphere2DDataConfig;
/////////////#endif
///////		*data = *i_vector.data;
///////		///*data_DE = *i_vector.data_DE;
///////		data_DE->op_setVector(*i_vector.data_DE);
///////		level = i_vector.level;
///////		return *this;
///////	};
///////
///////	sweet_BraidVector operator+(
///////			const sweet_BraidVector &i_vector
///////	)	const
///////	{
///////		sweet_BraidVector out(config, level);
/////////////#if SWEET_XBRAID_SCALAR
/////////////		sweet_BraidVector out(level);
/////////////#elif SWEET_XBRAID_CART2D
/////////////		sweet_BraidVector out(cart2DDataConfig, level);
/////////////#elif SWEET_XBRAID_SPHERE2D
/////////////		sweet_BraidVector out(sphere2DDataConfig, level);
/////////////#endif
///////		*out.data = *data;
///////		*out.data += *i_vector.data;
///////
///////		////*out.data_DE = *this->data_DE;
///////		out.data_DE->op_setVectorPlusVector(*this->data_DE, *i_vector.data_DE);
///////
///////		return out;
///////	}
///////
///////	sweet_BraidVector operator*(
///////			const double i_value
///////	)	const
///////	{
///////		sweet_BraidVector out(config, level);
/////////////#if SWEET_XBRAID_SCALAR
/////////////		sweet_BraidVector out(level);
/////////////#elif SWEET_XBRAID_CART2D
/////////////		sweet_BraidVector out(cart2DDataConfig, level);
/////////////#elif SWEET_XBRAID_SPHERE2D
/////////////		sweet_BraidVector out(sphere2DDataConfig, level);
/////////////#endif
///////		*out.data = *data;
///////		*out.data *= i_value;
///////
///////		///*out.data_DE = *this->data_DE;
///////		///*out.data_DE *= i_value;
///////		out.data_DE->op_setVector(*this->data_DE);
///////		out.data_DE->op_mulScalar(i_value);
///////
///////		return out;
///////	}
///////
///////
///////	void allocate_data()
///////	{
///////
///////		data_DE = new t_dataContainer;
///////		data_DE_tmp = new t_dataContainer;
///////
///////		data_DE->setup(this->config);
///////		data_DE_tmp->setup(this->config);
///////
/////////////#if SWEET_XBRAID_SCALAR
/////////////		{
/////////////			data = new sweet::DEPRECATED_pint::Parareal_GenericData_Scalar<N_vec>;
/////////////			data->allocate_data();
/////////////
/////////////			data_DE = new ODE_Generic::DataContainer;
/////////////			data_DE->setup(this->config);
/////////////			data_DE_tmp->setup(this->config);
/////////////		}
/////////////
/////////////#elif SWEET_XBRAID_CART2D
/////////////		{
/////////////			data = new sweet::DEPRECATED_pint::Parareal_GenericData_Cart2DData_Spectral<N_vec>;
/////////////			data->setup_data_config(cart2DDataConfig);
/////////////			data->allocate_data();
/////////////
/////////////			data_DE = new PDE_SWECart2D::DataContainer;
/////////////			data_DE->setup(this->config);
/////////////			data_DE_tmp->setup(this->config);
/////////////		}
/////////////
/////////////#elif SWEET_XBRAID_SPHERE2D
/////////////		{
/////////////			data = new sweet::DEPRECATED_pint::Parareal_GenericData_Sphere2DData_Spectral<N_vec>;
/////////////			data->setup_data_config(sphere2DDataConfig);
/////////////			data->allocate_data();
/////////////
/////////////			data_DE = new PDE_SWESphere::DataContainer;
/////////////			data_DE->setup(this->config);
/////////////			data_DE_tmp->setup(this->config);
/////////////		}
/////////////#endif
///////	}
///////
///////
///////
///////};
///////
///////
////////**
/////// * operator to support operations such as:
/////// *
/////// * 1.5 * arrayData;
/////// *
/////// * Otherwise, we'd have to write it as arrayData*1.5
/////// *
/////// */
///////inline
///////static
///////sweet_BraidVector operator*(
///////		const double i_value,
///////		const sweet_BraidVector &i_vector
///////)
///////{
///////	return i_vector * i_value;
///////}



// Wrapper for BRAID's App object·
// --> Put all time INDEPENDENT information here
class sweet_BraidApp
			: public BraidApp, sweet::DEPRECATED_pint::PInT_Common
{

public:

	sweet::Shacks::Dictionary*			shackDict;
	sweet::XBraid::Shack*			shackXBraid;
	sweet::Parallelization::Shack*		shackParallelization;
	sweet::TimeTree::Shack*		shackTimestepControl;
	t_ShackBenchmarks*			shackBenchmarks;
	t_ShackDataOps*			shackDataOps;

	// Shacks corresponding to each level
	std::vector<t_ShackModel*>			shacksModel_levels;
	std::vector<t_ShackTimeDiscretization*>		shacksTimeDisc_levels;
	///std::vector<t_ShackBenchmarks*>		shacksBenchmarks_levels;
	std::vector<sweet::Shacks::Dictionary*>		shacksDict_levels;
	std::vector<sweet::TimeTree::Shack*>	shacksTimestepControl_levels;

	std::vector<t_ShackDataOps*>			shacksDataOps_levels;
/////#if SWEET_XBRAID_CART2D
/////	std::vector<sweet::Data::Cart2D::Shack*>		shacksCart2DDataOps_levels;
/////#elif SWEET_XBRAID_SPHERE2D
/////	std::vector<sweet::Data::Sphere2D::Shack*>		shacksSphere2DDataOps_levels;
/////#endif

	double				dt;
	std::vector<t_tsmType*>		timeSteppers;
	std::vector<t_tsmType_NewTS*>	timeSteppers_NewTS;
	std::vector<t_dataContainer*>	data_container;

#if SWEET_GUI
	// single vector to replace vectors above
	std::vector<SimulationGuiCallBacks*> levels_simulations;
#endif


	int			size_buffer;		// overestimated

	int rank;

	// Reference solutions (online error computation)
	std::vector<sweet_BraidVector*> xbraid_data_ref_exact;
	std::vector<sweet_BraidVector*> xbraid_data_fine_exact;

	// Solution from previous timestep (for SL)
	std::vector<std::vector<sweet_BraidVector*>> sol_prev; // sol_prev[level][timestep]
	std::vector<std::vector<int>> sol_prev_iter; // store iteration in which the solution has been stored
	std::vector<int> first_timeid_level; // store first timestep (time id) in this level and processor
	std::vector<int> last_timeid_level; // store last timestep (time id) in this level and processor

	// Timestepping method and orders for each level
	std::vector<std::string> tsms;
	std::vector<int> tsos;
	std::vector<int> tsos2;
	std::vector<bool> is_SL;
	bool contains_SL = false;

	// Viscosity for each level
	std::vector<int> viscosity_orders;
	std::vector<double> viscosity_coefficients;

	// Smaller spectral resolution among levels
	int min_spectral_size = INT_MAX;

	// Custom time grid
	std::vector<double> custom_time_steps = {};

	// Effective number of levels
	int nlevels = -1;

	// Use new timesteppers (...,...,...)
	bool useNewTimeSteppers = false;

	t_Config*			base_config = nullptr;
	t_Operators*			base_op = nullptr;
	t_OperatorsComplex*		base_op_complex = nullptr;

	std::vector<t_Config*>			config_levels = {};
	std::vector<t_Operators*>		op_levels = {};
	std::vector<t_OperatorsComplex*>	op_complex_levels = {};


	t_dataContainer		prog_initial_solution;
	t_dataContainer		prog_initial_guess;

public:

	// Constructor·
	sweet_BraidApp( MPI_Comm		i_comm_t,
			int			i_rank,
			double			i_tstart,
			double			i_tstop,
			int 			i_ntime,
			t_Config* i_config,
			t_Operators* i_op,
			t_OperatorsComplex* i_op_complex
/////#if SWEET_XBRAID_CART2D
/////			,
/////			sweet::Data::Cart2D::Config* i_cart2DDataConfig,
/////			sweet::Data::Cart2D::Operators* i_op_cart2d
/////#elif SWEET_XBRAID_SPHERE2D
/////			,
/////			sweet::Data::Sphere2D::Config* i_sphere2DDataConfig,
/////			sweet::Data::Sphere2D::Operators* i_op_sphere2D
/////			sweet::Data::Sphere2D::SphereOperatorsComplex* i_op_sphere2D_complex
/////#endif
			)
		:
			BraidApp(i_comm_t, i_tstart, i_tstop, i_ntime),
			rank(i_rank)
	{

			base_config = i_config;
			base_op = i_op;

			// get min_spectral size
			for (std::size_t i = 0; i < config_levels.size(); i++)
				min_spectral_size = std::min(	min_spectral_size,
									std::min((int)config_levels[i]->spectral_data_size[0], (int)config_levels[i]->spectral_data_size[1])
								);


////////#if SWEET_XBRAID_CART2D
////////			base_cart2DDataConfig = i_cart2DDataConfig;
////////			base_op_cart2d = i_op_cart2d;
////////
////////			// get min_spectral size
////////			for (std::size_t i = 0; i < cart2DDataConfig.size(); i++)
////////				min_spectral_size = std::min(	min_spectral_size,
////////									std::min((int)cart2DDataConfig[i]->spectral_data_size[0], (int)cart2DDataConfig[i]->spectral_data_size[1])
////////								);
////////
////////#elif SWEET_XBRAID_SPHERE2D
////////			base_sphere2DDataConfig = i_sphere2DDataConfig;
////////			base_op_sphere2D = i_op_sphere2D;
////////			base_op_sphere_complex = i_op_sphere2D_complex;
////////
////////			// get min_spectral size
////////			for (std::size_t i = 0; i < sphere2DDataConfig.size(); i++)
////////				min_spectral_size = std::min(	min_spectral_size,
////////									std::min(sphere2DDataConfig[i]->spectral_modes_m_max, sphere2DDataConfig[i]->spectral_modes_n_max)
////////								);
////////#endif
	}

	virtual ~sweet_BraidApp()
	{

		for (std::vector<t_tsmType*>::iterator it = timeSteppers.begin();
							it != timeSteppers.end();
							it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<t_tsmType_NewTS*>::iterator it = timeSteppers_NewTS.begin();
							it != timeSteppers_NewTS.end();
							it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<t_dataContainer*>::iterator it = data_container.begin();
							it != data_container.end();
							it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<sweet_BraidVector*>::iterator it = xbraid_data_ref_exact.begin();
								it != xbraid_data_ref_exact.end();
								it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<sweet_BraidVector*>::iterator it = xbraid_data_fine_exact.begin();
								it != xbraid_data_fine_exact.end();
								it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}

		for (std::vector<std::vector<sweet_BraidVector*>>::iterator it = sol_prev.begin();
										it != sol_prev.end();
										it++)
			for (std::vector<sweet_BraidVector*>::iterator it2 = it->begin();
										it2 != it->end();
										it2++)
				if (*it2)
				{
					delete *it2;
					*it2 = nullptr;
				}

#if SWEET_GUI
		for (std::vector<sweet::SimulationGuiCallbacks*>::iterator it =
				levels_simulations.begin();
				it != levels_simulations.end();
			it++)
			if (*it)
			{
				delete *it;
				*it = nullptr;
			}
#endif

	}

	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		std::cout << "GLOBAL REGISTRATION" << std::endl;
		shackDict = &io_shackDict;
		shackIOData = io_shackDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimestepControl = io_shackDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackBenchmarks = io_shackDict.getAutoRegistration<t_ShackBenchmarks>();
		shackXBraid = io_shackDict.getAutoRegistration<sweet::XBraid::Shack>();
		shackDataOps = io_shackDict.getAutoRegistration<t_ShackDataOps>();
///////#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
///////		shackCart2DDataOps = io_shackDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
///////#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
///////		shackSphere2DDataOps = io_shackDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
///////		shackPDESWESphere2D = io_shackDict.getAutoRegistration<PDE_SWESphere2D::Shack>();
///////#endif

		shackRegistrationLevels(&io_shackDict);

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(io_shackDict);

		return true;
	}

public:
	bool shackRegistrationLevels(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{

		std::cout << "LOCAL REGISTRATION" << std::endl;
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			std::cout << " --> REGISTRATION LEVEL " << level << std::endl;
			sweet::Shacks::Dictionary* shackDict_level = new sweet::Shacks::Dictionary;
			*shackDict_level = *io_shackDict;

			sweet::TimeTree::Shack* shackTimestepControl_level = nullptr;
			t_ShackTimeDiscretization* shackTimeDisc_level = nullptr;
			///t_ShackBenchmarks* shackBenchmark_level = nullptr;
			t_ShackModel* shackModel_level = nullptr;
			t_ShackDataOps* shackDataOps_level = nullptr;
/////#if SWEET_XBRAID_CART2D
/////			sweet::Data::Cart2D::Shack* shackCart2DDataOps_level = nullptr;
/////#elif SWEET_XBRAID_SPHERE2D
/////			sweet::Data::Sphere2D::Shack* shackSphere2DDataOps_level = nullptr;
/////#endif

			shackTimestepControl_level = shackDict_level->getAutoRegistration<sweet::TimeTree::Shack>();
			shackTimeDisc_level = shackDict_level->getAutoRegistration<t_ShackTimeDiscretization>();
			////shackBenchmark_level = shackDict_level->getAutoRegistration<t_ShackBenchmarks>();
			shackModel_level = shackDict_level->getAutoRegistration<t_ShackModel>();
			shackDataOps_level = shackDict_level->getAutoRegistration<t_ShackDataOps>();
#if SWEET_XBRAID_CART2D
			shackCart2DDataOps_level = shackDict_level->getAutoRegistration<sweet::Data::Cart2D::Shack>();
#elif SWEET_XBRAID_SPHERE2D
			shackSphere2DDataOps_level = shackDict_level->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
#endif


			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackDict_level);

			shacksDict_levels.push_back(shackDict_level);
			shacksTimestepControl_levels.push_back(shackTimestepControl_level);
			shacksTimeDisc_levels.push_back(shackTimeDisc_level);
			shacksModel_levels.push_back(shackModel_level);
			shacksDataOps_levels.push_back(shackDataOps_level);
//////#if SWEET_XBRAID_CART2D
//////			shacksCart2DDataOps_levels.push_back(shackCart2DDataOps_level);
//////#elif SWEET_XBRAID_SPHERE2D
//////			shacksSphere2DDataOps_levels.push_back(shackSphere2DDataOps_level);
//////#endif
		}

		return true;
	}


public:
	void setup(
			BraidCore& i_core,
			ODE_Generic::DE_Dahlquist::DataContainer::Simulation i_prog_initial_solution,
			ODE_Generic::DE_Dahlquist::DataContainer::Simulation i_prog_initial_guess
	)
	{

		/////////////////////////////////////////////////
		// get parameters from simVars and set to Core //
		/////////////////////////////////////////////////

		i_core.SetMaxLevels(shackXBraid->xbraid_max_levels);
		/////i_core.SetIncrMaxLevels();

		i_core.SetSkip(shackXBraid->xbraid_skip);

		i_core.SetMinCoarse(shackXBraid->xbraid_min_coarse);

		///i_core.SetRelaxOnlyCG(shackXBraid->xbraid_relax_only_cg);

		i_core.SetNRelax(-1, shackXBraid->xbraid_nrelax);
		if (shackXBraid->xbraid_nrelax0 > -1)
			i_core.SetNRelax(0, shackXBraid->xbraid_nrelax0);

		i_core.SetAbsTol(shackXBraid->xbraid_tol);
		i_core.SetRelTol(shackXBraid->xbraid_tol);

		i_core.SetTemporalNorm(shackXBraid->xbraid_tnorm);

		i_core.SetCFactor(-1, shackXBraid->xbraid_cfactor);
		if (shackXBraid->xbraid_cfactor0 > -1)
			i_core.SetCFactor(0, shackXBraid->xbraid_cfactor0);

		///i_core.SetPeriodic(shackXBraid->xbraid_periodic);

		////i_core.SetResidual();

		i_core.SetMaxIter(shackXBraid->xbraid_max_iter);

		i_core.SetPrintLevel(shackXBraid->xbraid_print_level);

		i_core.SetSeqSoln(shackXBraid->xbraid_use_seq_soln);

		i_core.SetAccessLevel(shackXBraid->xbraid_access_level);

		i_core.SetNFMG(shackXBraid->xbraid_fmg);
		if (shackXBraid->xbraid_fmg)
			i_core.SetFMG();

		i_core.SetNFMGVcyc(shackXBraid->xbraid_fmg_vcyc);

		i_core.SetStorage(shackXBraid->xbraid_storage);

		//i_core.SetRevertedRanks(shackXBraid->xbraid_reverted_ranks);

		////i_core.SetRefine(shackXBraid->xbraid_refine);
		///i_core.SetMaxRefinements(shackXBraid->xbraid_max_Refinements);

		i_core.SetTimeGrid(sweet_BraidApp::sweet_TimeGrid);


		prog_initial_solution = i_prog_initial_solution;
		prog_initial_guess = i_prog_initial_guess;
		////shackTimestepControl = i_shackTimestepControl;
		////shackBenchmarks = i_shackBenchmarks;

		///setup_timesteppers();
		setup();
	}






public:
	/*
	 * Setup timesteppers for each level.
	 * IMPORTANT: this function must be called after setting up initial conditions, since the benchmark may modify simulation parameters
	 *            call it in the first call of Step function.
	 */
	void setup_timesteppers()
	{

		////////////////////////////////////
		// get tsm and tso for each level //
		////////////////////////////////////
		tsms = getLevelParameterFromParameters<std::string>("timestepping_method");
		if (! this->useNewTimeSteppers)
		{
			tsos = getLevelParameterFromParameters<int>("timestepping_order", 1);
			tsos2 = getLevelParameterFromParameters<int>("timestepping_order", 2);
		}
		viscosity_orders = getLevelParameterFromParameters<int>("viscosity_order");
		viscosity_coefficients = getLevelParameterFromParameters<double>("viscosity_coefficient");

		// create a timeSteppers instance for each level
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{


			// Set tsm and tso to instance of shackTimeDiscretization
			shacksTimeDisc_levels[level]->timestepping_method = tsms[level];
			if (!this->useNewTimeSteppers)
			{
				shacksTimeDisc_levels[level]->timestepping_order = tsos[level];
				shacksTimeDisc_levels[level]->timestepping_order2 = tsos2[level];
			}

			// Configure timesteppers with the correct timestep for this level
			shacksTimestepControl_levels[level]->printShack();
			shacksTimestepControl_levels[level]->currentTimestepSize = shacksTimestepControl_levels[level]->currentTimestepSize *
												std::pow(shackXBraid->xbraid_cfactor, level);
			shacksTimestepControl_levels[level]->printShack();

			if (rank == 0)
				std::cout << "Timestep size at level " << level << " : " << shacksTimestepControl_levels[level]->currentTimestepSize << std::endl;

			if (this->useNewTimeSteppers)
			{
				t_tsmType_NewTS* tsm = new t_tsmType_NewTS;

				tsm->setup_1_registerAllTimesteppers();
				tsm->setup_2_shackRegistration(this->shacksDict_levels[level]);
				ERROR_FORWARD(*tsm);

				sweet::Shacks::ProgramArgumentsDictionary* shackProgArgDict_level = (sweet::Shacks::ProgramArgumentsDictionary*) this->shacksDict_levels[level];
				tsm->setup_3_timestepper(
								this->tsms[level],
								shackProgArgDict_level,
								this->op_levels[level],
								this->op_complex_levels[level],
								*this->data_container[level]
							);

				ERROR_FORWARD(*tsm);

				tsm->timeIntegrator->setTimeStepSize(this->shacksTimestepControl_levels[level]->currentTimestepSize);

				ERROR_FORWARD(*tsm);

				this->timeSteppers_NewTS.push_back(tsm);
			}


///////////////#if SWEET_XBRAID_SCALAR
///////////////			ODEScalarTimeSteppers* tsm = new ODEScalarTimeSteppers;
///////////////
///////////////			tsm->shackRegistration(*shacksDict_levels[level]);
///////////////			ERROR_FORWARD(*tsm);
///////////////
///////////////			tsm->shackTimeDisc = shacksTimeDisc_levels[level];
///////////////			tsm->setup(
///////////////					*shacksDict_levels[level]
///////////////				);
///////////////
///////////////			////tsm->setup(
///////////////			////		*simVars_levels[level]
///////////////			////	);
///////////////#elif SWEET_XBRAID_CART2D
///////////////	#if SWEET_XBRAID_CART2D_SWE
///////////////			PDE_SWECart2D::TimeSteppers* tsm = new PDE_SWECart2D::TimeSteppers;
///////////////
///////////////			tsm->shackRegistration(*shacksDict_levels[level]);
///////////////			ERROR_FORWARD(*tsm);
///////////////
///////////////			tsm->setup(
///////////////					///tsms[level],
///////////////					///tsos[level],
///////////////					///tsos2[level],
///////////////					op_cart2d[level],
///////////////					shacksDict_levels[level]
///////////////				);
///////////////	#elif SWEET_XBRAID_CART2D_BURGERS
///////////////			Burgers_Cart2D_TimeSteppers* tsm = new Burgers_Cart2D_TimeSteppers;
///////////////
///////////////			tsm->shackRegistration(*shacksDict_levels[level]);
///////////////			ERROR_FORWARD(*tsm);
///////////////
///////////////			tsm->setup(
///////////////					///tsms[level],
///////////////					///tsos[level],
///////////////					///tsos2[level],
///////////////					op_cart2d[level],
///////////////					shacksDict_levels[level]
///////////////				);
///////////////	#endif
///////////////#elif SWEET_XBRAID_SPHERE2D
///////////////
///////////////			if (this->useNewTimeSteppers)
///////////////			{
///////////////				t_tsmType_NewTS* tsm = new t_tsmType_NewTS;
///////////////
///////////////				tsm->setup_1_registerAllTimesteppers();
///////////////				tsm->setup_2_shackRegistration(this->shacksDict_levels[level]);
///////////////				ERROR_FORWARD(*tsm);
///////////////
///////////////				sweet::Shacks::ShackProgArgDictionary* shackProgArgDict_level = (sweet::Shacks::ShackProgArgDictionary*) this->shacksDict_levels[level];
///////////////				tsm->setup_3_timestepper(
///////////////								this->tsms[level],
///////////////								///this->shacksDict_levels[level],
///////////////								shackProgArgDict_level,
///////////////								this->op_sphere[level],
///////////////								this->op_sphere_complex[level],
///////////////								*this->data_container[level]
///////////////							);
///////////////
///////////////				ERROR_FORWARD(*tsm);
///////////////
///////////////				tsm->timeIntegrator->setTimeStepSize(this->shacksTimestepControl_levels[level]->current_timestepSize);
///////////////
///////////////				ERROR_FORWARD(*tsm);
///////////////
///////////////				this->timeSteppers_NewTS.push_back(tsm);
///////////////			}
///////////////			else
///////////////			{
///////////////				PDESWESphere_TimeSteppers* tsm = new PDESWESphere_TimeSteppers;
///////////////
///////////////				tsm->setup_1_registerAllTimesteppers();
///////////////
///////////////				tsm->setup_2_shackRegistration(this->shacksDict_levels[level]);
///////////////				ERROR_FORWARD(*tsm);
///////////////
///////////////				tsm->setup_3_timestepper(
///////////////								this->tsms[level],
///////////////								this->shacksDict_levels[level],
///////////////								this->op_sphere[level]
///////////////							);
///////////////
///////////////				this->timeSteppers.push_back(tsm);
///////////////			}
///////////////#endif

			// check if tsm is SL
			if ( std::find(SL_tsm.begin(), SL_tsm.end(), tsms[level]) == SL_tsm.end())
				is_SL.push_back(false);
			else
			{
				is_SL.push_back(true);
				contains_SL = true; // requires extra communication
////#if SWEET_XBRAID_CART2D
////				actual_size_buffer += N * cart2DDataConfig[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
////#elif SWEET_XBRAID_SPHERE2D
////				actual_size_buffer += N * sphere2DDataConfig[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
////#endif
			}
		}

	}

public:
	void setup()
	{

		PInT_Common::setup(
					shackIOData,
					shackDataOps
//////#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
//////					,
//////					shackCart2DDataOps
//////#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
//////					,
//////					shackSphere2DDataOps
//////#endif
		);


		// create vectors for storing solutions from previous timestep (SL)
		for (int i = 0; i < shackXBraid->xbraid_max_levels; i++)
		{
			std::vector<sweet_BraidVector*> v = {};
			sol_prev.push_back(v);

			std::vector<int> w = {};
			sol_prev_iter.push_back(w);

			first_timeid_level.push_back(INT_MAX);
			last_timeid_level.push_back(-1);
		}


		// set custom time grid
		double t = 0;
		while (t < shackTimestepControl->maxSimulationTime - 1e-10)
		{
			double dt = shackTimestepControl->currentTimestepSize;
			double dt2;
			if ( t + dt < shackTimestepControl->maxSimulationTime - 1e-10)
				dt2 = dt;
			else
				dt2 = shackTimestepControl->maxSimulationTime - t;
			custom_time_steps.push_back(dt2);
			t += dt2;
			////std::cout << "TIME STEP " << dt2 << std::endl;
		}




		// set spectral resolution and Cart2D/Sphere2D ops for each level
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			if (shackXBraid->xbraid_spatial_coarsening)
			{
				int N_physical[2] = {-1, -1};
				int N_spectral[2];
				for (int j = 0; j < 2; j++)
				{
					// proportional to time step
					if (shackXBraid->xbraid_spatial_coarsening == 1)
						N_spectral[j] = std::max(4, int(shackDataOps->space_res_spectral[j] / std::pow(shackXBraid->xbraid_cfactor, level)));
	////#if SWEET_XBRAID_CART2D
	////					N_spectral[j] = std::max(4, int(shackCart2DDataOps->space_res_spectral[j] / std::pow(shackXBraid->xbraid_cfactor, level)));
	////#else
	////					N_spectral[j] = std::max(4, int(shackSphere2DDataOps->space_res_spectral[j] / std::pow(shackXBraid->xbraid_cfactor, level)));
	////#endif
					else if (shackXBraid->xbraid_spatial_coarsening > 1)
					{
						if (level == 0)
							N_spectral[j] = std::max(4, shackDataOps->space_res_spectral[j]);
	////#if SWEET_XBRAID_CART2D
	////						N_spectral[j] = std::max(4, shackCart2DDataOps->space_res_spectral[j]);
	////#else
	////						N_spectral[j] = std::max(4, shackSphere2DDataOps->space_res_spectral[j]);
	////#endif
						else
							N_spectral[j] = std::max(4, shackXBraid->xbraid_spatial_coarsening);
					}
					else
						SWEETErrorFatal("Invalid parameter xbraid_spatial_coarsening");
				}

				config_levels.push_back(new t_Config);
				config_levels.back()->setupAuto(shacksDataOps_levels[level]);
				op_levels.push_back(new t_Operators(config_levels.back(), shacksDataOps_levels[level]));
				op_complex_levels.push_back(new t_OperatorsComplex(config_levels.back(), shacksDataOps_levels[level]));

				data_container.back()->setup(config_levels.back());

	//////#if SWEET_XBRAID_CART2D
	//////			cart2DDataConfig.push_back(new sweet::Data::Cart2D::PlaneData_Config);
	//////			cart2DDataConfig.back()->setupAuto(N_physical, N_spectral, shackCart2DDataOps->reuse_spectral_transformation_plans);
	//////			op_cart2d.push_back(new sweet::Data::Cart2D::PlaneOperators(cart2DDataConfig.back(), shacksCart2DDataOps_levels[level]));

	//////			data_container.back()->setup(planeDataConfig.back());
	//////#else
	//////			sphereDataConfig.push_back(new sweet::Data::Sphere2D::SphereData_Config);
	//////			sphereDataConfig.back()->setupAuto(shacksSphere2DDataOps_levels[level]);
	//////			op_sphere.push_back(new sweet::Data::Sphere2D::SphereOperators(sphere2DDataConfig.back(), shacksSphere2DDataOps_levels[level]));
	//////			op_sphere_complex.push_back(new sweet::Data::Sphere2D::SphereOperatorsComplex(sphere2DDataConfig.back(), shacksSphere2DDataOps_levels[level]));

	//////			data_container.back()->setup(sphereDataConfig.back());

	//////#endif

				std::cout << "Spectral resolution at level " << level << " : " << N_spectral[0] << " " << N_spectral[1] << std::endl;
			}
			else
			{
				config_levels.push_back(base_config);
				op_levels.push_back(base_op);
				op_complex_levels.push_back(base_op_complex);
	////#if SWEET_XBRAID_CART2D
	////			cart2DDataConfig.push_back(base_cart2DDataConfig);
	////			op_cart2d.push_back(base_op_cart2d);
	////#else
	////			sphere2DDataConfig.push_back(base_sphere2DDataConfig);
	////			op_sphere2D.push_back(base_op_sphere2D);
	////			op_sphere2D_complex.push_back(base_op_sphere2D_complex);
	////#endif
			}
		}



///////#if SWEET_XBRAID_SCALAR
///////		data_container.back()->setup();
///////#elif SWEET_XBRAID_CART2D || SWEET_XBRAID_SPHERE2D
///////		// set spectral resolution and Cart2D/Sphere2D ops for each level
///////		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
///////		{
///////			if (shackXBraid->xbraid_spatial_coarsening)
///////			{
///////				int N_physical[2] = {-1, -1};
///////				int N_spectral[2];
///////				for (int j = 0; j < 2; j++)
///////				{
///////					// proportional to time step
///////					if (shackXBraid->xbraid_spatial_coarsening == 1)
///////	#if SWEET_XBRAID_CART2D
///////						N_spectral[j] = std::max(4, int(shackCart2DDataOps->space_res_spectral[j] / std::pow(shackXBraid->xbraid_cfactor, level)));
///////	#else
///////						N_spectral[j] = std::max(4, int(shackSphere2DDataOps->space_res_spectral[j] / std::pow(shackXBraid->xbraid_cfactor, level)));
///////	#endif
///////					else if (shackXBraid->xbraid_spatial_coarsening > 1)
///////					{
///////						if (level == 0)
///////	#if SWEET_XBRAID_CART2D
///////							N_spectral[j] = std::max(4, shackCart2DDataOps->space_res_spectral[j]);
///////	#else
///////							N_spectral[j] = std::max(4, shackSphere2DDataOps->space_res_spectral[j]);
///////	#endif
///////						else
///////							N_spectral[j] = std::max(4, shackXBraid->xbraid_spatial_coarsening);
///////					}
///////					else
///////						SWEETErrorFatal("Invalid parameter xbraid_spatial_coarsening");
///////				}
///////	#if SWEET_XBRAID_CART2D
///////				cart2DDataConfig.push_back(new sweet::Data::Cart2D::PlaneData_Config);
///////				cart2DDataConfig.back()->setupAuto(N_physical, N_spectral, shackCart2DDataOps->reuse_spectral_transformation_plans);
///////				op_cart2d.push_back(new sweet::Data::Cart2D::PlaneOperators(cart2DDataConfig.back(), shacksCart2DDataOps_levels[level]));
///////
///////				data_container.back()->setup(planeDataConfig.back());
///////	#else
///////				sphereDataConfig.push_back(new sweet::Data::Sphere2D::SphereData_Config);
///////				///sphereDataConfig.back()->setupAuto(N_physical, N_spectral, shackSphereDataOps->reuse_spectral_transformation_plans);
///////				sphereDataConfig.back()->setupAuto(shacksSphere2DDataOps_levels[level]);
///////				op_sphere.push_back(new sweet::Data::Sphere2D::SphereOperators(sphere2DDataConfig.back(), shacksSphere2DDataOps_levels[level]));
///////				op_sphere_complex.push_back(new sweet::Data::Sphere2D::SphereOperatorsComplex(sphere2DDataConfig.back(), shacksSphere2DDataOps_levels[level]));
///////
///////				data_container.back()->setup(sphereDataConfig.back());
///////
///////	#endif
///////
///////				std::cout << "Spectral resolution at level " << level << " : " << N_spectral[0] << " " << N_spectral[1] << std::endl;
///////			}
///////			else
///////			{
///////	#if SWEET_XBRAID_CART2D
///////				cart2DDataConfig.push_back(base_cart2DDataConfig);
///////				op_cart2d.push_back(base_op_cart2d);
///////	#else
///////				sphere2DDataConfig.push_back(base_sphere2DDataConfig);
///////				op_sphere2D.push_back(base_op_sphere2D);
///////				op_sphere2D_complex.push_back(base_op_sphere2D_complex);
///////	#endif
///////			}
///////		}
///////#endif


		// get min_spectral size
		for (std::size_t i = 0; i < config_levels.size(); i++)
			min_spectral_size = std::min(	min_spectral_size,
								std::min((int)config_levels[i]->spectral_data_size[0], (int)config_levels[i]->spectral_data_size[1])
							);
/////#if SWEET_XBRAID_CART2D
/////		// get min_spectral size
/////		for (std::size_t i = 0; i < cart2DDataConfig.size(); i++)
/////			min_spectral_size = std::min(	min_spectral_size,
/////								std::min((int)cart2DDataConfig[i]->spectral_data_size[0], (int)cart2DDataConfig[i]->spectral_data_size[1])
/////							);
/////
/////#elif SWEET_XBRAID_SPHERE2D
/////		// get min_spectral size
/////		for (std::size_t i = 0; i < sphere2DDataConfig.size(); i++)
/////			min_spectral_size = std::min(	min_spectral_size,
/////								std::min(sphere2DDataConfig[i]->spectral_modes_m_max, sphere2DDataConfig[i]->spectral_modes_n_max)
/////							);
/////#endif

		// get buffer size
		////! To be updated depending on the tsm
		// Overestimated
		size_buffer = 0;
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
			size_buffer += N_vec * config_levels[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);

////#if SWEET_XBRAID_SCALAR
////		size_buffer = N_vec * sizeof(double);
////#elif SWEET_XBRAID_CART2D
////		////! To be updated depending on the tsm
////		// Overestimated
////		size_buffer = 0;
////		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
////			size_buffer += N_vec * cart2DDataConfig[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
////#elif SWEET_XBRAID_SPHERE2D
////		////! To be updated depending on the tsm
////		// Overestimated
////		size_buffer = 0;
////		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
////			size_buffer += N_vec * sphere2DDataConfig[level]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
////#endif


	}


	/*
	 * Get specific parameters for each of the N levels
	 * Input string must contain 1, 2 or N orders separated by comma:
	 *  - 1 tso: same tso for all levels;
	 *  - 2 tso: first one is for level 0 (finest); all coarse levels use the second one
	 *  - N tso: i-th level uses the i-th tso.
	 */
	template<typename T>
	std::vector<T> getLevelParameterFromParameters(std::string i_param_name, int i_order = 0)
	{

		if ( ! (
			i_param_name == "timestepping_method" ||
			i_param_name == "timestepping_order" ||
			i_param_name == "viscosity_order" ||
			i_param_name == "viscosity_coefficient"
			))
			SWEETErrorFatal("Wrong param_name " + i_param_name);

		if (i_param_name == "timestepping_order")
			SWEET_ASSERT(i_order == 1 || i_order == 2);

		std::vector<T> out = {};
		std::stringstream all_param;

		if (i_param_name == "timestepping_method")
		{
			all_param = std::stringstream(this->shackXBraid->xbraid_timestepping_method);
			if (this->shackXBraid->xbraid_timestepping_method.find('(') != std::string::npos)
				this->useNewTimeSteppers = true;
			else
				this->useNewTimeSteppers = false;
		}
		else if (i_param_name == "timestepping_order")
		{
			if (i_order == 1)
				all_param = std::stringstream(shackXBraid->xbraid_timestepping_order);
			else if (i_order == 2)
				all_param = std::stringstream(shackXBraid->xbraid_timestepping_order2);
		}
		else if (i_param_name == "viscosity_order")
			all_param = std::stringstream(shackXBraid->xbraid_viscosity_order);
		else if (i_param_name == "viscosity_coefficient")
			all_param = std::stringstream(shackXBraid->xbraid_viscosity_coefficient);

		if ( ( i_param_name != "timestepping_method" ) || !this->useNewTimeSteppers  )
		{
			while (all_param.good())
			{
				std::string str;
				getline(all_param, str, ',');
				std::stringstream ss(str);
				T conv;
				if (ss >> conv)
					out.push_back(conv);
				else
					SWEETErrorFatal("Unable to convert parameter: " + str);
			}
		}
		else
		{
			// string should be: "(...,...,...);(...,...,...);..."
			while (all_param.good())
			{
				std::string str;
				getline(all_param, str, ';');
				std::stringstream ss(str);
				T conv;
				if (ss >> conv)
					out.push_back(conv);
				else
					SWEETErrorFatal("Unable to convert parameter: " + str);
			}
		}

		if ( ! (out.size() == 1 || out.size() == 2 || (int)out.size() == shackXBraid->xbraid_max_levels ) )
			SWEETErrorFatal("xbraid_" + i_param_name +  "must contain 1, 2 or N timestepping orders.");

		// all levels use same param
		if (out.size() == 1)
			for (int level = 1; level < shackXBraid->xbraid_max_levels; level++)
				out.push_back(out[0]);

		// all coarse levels use same tso
		if (out.size() == 2)
			for (int level = 2; level < shackXBraid->xbraid_max_levels; level++)
				out.push_back(out[1]);

		return out;
	}



public:
	sweet_BraidVector* create_new_vector(int i_level)
	{
		sweet_BraidVector* U = new sweet_BraidVector(config_levels[i_level], i_level);
/////#if SWEET_XBRAID_SCALAR
/////		sweet_BraidVector* U = new sweet_BraidVector(i_level);
/////#elif SWEET_XBRAID_CART2D
/////		sweet_BraidVector* U = new sweet_BraidVector(cart2DDataConfig[i_level], i_level);
/////#elif SWEET_XBRAID_SPHERE2D
/////		sweet_BraidVector* U = new sweet_BraidVector(sphere2DDataConfig[i_level], i_level);
/////#endif
		return U;
	}

private:
	void store_prev_solution(
					sweet_BraidVector* i_U,
					int i_time_id,
					int i_level,
					int iter
				)
	{
		// if not SL scheme: nothing to do
		//if ( std::find(SL_tsm.begin(), SL_tsm.end(), tsms[i_level]) == SL_tsm.end())
		if ( ! is_SL[i_level] )
			return;

		// if solution has already been stored in this iteration: nothing to do
		if ( sol_prev_iter[i_level][i_time_id] == iter )
			return;
		///SWEET_ASSERT(sol_prev_iter[i_level][i_time_id] == iter - 1);

		// create vector if necessary
		if ( ! sol_prev[i_level][i_time_id] )
			sol_prev[i_level][i_time_id] = create_new_vector(i_level);

		// set solution
		*sol_prev[i_level][i_time_id] = *i_U;
		sol_prev_iter[i_level][i_time_id] = iter;
		first_timeid_level[i_level] = std::min(first_timeid_level[i_level], i_time_id);
		last_timeid_level[i_level] = std::max(last_timeid_level[i_level], i_time_id);
	}

	void set_prev_solution(
					sweet_BraidVector* i_U,
					int i_time_id,
					int i_level
				)
	{

		// if not SL scheme: nothing to do
		///if (  std::find(SL_tsm.begin(), SL_tsm.end(), tsms[i_level]) == SL_tsm.end())
		if ( ! is_SL[i_level] )
			return;

		// if t == 0 or prev solution not available
		// then: prev_solution = solution
		bool prev_sol_exists = true;

		if ( i_time_id == 0 )
			prev_sol_exists = false;
		if ( prev_sol_exists && (!sol_prev[i_level][i_time_id - 1]) )
			prev_sol_exists = false;

		// only store prev solution if it is not the first time step inside a coarse slice
		//if (i_time_id % shackXBraid->xbraid_cfactor == 0)
		if (i_level < nlevels - 1)
			prev_sol_exists = false;

		if (prev_sol_exists)
			timeSteppers[i_level]->timestepper->set_previous_solution(sol_prev[i_level][i_time_id - 1]->data);
		else
			timeSteppers[i_level]->timestepper->set_previous_solution(i_U->data);
	}

public:
	/* --------------------------------------------------------------------
	 * Time integrator routine that performs the update
	 *   u_i = Phi_i(u_{i-1}) + g_i 
	 * 
	 * When Phi is called, u is u_{i-1}.
	 * The return value is that u is set to u_i upon completion
	 *
	 * Always receives and returns a solution defined on the finest spatial grid
	 * Spatial interpolation is performed if necessary
	 * -------------------------------------------------------------------- */
	braid_Int
	Step(
			braid_Vector		io_U,
			braid_Vector		i_ustop,
			braid_Vector		i_fstop,
			BraidStepStatus&	io_status
			)
	{

		// First call of this function: create timesteppers
		// Benchmark::setup_initial_conditions has already been called
		if (timeSteppers.size() == 0)
		{

			// if rank > 0, call setup_init_conditions to ensure that parameters from benchmark are set
			if (rank > 0)
			{
				braid_Vector dummy;
				Init(0., &dummy);
			}

			setup_timesteppers();

			io_status.GetNLevels(&nlevels);
		}


		// Vector defined in the finest level
		sweet_BraidVector* U = (sweet_BraidVector*) io_U;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;
		int time_id;
		int iter;

		/* Grab status of current time step */
		io_status.GetTstartTstop(&tstart, &tstop);
		io_status.GetLevel(&level);
		io_status.GetTIndex(&time_id);
		io_status.GetIter(&iter);

		// Vector defined in the current level (defined via interpolation if necessary)
		sweet_BraidVector* U_level = create_new_vector(level);

		// Interpolate to coarser grid in space if necessary
		if (shackXBraid->xbraid_spatial_coarsening && level > 0)
		/////if (shackXBraid->xbraid_spatial_coarsening)
			U_level->data->restrict(*U->data);
		else
			*U_level->data = *U->data;


		// create containers for prev solution
		if (sol_prev[level].size() == 0)
		{
			// store nt solutions (overestimated for coarse levels!)
			int nt;
			io_status.GetNTPoints(&nt);
			for (int i = 0; i < nt + 1; i++)
			{
				sol_prev[level].push_back(nullptr);
				sol_prev_iter[level].push_back(-1);
			}
				///sol_prev[level].push_back(create_new_vector());
		}

		// store solution for SL
		if (time_id == 0)
			store_prev_solution(U_level, time_id, level, iter);

		// set prev solution for SL
		set_prev_solution(U_level, time_id, level);

		// TODO: check if this is thread safe
		/////simVars->timecontrol.current_simulation_time = tstart;
		/////simVars->timecontrol.current_timestepSize = tstop - tstart;
		// TODO

		//////std::cout << rank << " " << iter << " " << level << " " << tstart << " " << tstop << std::endl;
		if (this->useNewTimeSteppers)
		{
			//this->timeSteppers_NewTS[level]->timeIntegrator->_eval_integration(
			this->timeSteppers_NewTS[level]->runIntegration(
					///dataConfigOps.prog,
					///dataConfigOps.progTmp,
					*U_level->data_DE,
					*U_level->data_DE_tmp,
					tstart
				);
			///dataConfigOps.prog.swap(dataConfigOps.progTmp);
			U_level->data_DE->swap(*U_level->data_DE_tmp);
		}
		else
			this->timeSteppers[level]->timestepper->runTimestep(
									U_level->data,
									tstop - tstart,
									tstart
			);

		// TODO: use DE term for viscosity!!!
		// Apply viscosity a posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral
		///if (simVars->sim.viscosity != 0 && simVars->misc.use_nonlinear_only_visc == 0)
#if SWEET_XBRAID_CART2D || SWEET_XBRAID_SPHERE2D
		if (viscosity_coefficients[level] != 0) /* && simVars->misc.use_nonlinear_only_visc == 0) */
		{
	#if SWEET_XBRAID_CART2D
			for (int i = 0; i < N_vec; i++)
			{
				///sweet::PlaneData_Spectral* field = U_level->data->get_pointer_to_data_PlaneData_Spectral()->simfields[i];
				sweet::Data::Cart2D::PlaneData_Spectral* field;
				if (this->useNewTimeStepper)
					field = &U_level->data_DE->data[i];
				else
					field = U_level->data->get_pointer_to_data_PlaneData_Spectral()->simfields[i];
				*field = this->op_plane[level]->implicit_diffusion(	*field,
											(tstop - tstart) * this->viscosity_coefficients[level],
											this->viscosity_orders[level]);
			}
	#elif SWEET_XBRAID_SPHERE2D
			for (int i = 0; i < N_vec; i++)
			{
				///sweet::SphereData_Spectral* field = U_level->data->get_pointer_to_data_SphereData_Spectral()->simfields[i];
				sweet::Data::Sphere2D::SphereData_Spectral* field;
				if (this->useNewTimeSteppers)
					field = &U_level->data_DE->data[i];
				else
					field = U_level->data->get_pointer_to_data_SphereData_Spectral()->simfields[i];
				*field = this->op_sphere[level]->implicit_hyperdiffusion(	*field,
												(tstop - tstart) * this->viscosity_coefficients[level],
												this->viscosity_orders[level],
												this->shackSphereDataOps->sphere_radius);
			}
	#endif
		}
#endif

		store_prev_solution(U_level, time_id + 1, level, iter);

		// Interpolate to finest grid in space if necessary
		if (this->useNewTimeSteppers)
		{
			if (this->shackXBraid->xbraid_spatial_coarsening && level > 0)
				U->data_DE->pad_zeros(*U_level->data_DE);
			else
				U->data_DE->op_setVector(*U_level->data_DE);
		}
		else
		{
			if (this->shackXBraid->xbraid_spatial_coarsening && level > 0)
				U->data->pad_zeros(*U_level->data);
			else
				*U->data = *U_level->data;
		}

		/* Tell XBraid no refinement */
		io_status.SetRFactor(1);

		delete U_level;

		return 0;
	}

		/* --------------------------------------------------------------------
		 * -------------------------------------------------------------------- */

	virtual braid_Int
	Residual(
				braid_Vector			i_ustop,
				braid_Vector			o_r,
				BraidStepStatus&		io_status
		)
	{

		//braid_StepStatus& status = (braid_StepStatus&) io_status;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;

		/* Grab status of current time step */
		io_status.GetTstartTstop(&tstart, &tstop);
	
		/* Grab level */
		io_status.GetLevel(&level);
	
		/* Set the new dt in the user's manager*/
		dt = tstop - tstart;


		sweet_BraidVector* u = (sweet_BraidVector*) i_ustop;
		sweet_BraidVector* r = (sweet_BraidVector*) o_r;

		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Create a vector object for a given time point.
	 * This function is only called on the finest level.
	 * -------------------------------------------------------------------- */
	braid_Int
	Init(
			double		i_t,
			braid_Vector*	o_U
			)
	{

		sweet_BraidVector* U = create_new_vector(0);

	// Set correct resolution in shacks
		for (int level = 0; level < shackXBraid->xbraid_max_levels; level++)
		{
			for (int j = 0; j < 2; j++)
			{
				shacksDataOps_levels[level]->space_res_physical[j] = config_levels[level]->grid_res[j];
				shacksDataOps_levels[level]->space_res_spectral[j] = config_levels[level]->spectral_data_size[j];
			}
/////#if SWEET_XBRAID_CART2D
/////			for (int j = 0; j < 2; j++)
/////			{
/////				shacksCart2DDataOps_levels[level]->space_res_physical[j] = cart2DDataConfig[level]->grid_res[j];
/////				shacksCart2DDataOps_levels[level]->space_res_spectral[j] = cart2DDataConfig[level]->spectral_data_size[j];
/////			}
/////#elif SWEET_XBRAID_SPHERE2D
/////			shacksSphere2DDataOps_levels[level]->space_res_physical[0] = sphere2DDataConfig[level]->grid_num_lon;
/////			shacksSphere2DDataOps_levels[level]->space_res_physical[1] = sphere2DDataConfig[level]->grid_num_lat;
/////			shacksSphere2DDataOps_levels[level]->space_res_spectral[0] = sphere2DDataConfig[level]->spectral_modes_m_max;
/////			shacksSphere2DDataOps_levels[level]->space_res_spectral[1] = sphere2DDataConfig[level]->spectral_modes_n_max;
/////#endif
		}






		if (i_t == tstart)
			*U->data_DE = prog_initial_solution;
		else
			*U->data_DE = prog_initial_guess;

		*o_U = (braid_Vector) U;













/////////	#if SWEET_XBRAID_SCALAR
/////////		double u0;
/////////	#elif SWEET_XBRAID_CART2D
/////////		#if SWEET_XBRAID_CART2D_SWE
/////////		sweet::Data::Cart2D::DataSpectral t0_prog_h_pert(cart2DDataConfig[0]);
/////////		#endif
/////////		sweet::Data::Cart2D::DataSpectral t0_prog_u(cart2DDataConfig[0]);
/////////		sweet::Data::Cart2D::DataSpectral t0_prog_v(cart2DDataConfig[0]);
/////////	#elif SWEET_XBRAID_SPHERE2D
/////////		sweet::Data::Sphere2D::DataSpectral t0_prog_phi_pert(sphere2DDataConfig[0]);
/////////		sweet::Data::Sphere2D::DataSpectral t0_prog_vrt(sphere2DDataConfig[0]);
/////////		sweet::Data::Sphere2D::DataSpectral t0_prog_div(sphere2DDataConfig[0]);
/////////	#endif
/////////
/////////
/////////		if ( i_t == tstart )
/////////		{
/////////			////std::cout << "Setting initial solution " << std::endl;
/////////	#if SWEET_XBRAID_SCALAR
/////////			u0 = shackBenchmarks->u0;
/////////			////////U->data->dataArrays_2_GenericData_Scalar(u0);
/////////	
/////////	#elif SWEET_XBRAID_CART2D
/////////		#if SWEET_XBRAID_CART2D_SWE
/////////			PDE_SWECart2D::Benchmarks::BenchmarksCombined sweCart2DBenchmarks;
/////////			sweCart2DBenchmarks.shackRegistration(shacksDict_levels[0]);
/////////			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sweCart2DBenchmarks);
/////////			sweCart2DBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, op_cart2d[0], cart2DDataConfig[0]);
/////////
/////////			// Dummy initialization in coarse levels
/////////			// The only purpose is to call Operators setup from benchmark (sim parameters may change!)
/////////			for (size_t level = 1; level < op_cart2d.size(); level++)
/////////			{
/////////				sweet::Data::Cart2D::DataSpectral dummy1(cart2DDataConfig[level]);
/////////				sweet::Data::Cart2D::DataSpectral dummy2(cart2DDataConfig[level]);
/////////				sweet::Data::Cart2D::DataSpectral dummy3(cart2DDataConfig[level]);
/////////
/////////				PDE_SWECart2D::Benchmarks::BenchmarksCombined sweCart2DBenchmarks_dummy;
/////////				sweCart2DBenchmarks_dummy.shackRegistration(shacksDict_levels[level]);
/////////				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sweCart2DBenchmarks_dummy);
/////////				sweCart2DBenchmarks_dummy.setupInitialConditions(dummy1, dummy2, dummy3, op_cart2d[level], cart2DDataConfig[level]);
/////////			}
/////////
/////////		#elif SWEET_XBRAID_CART2D_BURGERS
/////////			sweet::Data::Cart2D::DataGrid t0_prog_u_phys(t0_prog_u.cart2DDataConfig[0]);
/////////			sweet::Data::Cart2D::DataGrid t0_prog_v_phys(t0_prog_v.cart2DDataConfig[0]);
/////////			if (simVars->disc.space_grid_use_c_staggering)
/////////			{
/////////				t0_prog_u_phys.grid_update_lambda_array_indices(
/////////							[&](int i, int j, double &io_data)
/////////					{
/////////						double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.cart2d_domain_size[0];
/////////						double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.cart2d_domain_size[1];
/////////						io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
/////////					}
/////////				);
/////////				t0_prog_v_phys.grid_update_lambda_array_indices(
/////////							[&](int i, int j, double &io_data)
/////////					{
/////////					io_data = 0.0;
/////////					}
/////////				);
/////////			}
/////////			else
/////////			{
/////////				t0_prog_u_phys.grid_update_lambda_array_indices(
/////////							[&](int i, int j, double &io_data)
/////////					{
/////////						double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.cart2d_domain_size[0];
/////////						double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.cart2d_domain_size[1];
/////////						io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
/////////					}
/////////				);
/////////
/////////				t0_prog_v_phys.grid_update_lambda_array_indices(
/////////							[&](int i, int j, double &io_data)
/////////					{
/////////					io_data = 0.0;
/////////					}
/////////				);
/////////			}
/////////			t0_prog_u.loadCart2DDataGrid(t0_prog_u_phys);
/////////			t0_prog_v.loadCart2DDataGrid(t0_prog_v_phys);
/////////		#endif
/////////
/////////	#elif SWEET_XBRAID_SPHERE2D
/////////			PDE_SWESphere2D::Benchmarks::BenchmarksCombined sphere2DBenchmarks;
/////////
/////////			sphere2DBenchmarks.setup_1_registerAllBenchmark();
/////////			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);
/////////
/////////			sphere2DBenchmarks.setup_2_shackRegistration(shacksDict_levels[0]);
/////////			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);
/////////			////sphere2DBenchmarks.setup(*simVars_levels[0], *op_sphere2D[0]);
/////////
/////////			sphere2DBenchmarks.setup_3_benchmarkDetection();
/////////			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks);
/////////
/////////			sphere2DBenchmarks.setup_4_benchmarkSetup_1_withoutOps();
/////////			sphere2DBenchmarks.setup_5_benchmarkSetup_2_withOps(op_sphere2D[0]);
/////////
/////////			sphere2DBenchmarks.benchmark->getInitialState(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
/////////
/////////			// Dummy initialization in coarse levels
/////////			// The only purpose is to call Operators setup from benchmark (sim parameters may change!)
/////////			for (size_t level = 1; level < op_sphere2D.size(); level++)
/////////			{
/////////				sweet::Data::Sphere2D::DataSpectral dummy1(sphere2DDataConfig[level]);
/////////				sweet::Data::Sphere2D::DataSpectral dummy2(sphere2DDataConfig[level]);
/////////				sweet::Data::Sphere2D::DataSpectral dummy3(sphere2DDataConfig[level]);
/////////
/////////				PDE_SWESphere2D::Benchmarks::BenchmarksCombined sphere2DBenchmarks_dummy;
/////////
/////////				sphere2DBenchmarks_dummy.setup_1_registerAllBenchmark();
/////////				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks_dummy);
/////////
/////////				sphere2DBenchmarks_dummy.setup_2_shackRegistration(shacksDict_levels[level]);
/////////				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks_dummy);
/////////				///sphere2DBenchmarks_dummy.setup(*simVars_levels[level], *op_sphere2D[level]);
/////////
/////////				sphere2DBenchmarks_dummy.setup_3_benchmarkDetection();
/////////				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DBenchmarks_dummy);
/////////
/////////				sphere2DBenchmarks_dummy.setup_4_benchmarkSetup_1_withoutOps();
/////////				sphere2DBenchmarks_dummy.setup_5_benchmarkSetup_2_withOps(op_sphere2D[level]);
/////////
/////////				sphere2DBenchmarks_dummy.benchmark->getInitialState(dummy1, dummy2, dummy3);
/////////			}
/////////	#endif
/////////
/////////
/////////		}
/////////		else if (shackXBraid->xbraid_use_rand)
/////////		{
/////////			///std::cout << "Setting random solution " << std::endl;
/////////	#if SWEET_XBRAID_SCALAR
/////////			u0 = ((double)braid_Rand())/braid_RAND_MAX;
/////////	#elif SWEET_XBRAID_CART2D
/////////
/////////		#if SWEET_XBRAID_CART2D_SWE
/////////			///Cart2DData_Spectral t0_prog_h_pert(cart2DDataConfig[0]);
/////////			sweet::Data::Cart2D::DataGrid t0_prog_h_phys(cart2DDataConfig[0]);
/////////		#endif
/////////
/////////			///Cart2DData_Spectral t0_prog_u(cart2DDataConfig[0]);
/////////			////Cart2DData_Spectral t0_prog_v(cart2DDataConfig[0]);
/////////
/////////			sweet::Data::Cart2D::DataGrid t0_prog_u_phys(cart2DDataConfig[0]);
/////////			sweet::Data::Cart2D::DataGrid t0_prog_v_phys(cart2DDataConfig[0]);
/////////
/////////		#if SWEET_XBRAID_CART2D_SWE
/////////			t0_prog_h_phys.grid_update_lambda_array_indices(
/////////						[&](int i, int j, double &io_data)
/////////				{
/////////					io_data = shacksModel_levels[0]->h0 + ((double)braid_Rand())/braid_RAND_MAX;
/////////				}
/////////			);
/////////		#endif
/////////			t0_prog_u_phys.grid_update_lambda_array_indices(
/////////						[&](int i, int j, double &io_data)
/////////				{
/////////					io_data = ((double)braid_Rand())/braid_RAND_MAX;
/////////				}
/////////			);
/////////			t0_prog_v_phys.grid_update_lambda_array_indices(
/////////						[&](int i, int j, double &io_data)
/////////				{
/////////					io_data = ((double)braid_Rand())/braid_RAND_MAX;
/////////				}
/////////			);
/////////
/////////		#if SWEET_XBRAID_CART2D_SWE
/////////			t0_prog_h_pert.loadCart2DDataGrid(t0_prog_h_phys);
/////////		#endif
/////////			t0_prog_u.loadCart2DDataGrid(t0_prog_u_phys);
/////////			t0_prog_v.loadCart2DDataGrid(t0_prog_v_phys);
/////////
/////////	#elif SWEET_XBRAID_SPHERE2D
/////////			sweet::Data::Sphere2D::DataGrid t0_prog_phi_pert_phys(sphere2DDataConfig[0]);
/////////			sweet::Data::Sphere2D::DataGrid t0_prog_vrt_phys(sphere2DDataConfig[0]);
/////////			sweet::Data::Sphere2D::DataGrid t0_prog_div_phys(sphere2DDataConfig[0]);
/////////
/////////			t0_prog_phi_pert_phys.grid_update_lambda_array(
/////////						[&](int i, int j, double &io_data)
/////////				{
/////////					io_data = shackPDESWESphere2D->h0 + ((double)braid_Rand())/braid_RAND_MAX;
/////////				}
/////////			);
/////////			t0_prog_vrt_phys.grid_update_lambda_array(
/////////						[&](int i, int j, double &io_data)
/////////				{
/////////					io_data = ((double)braid_Rand())/braid_RAND_MAX;
/////////				}
/////////			);
/////////			t0_prog_div_phys.grid_update_lambda_array(
/////////						[&](int i, int j, double &io_data)
/////////				{
/////////					io_data = ((double)braid_Rand())/braid_RAND_MAX;
/////////				}
/////////			);
/////////
/////////			t0_prog_phi_pert.loadSphere2DDataGrid(t0_prog_phi_pert_phys);
/////////			t0_prog_vrt.loadSphere2DDataGrid(t0_prog_vrt_phys);
/////////			t0_prog_div.loadSphere2DDataGrid(t0_prog_div_phys);
/////////	#endif
/////////		}
/////////		else
/////////		{
/////////			/* Sets U as an all zero vector*/
/////////	#if SWEET_XBRAID_SCALAR
/////////			u0 = 0;
/////////	#elif SWEET_XBRAID_CART2D
/////////		#if SWEET_XBRAID_CART2D_SWE
/////////			t0_prog_h_pert.spectral_setZero();
/////////		#endif
/////////			t0_prog_u.spectral_setZero();
/////////			t0_prog_v.spectral_setZero();
/////////	#elif SWEET_XBRAID_SPHERE2D
/////////			t0_prog_phi_pert.spectral_setZero();
/////////			t0_prog_vrt.spectral_setZero();
/////////			t0_prog_div.spectral_setZero();
/////////	#endif
/////////		}
/////////
/////////
/////////		if (this->useNewTimeSteppers)
/////////		{
/////////	#if SWEET_XBRAID_SCALAR
/////////			U->data_DE->setData(u0);
/////////	#elif SWEET_XBRAID_CART2D
/////////			U->data_DE->setData(
/////////			#if SWEET_XBRAID_CART2D_SWE
/////////											t0_prog_h_pert,
/////////			#endif
/////////											t0_prog_u,
/////////											t0_prog_v);
/////////	#elif SWEET_XBRAID_SPHERE
/////////			U->data_DE->setData(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
/////////	#endif
/////////		}
/////////		else
/////////		{
/////////	#if SWEET_XBRAID_SCALAR
/////////			U->data->dataArrays_2_GenericData_Scalar(u0);
/////////	#elif SWEET_XBRAID_PLANE
/////////			U->data->dataArrays_2_GenericData_PlaneData_Spectral(
/////////			#if SWEET_XBRAID_PLANE_SWE
/////////											t0_prog_h_pert,
/////////			#endif
/////////											t0_prog_u,
/////////											t0_prog_v);
/////////	#elif SWEET_XBRAID_SPHERE
/////////			U->data->dataArrays_2_GenericData_SphereData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
/////////	#endif
/////////		}
/////////
/////////		*o_U = (braid_Vector) U;

		return 0;
	}



	/* --------------------------------------------------------------------
	 * Create a copy of a vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Clone(
			braid_Vector	i_U,
			braid_Vector*	o_V
		)
	{
		sweet_BraidVector* U = (sweet_BraidVector*) i_U;
		sweet_BraidVector* V = create_new_vector(U->level);
		*V = *U;
		*o_V = (braid_Vector) V;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Destroy vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Free(
			braid_Vector	i_U)
	{
		sweet_BraidVector* U = (sweet_BraidVector*) i_U;
		delete U;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Compute vector sum y = alpha*x + beta*y.
	 * -------------------------------------------------------------------- */
	braid_Int
	Sum(
			double			i_alpha,
			braid_Vector		i_X,
			double			i_beta,
			braid_Vector		io_Y
		)
	{

		sweet_BraidVector* X = (sweet_BraidVector*) i_X;
		sweet_BraidVector* Y = (sweet_BraidVector*) io_Y;

		*Y = *X * i_alpha + *Y * i_beta;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * User access routine to spatial solution vectors and allows for user
	 * output.  The default XBraid parameter of access_level=1, calls 
	 * my_Access only after convergence and at every time point.
	 * -------------------------------------------------------------------- */
	braid_Int
	Access(
				braid_Vector		i_U,
				BraidAccessStatus&	io_astatus
			)
	{
		double     tstart         = (tstart);
		double     tstop          = (tstop);
		////int        nt             = (nt);
	
		double     rnorm, disc_err, t;
		int        iter, level, done, index, myid, it;
		char       filename[255], filename_mesh[255], filename_err[255], filename_sol[255];


		sweet_BraidVector *U = (sweet_BraidVector*) i_U;

		/* Retrieve current time from Status Object */
		////braid_AccessStatusGetT(astatus, &t);
		io_astatus.GetT(&t);
		io_astatus.GetTIndex(&it);
		io_astatus.GetIter(&iter);
		io_astatus.GetLevel(&level);

		/* Retrieve XBraid State Information from Status Object */
		///////////MPI_Comm_rank(app->comm_x, &myid);
		///////////braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
		///////////braid_AccessStatusGetResidual(astatus, &rnorm);


		////std::cout << "ACCESS " << rank << " " << level << " " << t << std::endl;

		if (level == 0 /*&& rank == 0*/)
		{


			// Decide whether to output or not
			bool do_output = false;
			double small = 1e-10;

			// output each time step if:
			// output_timestep < 0 (i.e. output every timestep)
			// t == 0
			// t == Tmax
			// t is a multiple of dt_output
			if (
				shackIOData->outputEachSimTime < 0 ||
				std::abs(t) < small ||
				std::abs(t - shackTimestepControl->maxSimulationTime) < small ||
				fmod(t, shackIOData->outputEachSimTime) < small ||
				std::abs(fmod(t, shackIOData->outputEachSimTime) - shackIOData->outputEachSimTime) < small
			)
				do_output = true;

			if (shackXBraid->xbraid_no_output)
				do_output = false;

			////if (do_output)
			////	std::cout << "AAA " << t << " " << it << " " <<  t * simVars->iodata.output_time_scale << " " << simVars->iodata.output_each_sim_seconds << " " << fmod(t, simVars->iodata.output_each_sim_seconds) << " " << do_output << std::endl;

			if (do_output)
			{
				// Output physical solution to file
				if (shackXBraid->xbraid_store_iterations)
				{
					if (this->useNewTimeSteppers)
						this->output_data_file(
									U->data_DE,
									iter,
									it,
									t
						);
					else
						this->output_data_file(
									U->data,
									iter,
									it,
									t
						);
				}
				// Compute and store errors w.r.t. ref solution
				if (shackXBraid->xbraid_load_ref_csv_files)
				{
					// create containers for ref solution
					if (xbraid_data_ref_exact.size() == 0)
					{
						int nt;
						io_astatus.GetNTPoints(&nt);
						for (int i = 0; i < nt + 1; i++)
							xbraid_data_ref_exact.push_back(create_new_vector(0));
					}

					if (it >= 0)
					{
						if (this->useNewTimeSteppers)
							this->store_pint_error(
											U->data_DE,
											this->xbraid_data_ref_exact[it]->data_DE,
											N_vec,
											iter /* + 1 */,
											it,
											t,
											this->shackXBraid->xbraid_path_ref_csv_files,
											"ref",
											"xbraid"
							);
						else
							this->store_pint_error(
											U->data,
											this->xbraid_data_ref_exact[it]->data,
											N_vec,
											iter /* + 1 */,
											it,
											t,
											this->shackXBraid->xbraid_path_ref_csv_files,
											"ref",
											"xbraid"
							);
					}
				}
				// Compute and store errors w.r.t. fine (serial) solution
				if (shackXBraid->xbraid_load_fine_csv_files)
				{
					// create containers for fine solution
					if (xbraid_data_fine_exact.size() == 0)
					{
						int nt;
						io_astatus.GetNTPoints(&nt);
						for (int i = 0; i < nt + 1; i++)
							xbraid_data_fine_exact.push_back(create_new_vector(0));
					}

					if (it >= 0)
					{
						if (this->useNewTimeSteppers)
						{
							this->store_pint_error(
											U->data_DE,
											this->xbraid_data_fine_exact[it]->data_DE,
											N_vec,
											iter /* + 1 */,
											it,
											t,
											this->shackXBraid->xbraid_path_fine_csv_files,
											"fine",
											"xbraid"
							);
						}
						else
						{
							this->store_pint_error(
											U->data,
											this->xbraid_data_fine_exact[it]->data,
											N_vec,
											iter /* + 1 */,
											it,
											t,
											this->shackXBraid->xbraid_path_fine_csv_files,
											"fine",
											"xbraid"
							);
						}
					}
				}
			}

			// Store residual (residual per iteration)
			if (it == 0) {
				double res;
				io_astatus.GetResidual(&res);
				output_residual_file(res,
								iter);
			}



		}


		// TODO: verify if convergence stagnates and stop simulation

		return 0;
	}



	/* --------------------------------------------------------------------
	 * Compute norm of a spatial vector 
	 * -------------------------------------------------------------------- */
	braid_Int
	SpatialNorm(
			braid_Vector  i_U,
			double*       o_norm)
	{
		sweet_BraidVector* U = (sweet_BraidVector*) i_U;

		// Compute residual in the coarsest level
		// Restrict solution in spectral space then compute residual in physical space

		int max_level = (int)this->config.size() - 1;
		sweet_BraidVector* U_level = this->create_new_vector(max_level);
		U_level->data_DE->restrict(*U->data_DE);
		*o_norm = U_level->data_DE->reduceMaxAbs();
		delete U_level;


		///////////////if (this->useNewTimeSteppers)
		///////////////{
		///////////////#if SWEET_XBRAID_SCALAR
		///////////////		*o_norm = U->data_DE->grid_reduce_maxAbs();
		///////////////#elif SWEET_XBRAID_CART2D || SWEET_XBRAID_SPHERE2D
		///////////////	#if SWEET_XBRAID_CART2D
		///////////////		int max_level = (int)this->cart2DDataConfig.size() - 1;
		///////////////	#elif SWEET_XBRAID_SPHERE2D
		///////////////		int max_level = (int)this->sphere2DDataConfig.size() - 1;
		///////////////	#endif
		///////////////		sweet_BraidVector* U_level = this->create_new_vector(max_level);
		///////////////		U_level->data_DE->restrict(*U->data_DE);
		///////////////		*o_norm = U_level->data_DE->reduceMaxAbs();
		///////////////		delete U_level;
		///////////////#endif
		///////////////}
		///////////////else
		///////////////{
		///////////////#if SWEET_XBRAID_SCALAR
		///////////////		*o_norm = U->data->physical_reduce_maxAbs();
		///////////////#elif SWEET_XBRAID_CART2D || SWEET_XBRAID_SPHERE2D
		///////////////	#if SWEET_XBRAID_CART2D
		///////////////		int max_level = (int)this->cart2DDataConfig.size() - 1;
		///////////////	#elif SWEET_XBRAID_SPHERE2D
		///////////////		int max_level = (int)this->sphere2DDataConfig.size() - 1;
		///////////////	#endif
		///////////////		sweet_BraidVector* U_level = this->create_new_vector(max_level);
		///////////////		U_level->data->restrict(*U->data);
		///////////////		*o_norm = U_level->data->grid_reduce_maxAbs();
		///////////////		delete U_level;
		///////////////#endif
		///////////////}

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Return buffer size needed to pack one spatial braid_Vector.  Here the
	 * vector contains one double at every grid point and thus, the buffer 
	 * size is the number of grid points.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufSize(
			int*			o_size,
			BraidBufferStatus&	o_status)
	{
		*o_size = size_buffer;
		return 0;
	}


	/* --------------------------------------------------------------------
	 * Pack a braid_Vector into a buffer.
	 *
	 * Issue concerning SL: the first time step in each processor needs to receive
	 *                      the penult time step from the previous processor;
	 *                      However, communication is made only in level 0
	 * Solution (possibly not optimal): if SL is used in at least one level,
	 *                                  each communication includes the previous
	 *                                  time step of all levels.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufPack(
			braid_Vector		i_U,
			void*			o_buffer,
			BraidBufferStatus&	o_status
		)
	{

		sweet_BraidVector* U = (sweet_BraidVector*) i_U;

#if SWEET_XBRAID_SCALAR
		double* dbuffer = (double*) o_buffer;
#else
		std::complex<double>* dbuffer = (std::complex<double>*) o_buffer;
#endif

		// get buffer size
#if SWEET_XBRAID_SCALAR
		int actual_size_buffer = N_vec * sizeof(double);
#elif SWEET_XBRAID_CART2D
		int actual_size_buffer = N_vec * cart2DDataConfig[0]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
#elif SWEET_XBRAID_SPHERE2D
		int actual_size_buffer = N_vec * sphere2DDataConfig[0]->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
#endif

///#if SWEET_XBRAID_SCALAR
#if 1
		if (this->useNewTimeSteppers)
			U->data_DE->serialize(dbuffer);
		else
			U->data->serialize(dbuffer);
// communication of prev solution (SL) is no longer needed (see function store_prev_solution)
#elif SWEET_XBRAID_CART2D || SWEET_XBRAID_SPHERE2D
		// no SL method is used: only communicate solution
		if ( ! contains_SL )
			U->data->serialize(dbuffer);
		// SL method is used: also communcate prev solution
		else
		{
			// store solution from level 0
			std::complex<double> *level_buffer_data = nullptr;
	#if SWEET_XBRAID_CART2D
			int s = N * cart2DDataConfig[0]->spectral_array_data_number_of_elements;
	#elif SWEET_XBRAID_SPHERE2D
			int s = N * sphere2DDataConfig[0]->spectral_array_data_number_of_elements;
	#endif
			int s2 = 0;
			level_buffer_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(s * sizeof(std::complex<double>));
			U->data->serialize(level_buffer_data);
			std::copy(&level_buffer_data[0], &level_buffer_data[s], &dbuffer[s2]);
			s2 += s;
			sweet::Memory::MemBlockAlloc::free(level_buffer_data, s * sizeof(std::complex<double>));

			// store prev solution from every level using SL
			for (size_t level = 0; level < is_SL.size(); ++level)
			{
				if (is_SL[level])
				{
	#if SWEET_XBRAID_CART2D
					s = N * cart2DDataConfig[level]->spectral_array_data_number_of_elements;
	#elif SWEET_XBRAID_SPHERE2D
					s = N * sphere2DDataConfig[level]->spectral_array_data_number_of_elements;
	#endif
					level_buffer_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(s * sizeof(std::complex<double>));
					int time_id = last_timeid_level[level];
					if (time_id < 0)
						continue;
					sol_prev[level][time_id]->data->serialize(level_buffer_data);
					std::copy(&level_buffer_data[0], &level_buffer_data[s], &dbuffer[s2]);
					s2 += s;
					actual_size_buffer += s * sizeof(std::complex<double>);
					sweet::Memory::MemBlockAlloc::free(level_buffer_data, s * sizeof(std::complex<double>));
				}
			}
		}
#endif

		o_status.SetSize( actual_size_buffer );
		return 0;
	}

	/* --------------------------------------------------------------------
	 * Unpack a buffer and place into a braid_Vector
	 * -------------------------------------------------------------------- */
	braid_Int
	BufUnpack(
			void*			i_buffer,
			braid_Vector*		o_U,
			BraidBufferStatus&	io_status
		)
	{

		int level = 0;
		///io_status.GetLevel(&level);

		sweet_BraidVector* U = create_new_vector(level);

#if SWEET_XBRAID_SCALAR
		double* dbuffer = (double*) i_buffer;
#else
		std::complex<double>* dbuffer = (std::complex<double>*) i_buffer;
#endif

///#if SWEET_XBRAID_SCALAR
#if 1
		if (this->useNewTimeSteppers)
			U->data_DE->deserialize(dbuffer);
		else
			U->data->deserialize(dbuffer);
// communication of prev solution (SL) is no longer needed (see function store_prev_solution)
#elif SWEET_XBRAID_CART2D || SWEET_XBRAID_SPHERE2D
		// no SL method is used: only communicate solution
		if ( ! contains_SL )
			U->data->deserialize(dbuffer);
		// SL method is used: also communcate prev solution
		else
		{
			// get solution for level 0
			std::complex<double> *level_buffer_data = nullptr;
	#if SWEET_XBRAID_CART2D
			int s = N * cart2DDataConfig[0]->spectral_array_data_number_of_elements;
	#elif SWEET_XBRAID_SPHERE2D
			int s = N * sphere2DDataConfig[0]->spectral_array_data_number_of_elements;
	#endif
			int s2 = 0;
			level_buffer_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(s * sizeof(std::complex<double>));
			std::copy(&dbuffer[0], &dbuffer[s], &level_buffer_data[0]);
			U->data->deserialize(level_buffer_data);
			s2 += s;
			sweet::Memory::MemBlockAlloc::free(level_buffer_data, s * sizeof(std::complex<double>));

			// get prev solution for every level using SL
			for (size_t level = 0; level < is_SL.size(); ++level)
			{
				if (is_SL[level])
				{
	#if SWEET_XBRAID_CART2D
					s = N * cart2DDataConfig[level]->spectral_array_data_number_of_elements;
	#elif SWEET_XBRAID_SPHERE2D
					s = N * sphere2DDataConfig[level]->spectral_array_data_number_of_elements;
	#endif
					level_buffer_data = sweet::Memory::MemBlockAlloc::alloc<std::complex<double>>(s * sizeof(std::complex<double>));
					std::copy(&dbuffer[s2], &dbuffer[s2 + s], &level_buffer_data[0]);
					int time_id = first_timeid_level[level];
					if (time_id == INT_MAX)
						continue;
					sol_prev[level][time_id - 1] = create_new_vector(level);
					sol_prev[level][time_id - 1]->data->deserialize(level_buffer_data);
					s2 += s;
					sweet::Memory::MemBlockAlloc::free(level_buffer_data, s * sizeof(std::complex<double>));
				}
			}
		}
#endif

		*o_U = (braid_Vector) U;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Define time grid
	 * -------------------------------------------------------------------- */
	static braid_Int
	sweet_TimeGrid(
				_braid_App_struct* i_app,
				braid_Real* i_ta,
				braid_Int* i_ilower,
				braid_Int* i_iupper
			)
	{

		sweet_BraidApp* app =  (sweet_BraidApp*) i_app;

		double tstart;
		int lower = *i_ilower;
		int upper = *i_iupper;

		/* Start from the global tstart to compute the local tstart */
		tstart = app->tstart;
		for (int i = 0; i < lower; i++)
			tstart += app->custom_time_steps[i];

		/* Assign time point values for local time point index values lower:upper */
		for (int i = lower; i <= upper; i++)
		{
			i_ta[i - lower] = tstart;
			tstart += app->custom_time_steps[i];
		}

		return 0;
	}

};

}}

#endif
