/*
 * GeometryDependentDefinitions.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

// All includes, definitions and preprocessor directives used by XBraid and depending on the geometry (Scalar / Cart2D / Sphere2D)
// should be in this file


#ifndef INCLUDE_SWEET_XBRAID_GEOMETRYDEPENDENTDEFINITIONS_HPP
#define INCLUDE_SWEET_XBRAID_GEOMETRYDEPENDENTDEFINITIONS_HPP


#if SWEET_XBRAID_SCALAR
	///////#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Scalar.hpp>
	#include <sweet/Data/Scalar/Shack.hpp>
	#include <sweet/Data/Scalar/Config.hpp>
	#include <sweet/Data/Scalar/Operators.hpp>
	#include <programs/ODE_Scalar/ODEScalarTimeSteppers.hpp>
	#include <programs/ODE_Scalar/ShackODEScalar.hpp>
	#include <programs/ODE_Scalar/time/ShackODEScalarTimeDiscretization.hpp>
	#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/TimeTreeIR.hpp>
	#include <programs/ODE_Generic/DE_Dahlquist/FileOutput.hpp>
	#include <programs/ODE_Generic/DE_Dahlquist/Shack.hpp>
	#include <programs/ODE_Generic/DE_Dahlquist/Benchmarks/Shack.hpp>

namespace sweet {
namespace XBraid {

	///typedef ODEScalarTimeSteppers t_tsmType;
	typedef ODEScalarTimeSteppers t_tsmType;
	typedef ODE_Generic::DE_Dahlquist::TimeTree t_tsmType_NewTS;
	typedef ODE_Generic::DE_Dahlquist::DataContainer::Simulation t_DataContainer;
	typedef ODE_Generic::DE_Dahlquist::Shack t_ShackModel;
	typedef ShackODEScalarTimeDiscretization t_ShackTimeDiscretization;
	typedef ODE_Generic::DE_Dahlquist::Benchmarks::Shack t_ShackBenchmarks;
	typedef sweet::Data::Scalar::Shack t_ShackDataOps;
	typedef sweet::Data::Scalar::Config	t_Config;
	typedef sweet::Data::Scalar::Operators	t_Operators;
	typedef sweet::Data::Scalar::Operators	t_OperatorsComplex;
	typedef ODE_Generic::DE_Dahlquist::FileOutput t_FileOutput;
	#define N_vec 1

	std::vector<std::string> SL_tsm = {};
}}

#elif SWEET_XBRAID_CART2D
	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Cart2DData_Spectral.hpp>
	#include <sweet/Data/Cart2D/Shack.hpp>
	#if SWEET_XBRAID_CART2D_BURGERS
		#include <programs/pde_burgersCart2D/PDEBurgersCart2D_TimeSteppers.hpp>
		#include <programs/pde_burgersCart2D/ShackPDEBurgersCart2D.hpp>
		#include <programs/pde_burgersCart2D/time/ShackPDEBurgersCart2DTimeDiscretization.hpp>

namespace sweet {
namespace XBraid {

		typedef BurgersCart2DTimeSteppers t_tsmType;
		typedef BurgersCart2DTimeSteppersNewTS t_tsmType_NewTS;
		typedef ShackPDEBurgersCart2D t_ShackModel;
		////typedef ShackPDEBurgersCart2DTimeDiscretization t_ShackTimeDiscretization;
		typedef ShackPDEBurgersCart2DBenchmarks t_ShackBenchmarks;
}}

		#define N_vec 2
	#elif SWEET_XBRAID_CART2D_SWE
		#include <programs/PDE_SWECart2D/TimeSteppers.hpp>
		#include <programs/PDE_SWECart2D/BenchmarksCombined.hpp>
		#include <programs/PDE_SWECart2D/Shack.hpp>
		#include <programs/PDE_SWECart2D/TimeOld/Shack.hpp>
		#include <programs/PDE_SWECart2D/TimeTree/TimeTreeIR.hpp>
		#include <programs/PDE_SWECart2D/FileOutput.hpp>

namespace sweet {
namespace XBraid {

		typedef PDE_SWECart2D::TimeSteppers t_tsmType;
		typedef PDE_SWECart2D::TimeTree::TimeTree t_tsmType_NewTS;
		typedef PDE_SWECart2D::Shack t_ShackModel;
		typedef PDE_SWECart2D::DataContainer::Simulation t_DataContainer;
		typedef PDE_SWECart2D::TimeDiscretization::Shack t_ShackTimeDiscretization;
		typedef PDE_SWECart2D::Benchmarks::Shack t_ShackBenchmarks;
		typedef sweet::Data::Cart2D::Shack t_ShackDataOps;
		typedef sweet::Data::Cart2D::Config	t_Config;
		typedef sweet::Data::Cart2D::Operators	t_Operators;
		typedef sweet::Data::Cart2DComplex::Operators	t_OperatorsComplex;
		typedef PDE_SWECart2D::FileOutput t_FileOutput;

		std::vector<std::string> SL_tsm = { "l_cn_na_sl_nd_settls",
						    "l_rexi_na_sl_nd_etdrk",
						    "l_rexi_na_sl_nd_settls"
				};


}}

		#define N_vec 3
	#endif
#elif SWEET_XBRAID_SPHERE2D
	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Sphere2DData_Spectral.hpp>
	#include <sweet/Data/Sphere2D/Shack.hpp>
	#include <programs/PDE_SWESphere2D/TimeOld/PDESWESphere2D_TimeSteppers.hpp>
	#include <programs/PDE_SWESphere2D/TimeTree/TimeTree.hpp>
	#include <programs/PDE_SWESphere2D/BenchmarksCombined.hpp>
	#include <programs/PDE_SWESphere2D/Shack.hpp>
	#include <programs/PDE_SWESphere2D/TimeOld/ShackTimeDiscretization.hpp>
	#include <programs/PDE_SWESphere2D/DataContainer.hpp>
	#include <programs/PDE_SWESphere2D/FileOutput.hpp>

namespace sweet {
namespace XBraid {
	typedef PDESWESphere_TimeSteppers t_tsmType;
	typedef PDE_SWESphere::TimeTree::TimeTree t_tsmType_NewTS;
	typedef PDE_SWESphere::Shack t_ShackModel;
	typedef ShackTimeDiscretization t_ShackTimeDiscretization;
	typedef PDE_SWESphere::Benchmarks::Shack t_ShackBenchmarks;
	typedef PDE_SWESphere::DataContainer t_DataContainer;
	typedef sweet::Data::Sphere2D::Shack t_ShackDataOps;
	typedef sweet::Data::Sphere2D::Config	t_Config;
	typedef sweet::Data::Sphere2D::Operators	t_Operators;
	typedef sweet::Data::Sphere2DComplex::Operators	t_OperatorsComplex;
	typedef PDE_SWESphere2D::FileOutput t_FileOutput;

	std::vector<std::string> SL_tsm = { "lg_exp_na_sl_lc_nr_etd_uv",
					    "l_irk_na_sl_nr_settls_uv_only",
					    "l_irk_na_sl_nr_settls_vd_only",
					    "l_irk_na_sl_settls_uv_only",
					    "l_irk_na_sl_settls_vd_only",
					    "ln_sl_exp_settls_uv",
					    "ln_sl_exp_settls_vd",
					    "lg_exp_na_sl_lc_nr_etdrk_uv"
					   };


}}
	#define N_vec 3
#endif

namespace sweet {
namespace XBraid {

#if SWEET_XBRAID_CART2D
	// Grid Mapping (staggered grid)
	sweet::Data::Cart2D::GridMapping gridMapping;
#endif

}}

#endif
