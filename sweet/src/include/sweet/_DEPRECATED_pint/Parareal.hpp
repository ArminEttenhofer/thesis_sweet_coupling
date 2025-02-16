/*
 * Parareal.hpp
 *
 *  Created on: 11 Apr 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_HPP




#if SWEET_PARAREAL==0

namespace sweet {
namespace DEPRECATED_pint {


/*
 * Create empty Parareal implementations?
 */
class PararealSimulation_Base{};
class PararealData{};
template <class T> class PararealDataInherited	: public T {};

}}

#elif SWEET_PARAREAL==1

#	include <sweet/_DEPRECATED_pint/Parareal_SimulationInstance.hpp>
#	include <sweet/_DEPRECATED_pint/Parareal_Controller.hpp>
#	include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>


#elif SWEET_PARAREAL==2

	#if !SWEET_MPI
		#error "SWEET_MPI must be activated to use SWEET_PARAREAL=2"
	#endif

#	include <sweet/_DEPRECATED_pint/Parareal_SimulationInstance.hpp>
#	include <sweet/_DEPRECATED_pint/Parareal_Controller.hpp>
#	include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>

/////#	error "Parareal with MPI implemented but still requires a complete validation."

#endif



#endif
