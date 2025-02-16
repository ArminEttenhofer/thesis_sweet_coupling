/*
 * namespaces.hpp
 *
 *  Created on: Jun 29, 2023
 *      Author: martin
 */

#ifndef INCLUDE_SWEET_NAMESPACES_HPP
#define INCLUDE_SWEET_NAMESPACES_HPP


/*!
 * \brief SWEET! Shallow Water Equation Environment for Tests, Awesome!
 */
namespace sweet
{

	/*!
	 * \brief Namespace for different data containers representing data on either primitives (Cartesian grid / sphere) or a vector
	 */
	namespace Data {}

	/*!
	 * \brief Simple dictionary to store and load data including a Python implementation
	 */
	namespace Dict {}

	/*!
	 * \brief Various helpers for error handling including stack tracing with gdb.
	 */
	namespace Error {}

	/*!
	 * \brief Exponential integration
	 */
	namespace ExpIntegration {}

	/*!
	 * \brief Helper classes for GUI
	 */
	namespace GUI {}

	/*!
	 * \brief I/O helper routines
	 */
	namespace IO {}

	/*!
	 * \brief Math helper functions
	 */
	namespace LibMath {}

	/*!
	 * \brief Memory management (allocation) helpers
	 */
	namespace Memory {}

	/*!
	 * \brief Parallelization helpers
	 */
	namespace Parallelization {}

	/*!
	 * \brief SDC time integration scheme
	 */
	namespace SDC {}

	/*!
	 * \brief Semi-Lagrangian methods
	 */
	namespace SemiLagrangian {}

	/*!
	 * \brief Main Shack implementation
	 */
	namespace Shacks {}

	/*!
	 * \brief Tree-composable time integration methods
	 */
	namespace TimeTree {}

	/*!
	 * \brief Helper tools (everything which doesn't fit to other directories)
	 */
	namespace Tools {}

	/*!
	 * \brief XBraid parallel-in-time framework support
	 */
	namespace XBraid {}
}


#endif
