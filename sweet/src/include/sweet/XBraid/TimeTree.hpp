/*
 * TimeTree.hpp
 *
 *  Created on: 19 Jul 2023
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef INCLUDE_SWEET_XBRAID_TIMETREE_HPP
#define INCLUDE_SWEET_XBRAID_TIMETREE_HPP

#include <xbraid/braid.hpp>

#if SWEET_GUI
#include<sweet/GUI/VisSweet.hpp>
#endif

#include <algorithm>

#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/XBraid/Shack.hpp>
#include <sweet/XBraid/Vector.hpp>

namespace sweet {
namespace XBraid {

/* --------------------------------------------------------------------
 * XBraid time stepper wrapper
 * -------------------------------------------------------------------- */
class TimeTree
{
public:


	TimeTree()
	{
	}

	virtual ~TimeTree()
	{
	}

	TimeTree(const TimeTree &i_timetree)
	{
	};


	virtual
	bool runIntegration(
				const sweet::XBraid::Vector* i_U,
				const sweet::XBraid::Vector* o_U,
				double i_simulationTime
	) = 0;

	virtual
	bool storePreviousSolution(
				const sweet::XBraid::Vector* i_U
	) = 0;

};

}}

#endif
