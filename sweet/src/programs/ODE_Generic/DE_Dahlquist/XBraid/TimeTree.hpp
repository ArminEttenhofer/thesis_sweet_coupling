#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_XBRAID_TIMETREE_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_XBRAID_TIMETREE_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/XBraid/TimeTree.hpp>
#include <sweet/XBraid/Vector.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/TimeTreeIR.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/DataContainer/Simulation.hpp>
#include <programs/ODE_Generic/XBraid/DataContainer.hpp>


namespace ODE_Generic {
namespace DE_Dahlquist {
namespace XBraid {

class TimeTree :
	public ODE_Generic::XBraid::TimeTree
{

public:

	TimeTree()
	{
		tsm = new ODE_Generic::DE_Dahlquist::TimeTree;
		config = new ODE_Generic::DE_Dahlquist::DataContainer::Config;
	}

public:
	~TimeTree()
	{
		clear();
	}

};

}}}

#endif
