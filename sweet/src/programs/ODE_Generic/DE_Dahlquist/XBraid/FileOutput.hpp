#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_XBRAID_FILEOUTPUT_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_XBRAID_FILEOUTPUT_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/FileOutput.hpp>
#include <programs/ODE_Generic/XBraid/DataContainer.hpp>
#include <programs/ODE_Generic/XBraid/TimeTree.hpp>

#include <programs/ODE_Generic/DE_Dahlquist/FileOutput.hpp>

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace XBraid {

class FileOutput :
	public ODE_Generic::XBraid::FileOutput
{

public:
	FileOutput()
	{
		file_output = new ODE_Generic::DE_Dahlquist::FileOutput;
	}

public:
	~FileOutput()
	{
	}

public:
	void fileSave(
			sweet::XBraid::Vector* i_vector
	) override
	{
		ODE_Generic::DE_Dahlquist::XBraid::DataContainer* U = (ODE_Generic::DE_Dahlquist::XBraid::DataContainer*) i_vector;
		for (int i = 0; i < U->data2->N; i++)
			file_output->fileSave(U->data2->data[i], U->data2->var_names[i]);
	}

};

}}}

#endif
