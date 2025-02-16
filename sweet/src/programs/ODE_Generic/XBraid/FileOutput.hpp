#ifndef PROGRAMS_ODE_GENERIC_XBRAID_FILEOUTPUT_HPP
#define PROGRAMS_ODE_GENERIC_XBRAID_FILEOUTPUT_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/FileOutput.hpp>
#include <programs/ODE_Generic/XBraid/DataContainer.hpp>
#include <programs/ODE_Generic/XBraid/TimeTree.hpp>

#include <programs/ODE_Generic/DE_Dahlquist/FileOutput.hpp>

namespace ODE_Generic {
namespace XBraid {

class FileOutput :
	public sweet::XBraid::FileOutput
{

protected:
	ODE_Generic::FileOutput::Base* file_output = nullptr;

	ODE_Generic::Shack* shackModel = nullptr;

public:
	FileOutput()
	{
	}

public:
	~FileOutput()
	{
		clear();
	}

public:
	void clear()
	{
		if (file_output)
		{
			delete file_output;
			file_output = nullptr;
		}

		shackModel = nullptr;
	}

public:
	bool setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl,
			ODE_Generic::Shack* i_shackModel
		)
	{

		sweet::XBraid::FileOutput::setup(i_shackIOData, i_shackTimestepControl);

		shackModel = i_shackModel;

		file_output->setup(shackIOData, shackTimestepControl);
		return true;
	}


public:
	virtual
	void fileSave(
			sweet::XBraid::Vector* i_vector
	) = 0;
	//{
	//	ODE_Generic::XBraid::DataContainer* U = (ODE_Generic::XBraid::DataContainer*) i_vector;
	//	for (int i = 0; i < U->data->N; i++)
	//		file_output->fileSave(U->data->data[i], U->data->var_names[i]);
	//}


};

}}

#endif
