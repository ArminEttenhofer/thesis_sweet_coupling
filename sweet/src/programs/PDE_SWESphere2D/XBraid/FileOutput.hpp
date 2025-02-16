#ifndef PROGRAMS_PDE_SWESPHERE2D_XBRAID_FILEOUTPUT_HPP
#define PROGRAMS_PDE_SWESPHERE2D_XBRAID_FILEOUTPUT_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/FileOutput.hpp>
#include <programs/PDE_SWESphere2D/XBraid/DataContainer.hpp>
#include <programs/PDE_SWESphere2D/XBraid/TimeTree.hpp>
#include <programs/PDE_SWESphere2D/FileOutput.hpp>

namespace PDE_SWESphere2D {
namespace XBraid {

class FileOutput :
	public sweet::XBraid::FileOutput
{

	PDE_SWESphere2D::FileOutput file_output;

	sweet::Data::Sphere2D::Shack* shackDataOps = nullptr;
	PDE_SWESphere2D::Shack* shackModel = nullptr;
	sweet::Data::Sphere2D::Config* config = nullptr;
	sweet::Data::Sphere2D::Operators* ops = nullptr;
	sweet::Data::Sphere2DComplex::Operators* opsComplex = nullptr;


public:
	FileOutput()
	{
	}

public:
	~FileOutput()
	{
	}

public:
	bool setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl,
			sweet::Data::Sphere2D::Shack* i_shackDataOps,
			PDE_SWESphere2D::Shack* i_shackModel,
			sweet::Data::Sphere2D::Config* i_config,
			sweet::Data::Sphere2D::Operators* i_ops,
			sweet::Data::Sphere2DComplex::Operators* i_opsComplex
		)
	{

		sweet::XBraid::FileOutput::setup(i_shackIOData, i_shackTimestepControl);

		config = i_config;
		shackDataOps = i_shackDataOps;
		shackModel = i_shackModel;
		ops = i_ops;
		opsComplex = i_opsComplex;

		file_output.setup(shackIOData, shackTimestepControl, shackModel);
		return true;
	}


public:
	void fileSave(
			sweet::XBraid::Vector* i_vector
	) override
	{
		PDE_SWESphere2D::XBraid::DataContainer* U = (PDE_SWESphere2D::XBraid::DataContainer*) i_vector;
		for (int i = 0; i < U->data->N; i++)
			file_output.fileSave(U->data->data[i], U->data->var_names[i]);
	}


};

}}

#endif
