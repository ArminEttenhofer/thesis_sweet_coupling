#ifndef PROGRAMS_PDE_SWESPHERE2D_XBRAID_ERROR_HPP
#define PROGRAMS_PDE_SWESPHERE2D_XBRAID_ERROR_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/Error.hpp>
#include <programs/PDE_SWESphere2D/XBraid/DataContainer.hpp>
#include <programs/PDE_SWESphere2D/XBraid/TimeTree.hpp>

namespace PDE_SWESphere2D {
namespace XBraid {

class Error :
	public sweet::XBraid::Error
{

	sweet::Data::Sphere2D::Config* config = nullptr;

public:
	Error()
	{
	}

public:
	~Error()
	{
		clear();
	}

public:
	void clear()
	{
		config = nullptr;
	}

public:
	void setup(
		sweet::Data::Sphere2D::Config* i_config,
		sweet::IO::Shack *i_shackIOData
	)
	{
		sweet::XBraid::Error:: setup(i_shackIOData);
		config = i_config;
	}

public:

	std::vector<size_t> getRnorms() override
	{
		std::vector<size_t> rnorms = {};

		for (int ip = 0; ip <= 5; ip++)
		{
			int rnorm = config->spectral_modes_m_max / std::pow(2, ip);
			if (rnorm >= 8)
				rnorms.push_back(rnorm);
		}

		return rnorms;
	}



};

}}

#endif
