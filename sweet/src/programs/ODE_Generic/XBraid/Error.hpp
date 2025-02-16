#ifndef PROGRAMS_ODE_GENERIC_XBRAID_ERROR_HPP
#define PROGRAMS_ODE_GENERIC_XBRAID_ERROR_HPP

#include <sweet/XBraid/Vector.hpp>
#include <sweet/XBraid/Error.hpp>
#include <programs/ODE_Generic/XBraid/DataContainer.hpp>
#include <programs/ODE_Generic/XBraid/TimeTree.hpp>

namespace ODE_Generic {
namespace XBraid {

class Error :
	public sweet::XBraid::Error
{

	//////sweet::Data::Cart2D::Config* config = nullptr;

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
	}

public:
	void setup(
		/////sweet::Data::Cart2D::Config* i_config,
		sweet::IO::Shack *i_shackIOData
	)
	{
		sweet::XBraid::Error:: setup(i_shackIOData);
		///////config = i_config;
	}

public:

	std::vector<size_t> getRnorms() override
	{
		std::vector<size_t> rnorms = {};

		//////for (int ip = 0; ip <= 5; ip++)
		//////{
		//////	int rnorm = config->spectral_data_size[0] / std::pow(2, ip);
		//////	if (rnorm >= 8)
		//////		rnorms.push_back(rnorm);
		//////}

		return rnorms;
	}



};

}}

#endif
