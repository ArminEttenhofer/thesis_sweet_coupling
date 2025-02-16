/*
 * Storage.hpp
 *
 *  Created on: 7 March 2023
 *      Author: Thibaut LUNET <thibaut.lunet@tuhh.de>
 */
#ifndef INCLUDE_SWEET_SDC_STORAGE_HPP
#define INCLUDE_SWEET_SDC_STORAGE_HPP

#include <sweet/Data/Sphere2D/DataSpectral.hpp>

namespace sweet {
namespace SDC {

/*
 * Class to store solution data at one node
 */
class SWE_VariableVector
{
public:
	sweet::Data::Sphere2D::DataSpectral phi;
	sweet::Data::Sphere2D::DataSpectral vrt;
	sweet::Data::Sphere2D::DataSpectral div;

public:

	// Copy constructor
	SWE_VariableVector(const SWE_VariableVector &i_value)	:
		phi(i_value.phi.sphere2DDataConfig),
		vrt(i_value.phi.sphere2DDataConfig),
		div(i_value.phi.sphere2DDataConfig)
	{
		phi = i_value.phi;
		vrt = i_value.vrt;
		div = i_value.div;
	}

	// Default constructor
	SWE_VariableVector()
	{
	}

	// Fill values for phi, vort and div
	SWE_VariableVector& operator=(const SWE_VariableVector& u) {
		phi = u.phi;
		vrt = u.vrt;
		div = u.div;
		return *this;
	}

	bool setup(const sweet::Data::Sphere2D::Config *sphere2d_data_config)
	{
		phi.setup(sphere2d_data_config);
		vrt.setup(sphere2d_data_config);
		div.setup(sphere2d_data_config);

		return true;
	}

	void swap(SWE_VariableVector &io_value)
	{
		phi.swap(io_value.phi);
		vrt.swap(io_value.vrt);
		div.swap(io_value.div);
	}
};


// Class to store all the solution data to each nodes and two iterations
class SDC_NodeStorage {
	std::vector<SWE_VariableVector> data;


public:
	SDC_NodeStorage()
	{
	}

public:
	void setup(
			const sweet::Data::Sphere2D::Config* sphere2DDataConfig,
			size_t num_nodes
	){
		data.resize(num_nodes);

		for (size_t i = 0; i < num_nodes; i++)
		{
			data[i].setup(sphere2DDataConfig);
		}
	}

	SWE_VariableVector& operator[](int i)
	{
		return data[i];
	}

	void swap(SDC_NodeStorage &i_value)
	{
		for (size_t i = 0; i < data.size(); i++)
			data[i].swap(i_value.data[i]);
	}
};

}}

#endif
