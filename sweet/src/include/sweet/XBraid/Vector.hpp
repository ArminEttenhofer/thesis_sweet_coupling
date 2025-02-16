/*
 * Vector.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef INCLUDE_SWEET_XBRAID_VECTOR_HPP
#define INCLUDE_SWEET_XBRAID_VECTOR_HPP

////#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
////#include <sweet/_DEPRECATED_pint/PInT_Common.hpp>

#include <xbraid/braid.hpp>

#if SWEET_GUI
#include<sweet/GUI/VisSweet.hpp>
#endif

#include <algorithm>

#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/XBraid/Shack.hpp>

#include <sweet/Data/GenericContainer/ConfigBase.hpp>

///#include "GeometryDependentDefinitions.hpp"

namespace sweet {
namespace XBraid {

/* --------------------------------------------------------------------
 * XBraid vector 
 * Stores the state of the simulation for a given time step
 * Define BraidVector, can contain anything, and be named anything
 * --> Put all time-dependent information here
 * -------------------------------------------------------------------- */
class Vector
{
public:

	int level;
	std::size_t N; // number of variables stored in data
	std::vector<std::string> var_names; // name of variables stored in data

	Vector()
	{
	}

	virtual ~Vector()
	{
	}

	virtual
	void clear() = 0;

	virtual
	void op_setVector(
			const sweet::XBraid::Vector& i_U
	) = 0;

	virtual
	void op_addVector(
			const sweet::XBraid::Vector &i_U
	) = 0;

	virtual
	void op_subVector(
			const sweet::XBraid::Vector &i_U
	) = 0;

	virtual
	void op_mulScalar(
			const double i_value
	) = 0;


	virtual
	void allocate_data() = 0;


	virtual
	void restrict(Vector* i_data) = 0;

	virtual
	void pad_zeros(Vector* i_data) = 0;

	virtual
	void op_setVector(Vector* i_data) = 0;


	virtual
	double reduceSum() const = 0;

	virtual
	double reduceMaxAbs() const = 0;

	virtual
	double reduceSum(int i) const = 0;

	virtual
	double reduceMaxAbs(int i) const = 0;

	virtual
	double reduceMaxAbs(int i, int rnorm) const = 0;

	virtual
	double reduceNormL1Grid(bool normalized = false) const = 0;

	virtual
	double reduceNormL2Grid(bool normalized = false) const = 0;

	virtual
	double reduceNormLinfGrid() const = 0;

	virtual
	double reduceNormL1Grid(int i, bool normalized = false) const = 0;

	virtual
	double reduceNormL2Grid(int i, bool normalized = false) const = 0;

	virtual
	double reduceNormLinfGrid(int i) const = 0;


	virtual
	void serialize(std::complex<double> *i_data) = 0;

	virtual
	void deserialize(std::complex<double> *i_data) = 0;

	virtual
	void fileLoad(
			const char* i_filename_template,
			std::string i_path,
			std::string i_output_file_mode,
			double i_t,
			double i_timescale
	) = 0;

};

}}

#endif
