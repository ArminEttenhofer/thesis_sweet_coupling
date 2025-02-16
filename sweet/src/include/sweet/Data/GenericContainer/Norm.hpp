/*
 * DESolver_DataContainer_Base.hpp
 */

#ifndef INCLUDE_SWEET_DATA_GENERICCONTAINER_NORM_HPP
#define INCLUDE_SWEET_DATA_GENERICCONTAINER_NORM_HPP

#if SWEET_MPI
#include <mpi.h>
#endif

#include <complex>

#include <sweet/Error/Fatal.hpp>

namespace sweet {
namespace Data {
namespace GenericContainer {

class Norm
{
public:
	Norm()
	{
	}

#if 1
	virtual
	Norm* getNewDataContainer() const = 0;
#endif

	virtual
	void clear() = 0;

	virtual
	~Norm() {}

	/*!
	 * Return the sum of all values
	 *
	 * This is helpful for debugging reasons;
	 */
public:
	virtual
	double reduceSum() const
	{
		SWEETErrorFatal("reduceSum(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceMaxAbs() const
	{
		SWEETErrorFatal("reduceMaxAbs(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceSum(int i) const
	{
		SWEETErrorFatal("reduceSum(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceMaxAbs(int i) const 
	{
		SWEETErrorFatal("reduceMaxAbs(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceMaxAbs(int i, int rnorm) const
	{
		SWEETErrorFatal("reduceMaxAbs(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceNormL1Grid(bool normalized = false) const
	{
		SWEETErrorFatal("reduceL1Grid(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceNormL2Grid(bool normalized = false) const
	{
		SWEETErrorFatal("reduceL2Grid(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceNormLinfGrid() const
	{
		SWEETErrorFatal("reduceLinfGrid(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceNormL1Grid(int i, bool normalized = false) const
	{
		SWEETErrorFatal("reduceL1Grid(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceNormL2Grid(int i, bool normalized = false) const
	{
		SWEETErrorFatal("reduceL2Grid(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

	virtual
	double reduceNormLinfGrid(int i) const
	{
		SWEETErrorFatal("reduceLinfGrid(...) is not implemented in sweet::Data::GenericContainer::Norm");
		return 0;
	}

};

}}}

#endif
