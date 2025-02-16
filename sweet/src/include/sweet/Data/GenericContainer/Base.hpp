/*
 * DESolver_DataContainer_Base.hpp
 */

#ifndef INCLUDE_SWEET_DATA_GENERICCONTAINER_BASE_HPP
#define INCLUDE_SWEET_DATA_GENERICCONTAINER_BASE_HPP

#if SWEET_MPI
#include <mpi.h>
#endif

#include <complex>

#include <sweet/Error/Fatal.hpp>
#include "MPI.hpp"
#include "Norm.hpp"

namespace sweet {
namespace Data {
namespace GenericContainer {

class Base:
	public MPI,
	public Norm
{
public:
	enum DATA_ENUM {
		DATA_SIMULATION = 0,
		DATA_SEMI_LAGRANGIAN_POSITIONS = 1,
		DATA_TOPOGRAPHY = 2,

		DATA_USER_DEFINED_0 = 1000,
		DATA_USER_DEFINED_1 = 1001,
		DATA_USER_DEFINED_2 = 1002,
		/* ... */

	};
public:
	Base()
	{
	}

public:
	virtual
	void swap(Base &i_U) = 0;

public:
	virtual
	void op_setZero() = 0;

public:
	virtual
	void op_setVector(
			const GenericContainer::Base &i_a
		) = 0;

public:
	virtual
	void op_addVector(
			const GenericContainer::Base &i_a
		) = 0;

public:
	virtual
	void op_subVector(
			const GenericContainer::Base &i_a
		) = 0;

public:
	virtual
	void op_setVectorPlusVector(
			const GenericContainer::Base &i_a,
			const GenericContainer::Base &i_b
	) = 0;

public:
	virtual
	void op_setVectorPlusScalarMulVector(
			const GenericContainer::Base &i_a,
			double i_scalar,
			const GenericContainer::Base &i_b
		) = 0;

public:
	virtual
	void op_addScalarMulVector(
			double i_scalar,
			const GenericContainer::Base &i_a
		) = 0;

public:
	virtual
	void op_mulScalar(
			double i_scalar
		) = 0;

public:
	virtual
	void setup_like(
			const GenericContainer::Base &i_a
	) = 0;

#if 1
	virtual
	Base* getNewDataContainer() const = 0;
#endif

	virtual
	void clear() = 0;

	virtual
	~Base() {}

public:
	virtual
	void pad_zeros(const sweet::Data::GenericContainer::Base &i_data)
	{
		SWEETErrorFatal("TODO: pad_zeros(...) is not implemented");
	}

	virtual
	void restrict(const sweet::Data::GenericContainer::Base &i_data)
	{
		SWEETErrorFatal("TODO: restrict(...) is not implemented");
	}


};

}}}

#endif
