/*
 * DESolver_DataContainer_Base.hpp
 */

#ifndef INCLUDE_SWEET_DATA_GENERICCONTAINER_MPI_HPP
#define INCLUDE_SWEET_DATA_GENERICCONTAINER_MPI_HPP

#if SWEET_MPI
#include <mpi.h>
#endif

#include <complex>

#include <sweet/Error/Fatal.hpp>

namespace sweet {
namespace Data {
namespace GenericContainer {

class MPI
{
public:
	MPI()
	{
	}

#if 1
	virtual
	MPI* getNewDataContainer() const = 0;
#endif

	virtual
	void clear() = 0;

	virtual
	~MPI() {}

#if SWEET_MPI
	virtual
	void mpiBcast(MPI_Comm &i_mpi_comm)
	{
		SWEETErrorFatal("TODO: mpiBcast(...) is not implemented");
	}

	virtual
	void mpiReduce(
			const GenericContainer::MPI &i_a,
			MPI_Comm &i_mpi_comm
	)
	{
		SWEETErrorFatal("TODO: mpiReduce(...) is not implemented");
	}

	virtual
	void mpiReduceAll(
			const GenericContainer::MPI &i_a,
			MPI_Comm &i_mpi_comm
	)
	{
		SWEETErrorFatal("TODO: mpiReduceAll(...) is not implemented");
	}
#endif

public:
	virtual
	void serialize(std::complex<double> *i_data)
	{
		SWEETErrorFatal("serialize(...) is not implemented in sweet::Data::GenericContainer::MPI");
	}

	virtual
	void deserialize(std::complex<double> *i_data)
	{
		SWEETErrorFatal("deserialize(...) is not implemented in sweet::Data::GenericContainer::MPI");
	}

};

}}}

#endif
