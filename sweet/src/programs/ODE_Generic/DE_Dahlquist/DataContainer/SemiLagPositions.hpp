#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DATACONTAINER_SEMILAGPOSITIONS_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DATACONTAINER_SEMILAGPOSITIONS_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <utility>

#if SWEET_MPI
#include <mpi.h>
#endif

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DataContainer {

class SemiLagPositions :
	public sweet::Data::GenericContainer::Base
{
public:
	// How many number of DoF arrays
	static constexpr int N = 1;

	sweet::Data::Vector::Vector<double> data;

public:
	SemiLagPositions()
	{
	}

public:
	~SemiLagPositions() override
	{
		clear();
	}

private:
	static inline
	SemiLagPositions& cast(sweet::Data::GenericContainer::Base &i_U)
	{
		return static_cast<SemiLagPositions&>(i_U);
	}

private:
	static inline
	const SemiLagPositions& cast(const sweet::Data::GenericContainer::Base &i_U)
	{
		return static_cast<const SemiLagPositions&>(i_U);
	}

private:
	static inline
	SemiLagPositions& cast_MPI(sweet::Data::GenericContainer::MPI &i_U)
	{
		return static_cast<SemiLagPositions&>(i_U);
	}

private:
	static inline
	const SemiLagPositions& cast_MPI(const sweet::Data::GenericContainer::MPI &i_U)
	{
		return static_cast<const SemiLagPositions&>(i_U);
	}

public:
	void swap(
		Base &i_U
	) override
	{
		data.swap(cast(i_U).data);
	}

public:
	void setup_like(
		const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const SemiLagPositions &i_d = static_cast<const SemiLagPositions&>(i_a);

		data.setup(i_d.data.size());
	}

#if 1
	Base* getNewDataContainer() const override
	{
		SemiLagPositions *retval = new SemiLagPositions;
		retval->setup_like(*this);

		return retval;
	}
#endif

public:
	void setup(
//		std::size_t i_numElements
	)
	{
		//data.setup(i_numElements);
		data.setup(1);
	}

public:
	void clear() override
	{
		data.clear();
	}

public:
	void op_setVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const SemiLagPositions &i_A = cast(i_a);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i];
	}

public:
	void op_addVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const SemiLagPositions &i_A = cast(i_a);
		for (int i = 0; i < N; i++)
			data[i] += i_A.data[i];
	}

public:
	void op_subVector(
			const sweet::Data::GenericContainer::Base &i_A_
	) override
	{
		const SemiLagPositions &i_A = cast(i_A_);
		for (int i = 0; i < N; i++)
			data[i] -= i_A.data[i];
	}

public:
	void op_setZero() override
	{
		for (int i = 0; i < N; i++)
			data[i] = 0;
	}

public:
	void op_setVectorPlusVector(
			const sweet::Data::GenericContainer::Base &i_A_,
			const sweet::Data::GenericContainer::Base &i_B_
	) override
	{
		const SemiLagPositions &i_A = cast(i_A_);
		const SemiLagPositions &i_B = cast(i_B_);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i] + i_B.data[i];
	}

public:
	void op_setVectorPlusScalarMulVector(
			const sweet::Data::GenericContainer::Base &i_A_,
			double i_scalar,
			const sweet::Data::GenericContainer::Base &i_B_
	) override
	{
		const SemiLagPositions &i_A = cast(i_A_);
		const SemiLagPositions &i_B = cast(i_B_);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i] + i_scalar*i_B.data[i];
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::Data::GenericContainer::Base &i_A_
		) override
	{
		const SemiLagPositions &i_A = cast(i_A_);
		for (int i = 0; i < N; i++)
			data[i] += i_scalar*i_A.data[i];
	}

public:
	void op_mulScalar(
			double i_scalar
		) override
	{
		for (int i = 0; i < N; i++)
			data[i] *= i_scalar;
	}


#if SWEET_MPI
	void mpiBcast(MPI_Comm &i_mpi_comm) override
	{
		for (int i = 0; i < N; i++)
		{
			MPI_Bcast(
					////data[i].data,
					////data[i].numberOfElements,
					&data[i],
					data.numberOfElements,
					MPI_DOUBLE,
					0,
					i_mpi_comm
				);
		}
	}

	void mpiReduce(
			const sweet::Data::GenericContainer::MPI &i_a,
			MPI_Comm &i_mpi_comm
	) override
	{
		const SemiLagPositions &i_A = cast_MPI(i_a);

		for (int i = 0; i < N; i++)
		{
			MPI_Reduce(
					///i_A.data[i].data,
					///data[i].data,
					///data[i].numberOfElements,
					&(i_A.data.data[i]),
					&data[i],
					data.numberOfElements,
					MPI_DOUBLE,
					MPI_SUM,
					0,
					i_mpi_comm
				);
		}
	}

	void mpiReduceAll(
		const sweet::Data::GenericContainer::MPI &i_a,
		MPI_Comm &i_mpi_comm
	) override
	{
		const SemiLagPositions &i_A = cast_MPI(i_a);

		for (int i = 0; i < N; i++)
		{
			MPI_Allreduce(
					////i_A.data[i].data,
					////data[i].data,
					////data[i].numberOfElements,
					&(i_A.data.data[i]),
					&(data[i]),
					data.numberOfElements,
					MPI_DOUBLE,
					MPI_SUM,
					i_mpi_comm
				);
		}
	}
#endif


	double reduceSum()	const override
	{
		std::complex<double> sum = 0;

		for (int i = 0; i < N; i++)
			sum += data[i];

		return std::abs(sum);
	}

};

}}}

#endif
