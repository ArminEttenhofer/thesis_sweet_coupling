#ifndef PROGRAMS_PDE_SWECART2D_DATACONTAINER_SEMILAGPOSITIONS_HPP
#define PROGRAMS_PDE_SWECART2D_DATACONTAINER_SEMILAGPOSITIONS_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <utility>

#if SWEET_MPI
#include <mpi.h>
#endif

namespace PDE_SWECart2D {
namespace DataContainer {

class SemiLagPositions :
	public sweet::Data::GenericContainer::Base
{
public:
	// How many number of DoF arrays
	static constexpr int N = 2;

	sweet::Data::Vector::Vector<double> data[N];

	// Note, that these references don't increase the size of the class
	sweet::Data::Vector::Vector<double> &pos_x = data[0];
	sweet::Data::Vector::Vector<double> &pos_y = data[1];

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
		for (int i = 0; i < N; i++)
			data[i].swap(cast(i_U).data[i]);
	}

public:
	void setup_like(
		const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const SemiLagPositions &i_d = static_cast<const SemiLagPositions&>(i_a);

		for (int i = 0; i < N; i++)
			data[i].setup(i_d.data[i].numberOfElements);
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
		const sweet::Data::Cart2D::Config *i_cart2DData_Config
	)
	{
		for (int i = 0; i < N; i++)
			data[i].setup(i_cart2DData_Config->grid_number_elements);
	}

public:
	void clear() override
	{
		for (int i = 0; i < N; i++)
			data[i].clear();
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
			data[i].set_zero();
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
					data[i].data,
					data[i].numberOfElements,
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
					i_A.data[i].data,
					data[i].data,
					data[i].numberOfElements,
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
					i_A.data[i].data,
					data[i].data,
					data[i].numberOfElements,
					MPI_DOUBLE,
					MPI_SUM,
					i_mpi_comm
				);
		}
	}
#endif


};

}}

#endif
