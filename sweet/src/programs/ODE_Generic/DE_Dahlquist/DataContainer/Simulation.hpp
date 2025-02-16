/*
 * DataContainer.hpp
 */

#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DATACONTAINER_SIMULATION_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_DATACONTAINER_SIMULATION_HPP

#include <utility>
#include <complex>
#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/GenericContainer/MPI.hpp>
////#include <sweet/Data/Scalar/Config.hpp>
#include <sweet/Dict/Dict.hpp>
#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/Error/Base.hpp>

///#include <programs/ODE_Generic/DE_Dahlquist/DataContainer/Config.hpp>
///
///class Config;

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace DataContainer {

/*!
 * \brief Data Container for Dahlquist ODE equation
 */
class Simulation :
	public sweet::Data::GenericContainer::Base
{
public:
	sweet::Error::Base error;

	static constexpr int N = 1;
	std::vector<std::string> var_names = {"prog_u"};

	typedef std::complex<double> T;

	T data[N];

	// Note, that these references don't increase the size of the class
	T& U = data[0];


public:
	Simulation()
	{
	}

public:
	~Simulation() override
	{
		clear();
	}

public:
	Simulation(
		const Simulation &i_src
	)
	{
		op_setVector(i_src);
	}

private:
	static inline
	Simulation& cast(sweet::Data::GenericContainer::Base &i_U)
	{
		return static_cast<Simulation&>(i_U);
	}

private:
	static inline
	const Simulation& cast(const sweet::Data::GenericContainer::Base &i_U)
	{
		return static_cast<const Simulation&>(i_U);
	}

private:
	static inline
	Simulation& cast_MPI(sweet::Data::GenericContainer::MPI &i_U)
	{
		return static_cast<Simulation&>(i_U);
	}

private:
	static inline
	const Simulation& cast_MPI(const sweet::Data::GenericContainer::MPI &i_U)
	{
		return static_cast<const Simulation&>(i_U);
	}


public:
	void swap(
		Base &i_U_
	) override
	{
		std::swap(cast(i_U_).U, U);
	}

public:
	void setup_like(
		const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		// Nothing really to do since we know the data size
		U = -666;
	}

	Base* getNewDataContainer() const override
	{
		SWEETErrorFatal("DEPRECATED: Data should be allocated via DESolver_DataContainer::Config_Base");
		return nullptr;
	}

public:
	void setup()
	{
		U = -1;
	}

//////#if SWEET_XBRAID
//////public:
//////	void setup(
//////			const sweet::Data::Scalar::Config *i_config
//////		)
//////	{
//////		setup();
//////	}
//////#endif

public:
	void clear() override
	{
		U = 0;
	}

public:
	void op_setVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const Simulation &i_A = cast(i_a);
		U = i_A.U;
	}

public:
	void op_addVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const Simulation &i_A = cast(i_a);
		U += i_A.U;
	}

public:
	void op_subVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const Simulation &i_A = cast(i_a);
		U -= i_A.U;
	}

public:
	void op_setZero() override
	{
		U = 0;
	}

public:
	void op_setVectorPlusVector(
			const sweet::Data::GenericContainer::Base &i_a,
			const sweet::Data::GenericContainer::Base &i_b
	) override
	{
		const Simulation &i_A = cast(i_a);
		const Simulation &i_B = cast(i_b);
		U = i_A.U + i_B.U;
	}

public:
	void op_setVectorPlusScalarMulVector(
			const sweet::Data::GenericContainer::Base &i_a,
			double i_scalar,
			const sweet::Data::GenericContainer::Base &i_b
	) override
	{
		const Simulation &i_A = cast(i_a);
		const Simulation &i_B = cast(i_b);
		U = i_A.U + i_scalar*i_B.U;
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::Data::GenericContainer::Base &i_a
		) override
	{
		const Simulation &i_A = cast(i_a);
		U += i_scalar*i_A.U;
	}

public:
	void op_mulScalar(
			double i_scalar
		) override
	{
		U *= i_scalar;
	}

public:
	void loadSolutionFromFile(
					const char *i_filename
	)
	{
	}



#if SWEET_MPI
	void mpiBcast(MPI_Comm &i_mpi_comm) override
	{
		MPI_Bcast(
				&U,
				2,
				MPI_DOUBLE,
				0,
				i_mpi_comm
			);
	}

	void mpiReduce(
			const sweet::Data::GenericContainer::MPI &i_a,
			MPI_Comm &i_mpi_comm
	) override
	{
		const Simulation &i_A = cast_MPI(i_a);

		MPI_Reduce(
				&(i_A.U),
				&U,
				2,
				MPI_DOUBLE,
				MPI_SUM,
				0,
				i_mpi_comm
			);
	}

	void mpiReduceAll(
		const sweet::Data::GenericContainer::MPI &i_a,
		MPI_Comm &i_mpi_comm
	) override
	{
		const Simulation &i_A = cast_MPI(i_a);

		MPI_Allreduce(
				&(i_A.U),
				&U,
				2,
				MPI_DOUBLE,
				MPI_SUM,
				i_mpi_comm
			);
	}
#endif


public:
	double reduceSum()	const override
	{
		return std::abs(U);
	}

	double reduceMaxAbs()	const override
	{
		return std::abs(U);
	}

	double reduceSum(int i)		const override
	{
		return std::abs(U);
	}

	double reduceMaxAbs(int i)	const override
	{
		return std::abs(U);
	}

	double reduceMaxAbs(int i, int rnorm)	const override
	{
		return std::abs(U);
	}

	double reduceNormL1Grid(bool normalized = false) const override
	{
		double norm = 0;
		for (int i = 0; i < N; i++)
			norm += reduceNormL1Grid(i, normalized);
		return norm;
	}

	double reduceNormL2Grid(bool normalized = false) const override
	{
		double norm = 0;
		for (int i = 0; i < N; i++)
		{
			double a = reduceNormL2Grid(i, normalized);
			norm += a * a;
		}
		return std::sqrt(norm);
	}

	double reduceNormLinfGrid() const override
	{
		double norm = 0;
		for (int i = 0; i < N; i++)
			norm = std::max(norm, reduceNormLinfGrid(i));
		return norm;
	}

	double reduceNormL1Grid(int i, bool normalized = false) const override
	{
		return std::abs(U);
	}

	double reduceNormL2Grid(int i, bool normalized = false) const override
	{
		return std::abs(U);
	}

	double reduceNormLinfGrid(int i) const override
	{
		return std::abs(U);
	}

	void pad_zeros(const sweet::Data::GenericContainer::Base &i_data_) override
	{
		op_setVector(i_data_);
	}

	void restrict(const sweet::Data::GenericContainer::Base &i_data_) override
	{
		op_setVector(i_data_);
	}

	void serialize(std::complex<double> *i_data) override
	{
		for (int i = 0; i < N; i++)
			std::copy(&data[i], &data[i + 1], &i_data[i]);
	};

	void deserialize(std::complex<double> *i_data) override
	{
		for (int i = 0; i < N; i++)
			std::copy(&i_data[i], &i_data[i + 1], &data[i]);
	};

/////	bool fileSave(
/////			const std::string &i_filename
/////	)	const
/////	{
/////		//SWEET_ASSERT_MSG(i_splitFiles == false, "Splitting of files not supported so far");
/////
/////		sweet::Dict::Dict dict;
/////
/////		// Setup Header
/////		dict.set("sweetMagicCode", "SWEET1505");
/////		dict.set("dataType", "VectorComplex");
/////
/////		/*
/////		 * Instead of storing just a single scalar, we store it as a vector
/////		 * This makes it more generic
/////		 */
/////		sweet::Dict::TypesArrayND<1,std::complex<double>> vector;
/////		vector.resize(1);
/////		vector.set(0, U);
/////
/////		// Store data
/////		dict.set("data", vector);
/////		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_BOOLEAN(dict);
/////
/////		// Write to file
/////		dict.fileSave(i_filename);
/////		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_BOOLEAN(dict);
/////
/////		return true;
/////	}

public:
	bool fileLoad(const std::string &i_filename)
	{

		sweet::Dict::Dict dict;
		dict.set("sweetMagicCode", "SWEET1505");
		dict.set("dataType", "VectorComplex");

		dict.fileLoad(i_filename);

		sweet::Dict::TypesArrayND<1,std::complex<double>> vector;
		vector.resize(1);

		dict.get("data", vector);

		U = vector.get(0);

		return true;
	}



};

}}}

#endif
