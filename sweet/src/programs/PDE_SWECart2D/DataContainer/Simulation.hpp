#ifndef PROGRAMS_PDE_SWECART2D_DATACONTAINER_SIMULATION_HPP
#define PROGRAMS_PDE_SWECART2D_DATACONTAINER_SIMULATION_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <utility>

#if SWEET_MPI
#include <mpi.h>
#endif

namespace PDE_SWECart2D {
namespace DataContainer {

class Simulation :
	public sweet::Data::GenericContainer::Base
{
public:
	// How many number of DoF arrays
	static constexpr int N = 3;
	std::vector<std::string> var_names = {"prog_h_pert", "prog_u", "prog_v"};

	sweet::Data::Cart2D::DataSpectral data[N];

	// Note, that these references don't increase the size of the class
	sweet::Data::Cart2D::DataSpectral& h_pert = data[0];
	sweet::Data::Cart2D::DataSpectral& u = data[1];
	sweet::Data::Cart2D::DataSpectral& v = data[2];

public:
	Simulation()
	{
	}

public:
	~Simulation() override
	{
		clear();
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
		const Simulation &i_d = static_cast<const Simulation&>(i_a);

		for (int i = 0; i < N; i++)
			data[i].setup(i_d.data[i].cart2DDataConfig);
	}

#if 1
	Base* getNewDataContainer() const override
	{
		Simulation *retval = new Simulation;
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
			data[i].setup(i_cart2DData_Config);
	}

public:
	void clear() override
	{
		for (int i = 0; i < N; i++)
			data[i].clear();
	}

	void setData(
								sweet::Data::Cart2D::DataSpectral &i_h_pert,
								sweet::Data::Cart2D::DataSpectral &i_u,
								sweet::Data::Cart2D::DataSpectral &i_v
							)
	{
		data[0] = i_h_pert;
		data[1] = i_u;
		data[2] = i_v;
	}


public:
	void op_setVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const Simulation &i_A = cast(i_a);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i];
	}

public:
	void op_addVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const Simulation &i_A = cast(i_a);
		for (int i = 0; i < N; i++)
			data[i] += i_A.data[i];
	}

public:
	void op_subVector(
			const sweet::Data::GenericContainer::Base &i_A_
	) override
	{
		const Simulation &i_A = cast(i_A_);
		for (int i = 0; i < N; i++)
			data[i] -= i_A.data[i];
	}

public:
	void op_setZero() override
	{
		for (int i = 0; i < N; i++)
			data[i].spectral_setZero();
	}

public:
	void op_setVectorPlusVector(
			const sweet::Data::GenericContainer::Base &i_A_,
			const sweet::Data::GenericContainer::Base &i_B_
	) override
	{
		const Simulation &i_A = cast(i_A_);
		const Simulation &i_B = cast(i_B_);
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
		const Simulation &i_A = cast(i_A_);
		const Simulation &i_B = cast(i_B_);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i] + i_scalar*i_B.data[i];
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::Data::GenericContainer::Base &i_A_
		) override
	{
		const Simulation &i_A = cast(i_A_);
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
					data[i].spectral_space_data,
					data[i].cart2DDataConfig->spectral_array_data_number_of_elements,
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
		const Simulation &i_A = cast_MPI(i_a);

		for (int i = 0; i < N; i++)
		{
			MPI_Reduce(
					i_A.data[i].spectral_space_data,
					data[i].spectral_space_data,
					data[i].cart2DDataConfig->spectral_array_data_number_of_elements,
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
		const Simulation &i_A = cast_MPI(i_a);

		for (int i = 0; i < N; i++)
		{
			MPI_Allreduce(
					i_A.data[i].spectral_space_data,
					data[i].spectral_space_data,
					data[i].cart2DDataConfig->spectral_array_data_number_of_elements,
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
			sum += data[i].spectral_reduce_sum_quad();

		return std::abs(sum);
	}

	double reduceMaxAbs()	const override
	{
		double e = -1;
		for (int i = 0; i < N; i++)
			e = std::max( e,
					data[i].spectral_reduce_max_abs());
		return e;
	}

	double reduceSum(int i) const override
	{
		return std::abs(data[i].spectral_reduce_sum_quad());
	}

	double reduceMaxAbs(int i) const override
	{
		return data[i].spectral_reduce_max_abs();
	}

	double reduceMaxAbs(int i, int rnorm) const override
	{
		return data[i].spectral_reduce_max_abs(rnorm);
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
		double d = data[i].toGrid().grid_reduce_norm1();

		if (normalized)
		{
			std::size_t resx_data = data[i].cart2DDataConfig->grid_data_size[0];
			std::size_t resy_data = data[i].cart2DDataConfig->grid_data_size[1];
			d /= (resx_data * resy_data);
		}

		return d;
	}

	double reduceNormL2Grid(int i, bool normalized = false) const override
	{
		double d = data[i].toGrid().grid_reduce_norm2();

		if (normalized)
		{
			std::size_t resx_data = data[i].cart2DDataConfig->grid_data_size[0];
			std::size_t resy_data = data[i].cart2DDataConfig->grid_data_size[1];
			d /= std::sqrt(resx_data * resy_data);
		}

		return d;
	}

	double reduceNormLinfGrid(int i) const override
	{
		///return reduceMaxAbs(i);
		return data[i].toGrid().grid_reduce_max_abs();
	}

	void pad_zeros(const sweet::Data::GenericContainer::Base &i_data_) override
	{
		const Simulation &i_data = cast(i_data_);
		for (int i = 0; i < N; i++)
			data[i] = data[i].pad_zeros(i_data.data[i]);

	}

	void restrict(const sweet::Data::GenericContainer::Base &i_data_) override
	{
		const Simulation &i_data = cast(i_data_);
		for (int i = 0; i < N; i++)
			data[i] = data[i].restrict(i_data.data[i]);
	}

	void serialize(std::complex<double> *i_data) override
	{
		int s = data[0].cart2DDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&data[i].spectral_space_data[0], &data[i].spectral_space_data[s], &i_data[i * s]);
	};

	void deserialize(std::complex<double> *i_data) override
	{
		int s = data[0].cart2DDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&i_data[i * s], &i_data[(i + 1) * s], &data[i].spectral_space_data[0]);
	};


public:
	bool fileLoad(
		const std::string &i_filename,
		const std::string &i_filemode
	)
	{
		if (i_filename.find("prog_h_pert") != std::string::npos)
			h_pert.loadDataFromFile(i_filename, i_filemode);
		else if (i_filename.find("prog_u") != std::string::npos)
			u.loadDataFromFile(i_filename, i_filemode);
		else if (i_filename.find("prog_v") != std::string::npos)
			v.loadDataFromFile(i_filename, i_filemode);
		else
			SWEETErrorFatal("Invalid input filename: " + i_filename);

		return true;
	}

};

}}

#endif
