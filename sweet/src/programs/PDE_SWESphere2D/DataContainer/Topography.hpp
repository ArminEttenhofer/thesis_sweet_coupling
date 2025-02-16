#ifndef PROGRAMS_PDE_SWESPHERE2D_DATACONTAINER_TOPOGRAPHY_HPP
#define PROGRAMS_PDE_SWESPHERE2D_DATACONTAINER_TOPOGRAPHY_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <utility>

#if SWEET_MPI
#include <mpi.h>
#endif

namespace PDE_SWESphere2D {
namespace DataContainer {

class Topography:
	public sweet::Data::GenericContainer::Base
{
public:
	// How many number of DoF arrays
	static constexpr int N = 1;
	std::vector<std::string> var_names = {"topography"};

	sweet::Data::Sphere2D::DataSpectral data[N];
	sweet::Data::Sphere2D::DataGrid data_grid[N];

	// Note, that these references don't increase the size of the class
	sweet::Data::Sphere2D::DataSpectral& topography = data[0];
	sweet::Data::Sphere2D::DataGrid& topography_grid = data_grid[0];

public:
	Topography()
	{
	}

public:
	~Topography() override
	{
		clear();
	}

private:
	static inline
	Topography& cast(sweet::Data::GenericContainer::Base &i_U)
	{
		return static_cast<Topography&>(i_U);
	}

private:
	static inline
	const Topography& cast(const sweet::Data::GenericContainer::Base &i_U)
	{
		return static_cast<const Topography&>(i_U);
	}

private:
	static inline
	Topography& cast_MPI(sweet::Data::GenericContainer::MPI &i_U)
	{
		return static_cast<Topography&>(i_U);
	}

private:
	static inline
	const Topography& cast_MPI(const sweet::Data::GenericContainer::MPI &i_U)
	{
		return static_cast<const Topography&>(i_U);
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
		const Topography &i_d = static_cast<const Topography&>(i_a);

		for (int i = 0; i < N; i++)
			data[i].setup(i_d.data[i].sphere2DDataConfig);
	}

#if 1
	Base* getNewDataContainer() const override
	{
		Topography *retval = new Topography;
		retval->setup_like(*this);

		return retval;
	}
#endif

public:
	void setup(
		const sweet::Data::Sphere2D::Config *i_sphere2DData_Config
	)
	{
		for (int i = 0; i < N; i++)
		{
			data[i].setup(i_sphere2DData_Config);
			data_grid[i].setup(i_sphere2DData_Config);

			data[i].spectral_setZero();
			data_grid[i].grid_setZero();
		}
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
		const Topography &i_A = cast(i_a);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i];
	}

public:
	void op_addVector(
			const sweet::Data::GenericContainer::Base &i_a
	) override
	{
		const Topography &i_A = cast(i_a);
		for (int i = 0; i < N; i++)
			data[i] += i_A.data[i];
	}

public:
	void op_subVector(
			const sweet::Data::GenericContainer::Base &i_A_
	) override
	{
		const Topography &i_A = cast(i_A_);
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
		const Topography &i_A = cast(i_A_);
		const Topography &i_B = cast(i_B_);
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
		const Topography &i_A = cast(i_A_);
		const Topography &i_B = cast(i_B_);
		for (int i = 0; i < N; i++)
			data[i] = i_A.data[i] + i_scalar*i_B.data[i];
	}

public:
	void op_addScalarMulVector(
			double i_scalar,
			const sweet::Data::GenericContainer::Base &i_A_
		) override
	{
		const Topography &i_A = cast(i_A_);
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
					data[i].sphere2DDataConfig->spectral_array_data_number_of_elements*2,
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
		const Topography &i_A = cast_MPI(i_a);

		for (int i = 0; i < N; i++)
		{
			MPI_Reduce(
					i_A.data[i].spectral_space_data,
					data[i].spectral_space_data,
					data[i].sphere2DDataConfig->spectral_array_data_number_of_elements*2,
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
		const Topography &i_A = cast_MPI(i_a);

		for (int i = 0; i < N; i++)
		{
			MPI_Allreduce(
					i_A.data[i].spectral_space_data,
					data[i].spectral_space_data,
					data[i].sphere2DDataConfig->spectral_array_data_number_of_elements*2,
					MPI_DOUBLE,
					MPI_SUM,
					i_mpi_comm
				);
		}
	}
#endif

public:
	bool fileLoad(
		const std::string &i_filename,
		const std::string &i_filemode
	)
	{
		if (i_filename.find("topography") != std::string::npos)
			topography.loadDataFromFile(i_filename, i_filemode);
		else
			SWEETErrorFatal("Invalid input filename: " + i_filename);

		return true;
	}




};

}}

#endif
