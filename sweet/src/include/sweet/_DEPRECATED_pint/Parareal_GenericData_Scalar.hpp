/*
 * Parareal_GenericData_Scalar.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_SCALAR_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_SCALAR_HPP

#include <assert.h>
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>

namespace sweet {
namespace DEPRECATED_pint {


template <int N>
class Parareal_GenericData_Scalar :
		public sweet::DEPRECATED_pint::Parareal_GenericData
{
	class DataContainer_Scalar :
			public sweet::DEPRECATED_pint::Parareal_GenericData::DataContainer<double>
	{

	public:

		DataContainer_Scalar()
		{
			nb_fields = N;
			simfields = new double[N];
			for (int i = 0; i < N; i++)
				simfields[i] = 0.;
		};

		DataContainer_Scalar(
				double i_data
		)
		{
			nb_fields = N;
			simfields = new double[N];
			for (int i = 0; i < N; i++)
				simfields[i] = i_data;
		};

		DataContainer_Scalar(
				double* i_data
		)
		{
			nb_fields = N;
			simfields = new double[N];
			for (int i = 0; i < N; i++)
				simfields[i] = i_data[i];
		};


		DataContainer_Scalar(DataContainer_Scalar &i_data)
		{
			time = i_data.time;
			level = i_data.level;
			nb_fields = i_data.nb_fields;
			for (int i = 0; i < N; i++)
				simfields[i] = i_data.simfields[i];
		};

		~DataContainer_Scalar()
		{
			if (simfields)
			{
				delete [] simfields;
				simfields = nullptr;
			}
		}
	};


public:

	DataContainer<double>* data = nullptr;

public:
	DataContainer<double>* get_pointer_to_data_Scalar() const override
	{
		return data;
	};

	void dataArrays_2_GenericData_Scalar(
						double &u
						) override
	{
		get_pointer_to_data_Scalar()->simfields[0] = u;
	}

	void GenericData_Scalar_2_dataArrays(
						double &u
						) override
	{
		u = get_pointer_to_data_Scalar()->simfields[0];
	}


public:

	Parareal_GenericData_Scalar():
		Parareal_GenericData()
	{
		////allocate_data();
	}

	Parareal_GenericData_Scalar(Parareal_GenericData_Scalar &i_data)
	{
		*(data) = *(i_data.get_pointer_to_data_Scalar());
		for (int i = 0; i < N; i++)
			data->simfields[i] = i_data.get_pointer_to_data_Scalar()->simfields[i];
	};

	Parareal_GenericData_Scalar& operator=(const Parareal_GenericData &i_data)	override
	{
		*(data) = *(i_data.get_pointer_to_data_Scalar());
		for (int i = 0; i < N; i++)
			data->simfields[i] = i_data.get_pointer_to_data_Scalar()->simfields[i];
		return *this;
	};

	~Parareal_GenericData_Scalar()
	{
		free_data();
	};

	void allocate_data()	override
	{
		data = new DataContainer_Scalar();
	}
	
	void free_data()	override
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
	}
	

#if 0
#error "M@J: This assignment doesn't make sense"
	/**
	 * Setup data
	 */
	void setup(
			double i_data
	)
	{
		allocate_data();
		data = i_data;
	}
#endif

/////#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	// size in bytes (for MPI)
	// size of each simfield of data
	std::size_t size()	override
	{
		return 1;
	}

	void serialize(double *i_data)	override
	{
		for (int i = 0; i < N; i++)
		{
			///std::memcpy(data + i, &data->simfields[i], 1);
			std::copy(&data->simfields[i], &data->simfields[i + 1], i_data);
		}
	};

	void deserialize(double *i_data)	override
	{
		for (int i = 0; i < N; i++)
		{
			///std::memcpy(&data->simfields[i], data + i, 1);
			std::copy(&i_data[i], &i_data[i + 1], &data->simfields[i]);
		}
	};
#endif


	void set_time(double i_time)	override
	{
		data->set_time(i_time);
	}


	double spectral_reduce_maxAbs()	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					std::abs(data->simfields[k]));
		return e;
	}

	double spectral_reduce_maxAbs(std::size_t rnorm)	override
	{
		return spectral_reduce_maxAbs();
	}

	double grid_reduce_maxAbs()	override
	{
		return spectral_reduce_maxAbs();
	}


	double grid_reduce_norm1()	override
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += std::abs(data->simfields[k]);
		return e;
	}

	double grid_reduce_norm2()	override
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += data->simfields[k] * data->simfields[k];
		return std::sqrt(e);
	}


	bool check_for_nan()	override
	{
		bool found_nan = false;
		for (int i = 0; i < N; i++)
			if (std::isnan(data->simfields[i]))
			{
				found_nan = true;
				break;
			}
		return found_nan;
	}

	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)	override
	{
#if SWEET_PARAREAL
		SWEET_ASSERT(data->time == i_data.get_pointer_to_data_Scalar()->time);
#endif
		SWEET_ASSERT(data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		for (int i = 0; i < N; i++)
			data->simfields[i] += i_data.get_pointer_to_data_Scalar()->simfields[i];

		return *this;
	}


	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)	override
	{
#if SWEET_PARAREAL
		SWEET_ASSERT(data->time == i_data.get_pointer_to_data_Scalar()->time);
#endif
		SWEET_ASSERT(data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		for (int i = 0; i < N; i++)
			data->simfields[i] -= i_data.get_pointer_to_data_Scalar()->simfields[i];

		return *this;
	}

	Parareal_GenericData& operator*=( double v)	override
	{

		for (int i = 0; i < N; i++)
			data->simfields[i] *= v;

		return *this;
	}

	// Nothing to do
	void restrict(const Parareal_GenericData& i_data)	override
	{
	}

	// Nothing to do
	void pad_zeros(const Parareal_GenericData& i_data)	override
	{
	}

	void grid_print()	override
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << ": " << data->simfields[i] << std::endl;
		}
	}

	void spectral_print()	override
	{
		grid_print();
	}


};

}}

#endif
