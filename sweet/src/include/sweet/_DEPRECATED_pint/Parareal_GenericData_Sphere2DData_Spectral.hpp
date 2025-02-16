/*
 * Parareal_GenericData_Sphere2DData_Spectral.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_SPHERE2DDATA_SPECTRAL_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_SPHERE2DDATA_SPECTRAL_HPP

#include <assert.h>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>

namespace sweet {
namespace DEPRECATED_pint {

template <int N>
class Parareal_GenericData_Sphere2DData_Spectral :
		public sweet::DEPRECATED_pint::Parareal_GenericData
{
	class DataContainer_Sphere2DData_Spectral :
			public sweet::DEPRECATED_pint::Parareal_GenericData::DataContainer<sweet::Data::Sphere2D::DataSpectral*>
	{

	public:

		DataContainer_Sphere2DData_Spectral(const sweet::Data::Sphere2D::Config* i_sphere2DDataConfig)
		{
			nb_fields = N;
			simfields = new sweet::Data::Sphere2D::DataSpectral*[N];
			for (int i = 0; i < N; i++)
				simfields[i] = new sweet::Data::Sphere2D::DataSpectral(i_sphere2DDataConfig);
		};

		DataContainer_Sphere2DData_Spectral(
				sweet::Data::Sphere2D::DataSpectral* i_simfields[N]
		)
		{
			nb_fields = N;
			simfields = new sweet::Data::Sphere2D::DataSpectral*[N];
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_simfields[i]);
		};

		DataContainer_Sphere2DData_Spectral(DataContainer_Sphere2DData_Spectral &i_data)
		{
			nb_fields = N;
			level = i_data.level;
			simfields = new sweet::Data::Sphere2D::DataSpectral*[N];
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_data.simfields[i]);
		};

		DataContainer_Sphere2DData_Spectral& operator=(const DataContainer_Sphere2DData_Spectral &i_data)
		{
			nb_fields = N;
			level = i_data.level;
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_data.simfields[i]);
			return *this;
		};


		~DataContainer_Sphere2DData_Spectral()
		{
			for (int i = 0; i < N; ++i)
				if (simfields[i])
					delete simfields[i];
			if (simfields)
				delete [] simfields;
		}

	};

public:

	DataContainer<sweet::Data::Sphere2D::DataSpectral*>* data = nullptr;

public:
	DataContainer<sweet::Data::Sphere2D::DataSpectral*>* get_pointer_to_data_Sphere2DData_Spectral() const override
	{
		return data;
	};


public:

	Parareal_GenericData_Sphere2DData_Spectral():
		sweet::DEPRECATED_pint::Parareal_GenericData()
	{
		///allocate_data();
	}


	Parareal_GenericData_Sphere2DData_Spectral(Parareal_GenericData_Sphere2DData_Spectral &i_data)
	{
		*(data) = *(i_data.get_pointer_to_data_Sphere2DData_Spectral());
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) = *(i_data.get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]);
		sphere2DDataConfig = i_data.sphere2DDataConfig;
	};

	Parareal_GenericData_Sphere2DData_Spectral& operator=(const Parareal_GenericData &i_data)
	override
	{
		*(data) = *(i_data.get_pointer_to_data_Sphere2DData_Spectral());
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) = *(i_data.get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]);
		sphere2DDataConfig = i_data.sphere2DDataConfig;
		return *this;
	};


	~Parareal_GenericData_Sphere2DData_Spectral()
	{
		free_data();
	};

	void allocate_data()
	override
	{
		data = new DataContainer_Sphere2DData_Spectral(sphere2DDataConfig);
	}

	void free_data()
	override
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
		
	}


////#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	// size in bytes (for MPI)
	// size of each simfield of data
	std::size_t size()
	{
		return N * data->simfields[0]->sphere2DDataConfig->spectral_array_data_number_of_elements;
	}

	void serialize(std::complex<double> *i_data)
	{
		int s = data->simfields[0]->sphere2DDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&data->simfields[i]->spectral_space_data[0], &data->simfields[i]->spectral_space_data[s], &i_data[i * s]);
	};

	void deserialize(std::complex<double> *i_data)
	{
		int s = data->simfields[0]->sphere2DDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&i_data[i * s], &i_data[(i + 1) * s], &data->simfields[i]->spectral_space_data[0]);
	};

#endif


	void set_time(double i_time)
	override
	{
		data->set_time(i_time);
	}

	double spectral_reduce_maxAbs()
	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					data->simfields[k]->spectral_reduce_max_abs());
		return e;
	}

	double spectral_reduce_maxAbs(std::size_t rnorm)
	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					data->simfields[k]->spectral_reduce_max_abs(rnorm));
		return e;
	}

	double grid_reduce_maxAbs()
	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					data->simfields[k]->toGrid().grid_reduce_max_abs());
		return e;
	}


	double grid_reduce_norm1()
	override
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += data->simfields[k]->toGrid().grid_reduce_norm1();
		return e;
	}

	double grid_reduce_norm2()
	override
	{
		double e = 0;
		for (int k = 0; k < N; k++)
		{
			double n = data->simfields[k]->toGrid().grid_reduce_norm2();
			e += n * n;
		}
		return std::sqrt(e);
	}

	bool check_for_nan()
	override
	{
		bool found_nan = false;

		int size_n = data->simfields[0]->sphere2DDataConfig->spectral_modes_n_max;
		int size_m = data->simfields[0]->sphere2DDataConfig->spectral_modes_m_max;

		for (int i = 0; i < N; i++)
			if (!found_nan)
			{
				for (int m = 0; m < size_m; ++m)
				{
					if (!found_nan)
					{
						for (int n = m; n < size_n; ++n)
							if ( std::isnan(data->simfields[i]->spectral_get_(n, m).real()) ||
								std::isnan(data->simfields[i]->spectral_get_(n, m).imag()) )
							{
								found_nan = true;
								break;
							}
					}
				}
			}

		return found_nan;
	}


	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)
	override
	{
#if SWEET_PARAREAL
		SWEET_ASSERT(data->time == i_data.get_pointer_to_data_Sphere2DData_Spectral()->time);
#endif
		SWEET_ASSERT(data->nb_fields == i_data.get_pointer_to_data_Sphere2DData_Spectral()->nb_fields);
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) += *(i_data.get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]);

		return *this;
	}


	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)
	override
	{
#if SWEET_PARAREAL
		SWEET_ASSERT(data->time == i_data.get_pointer_to_data_Sphere2DData_Spectral()->time);
#endif
		SWEET_ASSERT(data->nb_fields == i_data.get_pointer_to_data_Sphere2DData_Spectral()->nb_fields);
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) -= *(i_data.get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]);

		return *this;
	}

	Parareal_GenericData& operator*=(const double v)
	override
	{

		for (int i = 0; i < N; i++)
			*(data->simfields[i]) *= v;

		return *this;
	}


	void grid_print()
	override
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			data->simfields[i]->toGrid().print();
		}
	}

	void spectral_print()
	override
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			data->simfields[i]->spectral_print();
		}
	}

	void dataArrays_2_GenericData_Sphere2DData_Spectral(
								sweet::Data::Sphere2D::DataSpectral &phi,
								sweet::Data::Sphere2D::DataSpectral &vrt,
								sweet::Data::Sphere2D::DataSpectral &div
							) override
	{
		*(data->simfields[0]) = phi;
		*(data->simfields[1]) = vrt;
		*(data->simfields[2]) = div;
	}

	void GenericData_Sphere2DData_Spectral_to_dataArrays(
								sweet::Data::Sphere2D::DataSpectral &phi,
								sweet::Data::Sphere2D::DataSpectral &vrt,
								sweet::Data::Sphere2D::DataSpectral &div
							) override
	{
		phi = *(data->simfields[0]);
		vrt = *(data->simfields[1]);
		div = *(data->simfields[2]);
	}

	void restrict(const Parareal_GenericData& i_data)
	override
	{
		for (int i = 0; i < N; i++)
			*data->simfields[i] = data->simfields[i]->restrict( *(i_data.get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]) );
	}

	void pad_zeros(const Parareal_GenericData& i_data)
	override
	{
		for (int i = 0; i < N; i++)
			*data->simfields[i] = data->simfields[i]->pad_zeros( *(i_data.get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]) );
	}


};

}}

#endif
