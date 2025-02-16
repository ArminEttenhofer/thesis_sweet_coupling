/*
 * Parareal_GenericData_Cart2DData_Spectral.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_CART2DDATA_SPECTRAL_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_CART2DDATA_SPECTRAL_HPP

#include <assert.h>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>

namespace sweet {
namespace DEPRECATED_pint {


template <int N>
class Parareal_GenericData_Cart2DData_Spectral :
		public sweet::DEPRECATED_pint::Parareal_GenericData
{
	class DataContainer_Cart2DData_Spectral:
			public sweet::DEPRECATED_pint::Parareal_GenericData::DataContainer<sweet::Data::Cart2D::DataSpectral*>
	{


	public:

		DataContainer_Cart2DData_Spectral(sweet::Data::Cart2D::Config* i_cart2DDataConfig)
		{
			nb_fields = N;
			simfields = new sweet::Data::Cart2D::DataSpectral*[N];
			for (int i = 0; i < N; i++)
				simfields[i] = new sweet::Data::Cart2D::DataSpectral(i_cart2DDataConfig);
		};

		DataContainer_Cart2DData_Spectral(
				sweet::Data::Cart2D::DataSpectral* i_simfields[N]
		)
		{
			nb_fields = N;
			simfields = new sweet::Data::Cart2D::DataSpectral*[N];
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_simfields[i]);
		};

		DataContainer_Cart2DData_Spectral(DataContainer_Cart2DData_Spectral &i_data)
		{
			nb_fields = N;
			level = i_data.level;
			simfields = new sweet::Data::Cart2D::DataSpectral*[N];
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_data.simfields[i]);
		};

		DataContainer_Cart2DData_Spectral& operator=(const DataContainer_Cart2DData_Spectral &i_data)
		{
			nb_fields = N;
			level = i_data.level;
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_data.simfields[i]);
			return *this;
		};


		~DataContainer_Cart2DData_Spectral()
		{
			for (int i = 0; i < N; ++i)
				if (simfields[i])
					delete simfields[i];
			if (simfields)
				delete [] simfields;
		}
	};

public:

	DataContainer<sweet::Data::Cart2D::DataSpectral*>* data = nullptr;

public:
	DataContainer<sweet::Data::Cart2D::DataSpectral*>* get_pointer_to_data_Cart2DData_Spectral() const override
	{
		return data;
	};

public:

	Parareal_GenericData_Cart2DData_Spectral():
		sweet::DEPRECATED_pint::Parareal_GenericData()
	{
		////allocate_data();
	}

	Parareal_GenericData_Cart2DData_Spectral(Parareal_GenericData_Cart2DData_Spectral &i_data)
	{
		*(data) = *(i_data.get_pointer_to_data_Cart2DData_Spectral());
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) = *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]);
		cart2DDataConfig = i_data.cart2DDataConfig;
	};

	Parareal_GenericData_Cart2DData_Spectral& operator=(const Parareal_GenericData &i_data)
	{
		*(data) = *(i_data.get_pointer_to_data_Cart2DData_Spectral());
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) = *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]);
		cart2DDataConfig = i_data.cart2DDataConfig;
		return *this;
	};


	~Parareal_GenericData_Cart2DData_Spectral()
	{
		free_data();
	};


	void set_time(double i_time)
	{
		data->set_time(i_time);
	}

	void allocate_data()
	{
		data = new DataContainer_Cart2DData_Spectral(cart2DDataConfig);
	}
	
	void free_data()
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
	}



	/**
	 * Setup data
	 */
	void setup(
				sweet::Data::Cart2D::DataSpectral* i_simfields[N]
	)
	{
		get_pointer_to_data_Cart2DData_Spectral()->simfields = i_simfields;
	}


////#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	// size in bytes (for MPI)
	// size of each simfield of data
	std::size_t size()
	{
		return N * data->simfields[0]->cart2DDataConfig->spectral_array_data_number_of_elements;
	}

	void serialize(std::complex<double> *i_data)
	{
		int s = data->simfields[0]->cart2DDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&data->simfields[i]->spectral_space_data[0], &data->simfields[i]->spectral_space_data[s], &i_data[i * s]);
	};

	void deserialize(std::complex<double> *i_data)
	{
		int s = data->simfields[0]->cart2DDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&i_data[i * s], &i_data[(i + 1) * s], &data->simfields[i]->spectral_space_data[0]);
	};

#endif

	double spectral_reduce_maxAbs()
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
				data->simfields[k]->spectral_reduce_max_abs());
		return e;
	}

	double spectral_reduce_maxAbs(std::size_t rnorm)
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
				data->simfields[k]->spectral_reduce_max_abs(rnorm));
		return e;
	}

	double grid_reduce_maxAbs()
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
				data->simfields[k]->toGrid().grid_reduce_max_abs());
		return e;
	}

	double grid_reduce_norm1()
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += data->simfields[k]->toGrid().grid_reduce_norm1();
		return e;
	}

	double grid_reduce_norm2()
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
	{
		bool found_nan = false;

		int grid_size_x = data->simfields[0]->cart2DDataConfig->grid_data_size[0];
		int grid_size_y = data->simfields[0]->cart2DDataConfig->grid_data_size[1];

		for (int i = 0; i < N; i++)
		{
			sweet::Data::Cart2D::DataGrid data_phys = data->simfields[i]->toGrid();
			if (!found_nan)
			{
				for (int ix = 0; ix < grid_size_x; ++ix)
				{
					if (!found_nan)
					{
						for (int iy = 0; iy < grid_size_y; ++iy)
							if ( std::isnan(data_phys.grid_get(ix, iy)))
							{
								found_nan = true;
								break;
							}
					}
				}
			}
		}

		return found_nan;
	}


	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)
	{
#if SWEET_PARAREAL
		SWEET_ASSERT(data->time == i_data.get_pointer_to_data_Cart2DData_Spectral()->time);
#endif
		SWEET_ASSERT(data->nb_fields == i_data.get_pointer_to_data_Cart2DData_Spectral()->nb_fields);

		for (int i = 0; i < N; i++)
			*(data->simfields[i]) += *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]);

		return *this;
	}

	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)
	{
#if SWEET_PARAREAL
		SWEET_ASSERT(data->time == i_data.get_pointer_to_data_Cart2DData_Spectral()->time);
#endif
		SWEET_ASSERT(data->nb_fields == i_data.get_pointer_to_data_Cart2DData_Spectral()->nb_fields);

		for (int i = 0; i < N; i++)
			*(data->simfields[i]) -= *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]);

		return *this;
	}

	Parareal_GenericData& operator*=(const double v)
	{

		for (int i = 0; i < N; i++)
			*(data->simfields[i]) *= v;

		return *this;
	}

	void grid_print()
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			data->simfields[i]->toGrid().print();
		}
	}

	void spectral_print()
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			data->simfields[i]->spectral_print();
		}
	}


	void dataArrays_2_GenericData_Cart2DData_Spectral(
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
								sweet::Data::Cart2D::DataSpectral &h,
	#endif
								sweet::Data::Cart2D::DataSpectral &u,
								sweet::Data::Cart2D::DataSpectral &v
							) override
	{
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
		*(data->simfields[0]) = h;
		*(data->simfields[1]) = u;
		*(data->simfields[2]) = v;
	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
		*(data->simfields[0]) = u;
		*(data->simfields[1]) = v;
	#endif
	}

	void GenericData_Cart2DData_Spectral_2_dataArrays(
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
								sweet::Data::Cart2D::DataSpectral &h,
	#endif
								sweet::Data::Cart2D::DataSpectral &u,
								sweet::Data::Cart2D::DataSpectral &v
							) override
	{
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
		h = *(data->simfields[0]);
		u = *(data->simfields[1]);
		v = *(data->simfields[2]);
	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
		u = *(data->simfields[0]);
		v = *(data->simfields[1]);
	#endif
	}


	void restrict(const Parareal_GenericData& i_data)
	{
		for (int i = 0; i < N; i++)
			*data->simfields[i] = data->simfields[i]->restrict( *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]) );
			///data->simfields[i]->restrict( *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]) );
	}

	void pad_zeros(const Parareal_GenericData& i_data)
	{
		for (int i = 0; i < N; i++)
			*data->simfields[i] = data->simfields[i]->pad_zeros( *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]) );
			///data->simfields[i]->pad_zeros( *(i_data.get_pointer_to_data_Cart2DData_Spectral()->simfields[i]) );
	}

};

}}

#endif
