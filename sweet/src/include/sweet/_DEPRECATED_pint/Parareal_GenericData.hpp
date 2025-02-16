/*
 * Parareal_GenericData.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_GENERICDATA_HPP

/**
 * Generic Parareal Data class allowing a common parareal interface for any kind of data (cart2d, sphere2D);
 * This class may be inherited and specialized to each data type
 */

#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
#include <sweet/Data/Cart2D/Cart2D.hpp>
#elif SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#endif

namespace sweet {
namespace DEPRECATED_pint {



class Parareal_GenericData
{

public:
	template <class t_dataType>
	class DataContainer
	{
	public:
		int level;
		double time;
		int nb_fields;

	public:
		t_dataType* simfields;

	public:

		DataContainer()
		{
		};

		DataContainer(int i_nb_fields, double i_time, int i_level = 0) :
			level(i_level),
			time(i_time),
			nb_fields(i_nb_fields)
		{
		};

		DataContainer(DataContainer &i_data) :
			level(i_data.level),
			time(i_data.time),
			nb_fields(i_data.nb_fields)
		{
		};

		DataContainer operator=(const DataContainer &i_data)
		{
			level = i_data.level;
			time = i_data.time;
			nb_fields = i_data.nb_fields;
			return *this;
		};


		virtual ~DataContainer()
		{
		}

		void set_time(double i_time)
		{
			time = i_time;
		}

	};

public:

#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
	sweet::Data::Cart2D::Config* cart2DDataConfig = nullptr;
#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
	sweet::Data::Sphere2D::Config* sphere2DDataConfig = nullptr;
#endif

public:
	// different interface functions to avoid template in Parareal_GenericData
	// these interfaces are overridden in the respective child classes
#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	virtual DataContainer<double>* get_pointer_to_data_Scalar() const
	{
		SWEETErrorFatal("This interface function should not be called");
		DataContainer<double>* dummy = nullptr;
		return dummy;
	};

	virtual void dataArrays_2_GenericData_Scalar(
							double &u
							)
	{
	};

	virtual void GenericData_Scalar_2_dataArrays(
							double &u
							)
	{
	};


#elif SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
	virtual DataContainer<sweet::Data::Cart2D::DataSpectral*>* get_pointer_to_data_Cart2DData_Spectral() const
	{
		SWEETErrorFatal("This interface function should not be called");
		DataContainer<sweet::Data::Cart2D::DataSpectral*>* dummy = nullptr;
		return dummy;
	};

	virtual void dataArrays_2_GenericData_Cart2DData_Spectral(
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
								sweet::Data::Cart2D::DataSpectral &h,
	#endif
								sweet::Data::Cart2D::DataSpectral &u,
								sweet::Data::Cart2D::DataSpectral &v
							)
	{
	};

	virtual void GenericData_Cart2DData_Spectral_2_dataArrays(
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
								sweet::Data::Cart2D::DataSpectral &h,
	#endif
								sweet::Data::Cart2D::DataSpectral &u,
								sweet::Data::Cart2D::DataSpectral &v
							)
	{
	};


#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
	virtual DataContainer<sweet::Data::Sphere2D::DataSpectral*>* get_pointer_to_data_Sphere2DData_Spectral() const
	{
		SWEETErrorFatal("This interface function should not be called");
		DataContainer<sweet::Data::Sphere2D::DataSpectral*>* dummy = nullptr;
		return dummy;
	};

	virtual void dataArrays_2_GenericData_Sphere2DData_Spectral(
								sweet::Data::Sphere2D::DataSpectral &phi,
								sweet::Data::Sphere2D::DataSpectral &vrt,
								sweet::Data::Sphere2D::DataSpectral &div
							)
	{
	};

	virtual void GenericData_Sphere2DData_Spectral_to_dataArrays(
								sweet::Data::Sphere2D::DataSpectral &phi,
								sweet::Data::Sphere2D::DataSpectral &vrt,
								sweet::Data::Sphere2D::DataSpectral &div
							)
	{
	};


#endif


#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
	void setup_data_config(sweet::Data::Cart2D::Config* i_cart2DDataConfig)
	{
		cart2DDataConfig = i_cart2DDataConfig;
	};

#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
	void setup_data_config(sweet::Data::Sphere2D::Config* i_sphere2DDataConfig)
	{
		sphere2DDataConfig = i_sphere2DDataConfig;
	};
#endif

public:

	Parareal_GenericData()
	{
	};

	Parareal_GenericData(double i_time, int i_level = 0)
	{
	};

	Parareal_GenericData(Parareal_GenericData &i_data)
	{
	};

	virtual Parareal_GenericData& operator=(const Parareal_GenericData &i_data) = 0;


	//Parareal_GenericData(Parareal_GenericData &&i_data){
	//};

	virtual ~Parareal_GenericData()
	{
	}


	virtual void set_time(double i_time)=0;

	virtual void allocate_data()=0;
	
	virtual void free_data() = 0;

///#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	virtual std::size_t size() = 0;
	////virtual void serialize(void *data) = 0;
	////virtual void deserialize(void *data) = 0;
	#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	virtual void serialize(double *data) = 0;
	virtual void deserialize(double *data) = 0;
	#else
	virtual void serialize(std::complex<double> *data) = 0;
	virtual void deserialize(std::complex<double> *data) = 0;
	#endif
#endif

	virtual double spectral_reduce_maxAbs()=0;
	virtual double spectral_reduce_maxAbs(std::size_t rnorm)=0;
	virtual double grid_reduce_maxAbs()=0;
	virtual double grid_reduce_norm1()=0;
	virtual double grid_reduce_norm2()=0;

	virtual bool check_for_nan()=0;

////	virtual Parareal_GenericData operator+(const Parareal_GenericData &i_data) = 0;  // --> not possible because it is an abstract class
	virtual Parareal_GenericData& operator+=(const Parareal_GenericData &i_data) = 0;

////	virtual Parareal_GenericData operator-(const Parareal_GenericData &i_data) = 0;
	virtual Parareal_GenericData& operator-=(const Parareal_GenericData &i_data) = 0;

	virtual Parareal_GenericData& operator*=(const double v) = 0;

	virtual void restrict(const Parareal_GenericData& i_data) = 0;
	virtual void pad_zeros(const Parareal_GenericData& i_data) = 0;

	virtual void grid_print() = 0;
	virtual void spectral_print() = 0;

};

}}

#endif
