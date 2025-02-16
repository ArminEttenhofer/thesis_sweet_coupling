#ifndef PROGRAMS_LIBPFASST_INTERFACE_SPHERE2DDATAVARS_HPP
#define PROGRAMS_LIBPFASST_INTERFACE_SPHERE2DDATAVARS_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>

// Class containing the prognostic Sphere2DDataSpectral variables phi_pert, vrt, and div

class Sphere2DDataVars {

public:
	Sphere2DDataVars(
			sweet::Data::Sphere2D::Config *sphere2DDataConfig,
			int i_level
	)
	: prog_phi_pert(sphere2DDataConfig),
	  prog_vrt(sphere2DDataConfig),
	  prog_div(sphere2DDataConfig),
	  flat_data_array(nullptr),
	  flat_data_array_size(0),
	  level(i_level)
	{}

	~Sphere2DDataVars()
	{
		if (flat_data_array != nullptr)
		{
			// release the memory
			sweet::Memory::MemBlockAlloc::free(
					flat_data_array,
					flat_data_array_size*sizeof(double)
			);
		}
	}

	void allocate_flat_data_array(int i_n_elems)
	{
		if (flat_data_array == nullptr)
		{
			// allocate the memory
			flat_data_array_size = i_n_elems;
			flat_data_array      = sweet::Memory::MemBlockAlloc::alloc<double>(
					i_n_elems*sizeof(double)
			);
		}
	}

	// getters for the Sphere2DDataSpectral variables
	const sweet::Data::Sphere2D::DataSpectral& get_phi_pert() const  {return prog_phi_pert;}
	sweet::Data::Sphere2D::DataSpectral&       get_phi_pert()        {return prog_phi_pert;}
	const sweet::Data::Sphere2D::DataSpectral& get_vrt() const {return prog_vrt;}
	sweet::Data::Sphere2D::DataSpectral&       get_vrt()       {return prog_vrt;}
	const sweet::Data::Sphere2D::DataSpectral& get_div() const  {return prog_div;}
	sweet::Data::Sphere2D::DataSpectral&       get_div()        {return prog_div;}

	// getters for the flat data array
	double*&         get_flat_data_array()            {return flat_data_array;}
	const int&       get_flat_data_array_size() const {return flat_data_array_size;}
	int&             get_flat_data_array_size()       {return flat_data_array_size;}

	// getters for the level
	const int&       get_level() const {return level;};

protected:

	sweet::Data::Sphere2D::DataSpectral prog_phi_pert;
	sweet::Data::Sphere2D::DataSpectral prog_vrt;
	sweet::Data::Sphere2D::DataSpectral prog_div;

	// flat data array vector (currently used to pack and unpack)
	double *flat_data_array;
	int     flat_data_array_size;

	// pfasst level
	const int level;

	// default constructor, copy constructor, and operator= are disabled
	Sphere2DDataVars();
	Sphere2DDataVars(const Sphere2DDataVars&);
	Sphere2DDataVars& operator=(const Sphere2DDataVars&);
};

#endif
