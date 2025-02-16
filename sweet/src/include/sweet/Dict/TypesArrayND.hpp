/*
 *
 * ArrayND.hpp
 *  Created on: Feb 13, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DICT_TYPESARRAYND_HPP
#define INCLUDE_SWEET_DICT_TYPESARRAYND_HPP

#include <sweet/Error/Fatal.hpp>
#include <array>
#include <vector>
#include <ostream>


namespace sweet {
namespace Dict {

/*!
 * SWEET's 1D/2D/3D array data class for dictionaries.
 *
 * This is not intended for HPC, but just to have
 * some D dimensional container for arbitrary types.
 */
template <int D, typename T>
class TypesArrayND
{
	std::array<int,D> _shape;
	std::array<int,D> _offset;
	std::size_t _size;
	std::vector<T> _data;

public:
	TypesArrayND()
	{
		setZero();
	}


	TypesArrayND(const std::array<int,D> &i_shape)
	{
		setZero();
		setup(i_shape);
	}


	TypesArrayND(const std::array<int,D> &i_shape, const T i_data[])
	{
		setZero();
		setup(i_shape);
		operator=(i_data);
	}

	TypesArrayND(int i_shape0, int i_shape1 = -1, int i_shape2 = -1)
	{
		setZero();
		resize(i_shape0, i_shape1, i_shape2);
	}

	void setZero()
	{
		_size = 0;

		for (int i = 0; i < D; i++)
		{
			_shape[i] = 0;
			_offset[i] = 0;
		}
	}

	void setup(const std::array<int,D> &i_shape)
	{
		if (D < 1 || D > 3)
			SWEETErrorFatal("Only 1D, 2D or 3D are supported!");

		_shape = i_shape;

		_size = 1;
		for (int i = 0; i < D; i++)
			_size *= _shape[i];

		_data.resize(_size);
	}


	void resize(int i_shape0, int i_shape1 = -1, int i_shape2 = -1)
	{
		_shape[0] = i_shape0;

		if (D >= 2)
			_shape[1] = i_shape1;

		if (D >= 3)
			_shape[2] = i_shape2;

		if (D > 3)
			SWEETErrorFatal("Only 1D, 2D or 3D are supported!");

		_size = 1;
		for (int i = 0; i < D; i++)
			_size *= _shape[i];

		_data.resize(_size);
	}


	void setOffset(int i_offset0, int i_offset1 = -1, int i_offset2 = -1)
	{
		_offset[0] = i_offset0;

		if (D >= 2)
			_offset[1] = i_offset1;

		if (D >= 3)
			_offset[2] = i_offset2;

		if (D > 3)
			SWEETErrorFatal("Only 1D, 2D or 3D are supported!");
	}


	void setup(const std::array<int,D> &i_shape, T i_data[])
	{
		setup(i_shape);
		operator=(i_data);
	}

	std::size_t size()	const
	{
		return _size;
	}

	T *data()
	{
		return _data.data();
	}

	const T *data() const
	{
		return _data.data();
	}


	std::vector<T>& dataVector()
	{
		return _data;
	}

	const std::vector<T>& dataVector() const
	{
		return _data;
	}

	const std::array<int,D> &shape()	const
	{
		return _shape;
	}


	inline
	void set(int i0, const T &i_value)
	{
		if (D != 1)
			SWEETErrorFatal("Only for 1D");

		_data[i0+_offset[0]] = i_value;
	}

	inline
	void set(int i0, int i1, const T &i_value)
	{
		if (D != 2)
			SWEETErrorFatal("Only for 2D");

		_data[(i0+_offset[0])*_shape[1] + (i1+_offset[1])] = i_value;
	}

	inline
	void set3(int i0, int i1, int i2, const T &i_value)
	{
		if (D != 3)
			SWEETErrorFatal("Only for 3D");

		_data[(i0+_offset[0])*_shape[1]*_shape[2] + (i1+_offset[1])*_shape[2] + (i2+_offset[2])] = i_value;
	}

	/*
	 * Special getter which is just constant
	 *
	 * This is required if called by other constant functions from this class.
	 */
	const T& getConst(int i0, int i1=-1, int i2=-1) const
	{
		if (D == 1)
		{
			return _data[(i0+_offset[0])];
		}
		else if (D == 2)
		{
			return _data[(i0+_offset[0])*_shape[1] + (i1+_offset[1])];
		}
		else if (D == 3)
		{
			return _data[(i0+_offset[0])*_shape[1]*_shape[2] + (i1+_offset[1])*_shape[2] + (i2+_offset[2])];
		}

		SWEETErrorFatal("Not supported!");
		return _data[0];	// Dummy return to avoid compiler warnings
	}

private:
	const T& _getConstNoOffset(int i0, int i1=-1, int i2=-1) const
	{
		if (D == 1)
		{
			return _data[i0];
		}
		else if (D == 2)
		{
			return _data[i0*_shape[1] + i1];
		}
		else if (D == 3)
		{
			return _data[i0*_shape[1]*_shape[2] + i1*_shape[2] + i2];
		}

		SWEETErrorFatal("Not supported!");
		return _data[0];	// Dummy return to avoid compiler warnings
	}

public:
	T& get(int i0, int i1=-1, int i2=-1)
	{
		SWEET_ASSERT(i0+_offset[0] >= 0);

		if (D >= 2)
		{
			SWEET_ASSERT(i1+_offset[1] >= 0);

			if (D >= 3)
			{
				SWEET_ASSERT(i2+_offset[2] >= 0);
			}
		}

		if (D == 1)
		{
			return _data[i0+_offset[0]];
		}
		else if (D == 2)
		{
			return _data[(i0+_offset[0])*_shape[1] + (i1+_offset[1])];
		}
		else if (D == 3)
		{
			return _data[(i0+_offset[0])*_shape[1]*_shape[2] + (i1+_offset[1])*_shape[2] + (i2+_offset[2])];
		}

		SWEETErrorFatal("Not supported!");
		return _data[0];	// Dummy return to avoid compiler warnings
	}


	/*
	 * This is just a convenience handler to use allow a very compact
	 * access to this array rather than using .get(...)
	 *
	 * The array operator[] is no option for us since this only supports one
	 * argument in C++11
	 */
	inline
	const T& operator()(int i0, int i1=-1, int i2=-1)	const
	{
		return getConst(i0, i1, i2);
	}
	inline
	T& operator()(int i0, int i1=-1, int i2=-1)
	{
		return get(i0, i1, i2);
	}

	inline
	TypesArrayND<D,T>& operator=(const T *i_values_flat)
	{
		for (int i = 0; i < D; i++)
			if (_shape[i] == 0)
				SWEETErrorFatal("Shape is 0, you need to resize array before assigning raw data!");

		// We simply hope that the data is properly allocated
		for (std::size_t i = 0; i < _size; i++)
			_data[i] = i_values_flat[i];

		return *this;
	}

	inline
	bool operator==(const TypesArrayND<D,T> &i_a)	const
	{
		for (std::size_t d = 0; d < D; d++)
			if (_shape[d] != i_a._shape[d])
				return false;

		for (std::size_t i = 0; i < _size; i++)
			if (_data[i] != i_a._data[i])
				return false;

		return true;
	}

	inline
	bool operator!=(const TypesArrayND<D,T> &i_a)	const
	{
		return !operator==(i_a);
	}

	inline
	TypesArrayND<D,T>& operator=(const TypesArrayND<D,T> &a)
	{
		setup(a._shape);

		for (std::size_t i = 0; i < _size; i++)
			_data[i] = a._data[i];

		return *this;
	}

private:
	void _printOffset(std::ostream& os)	const
	{
		os << "offset=(";
		for (int i = 0; i < D; i++)
		{
			os << _offset[i];
			if (i != D-1)
				os << ", ";
		}
		os << ")";
	}

public:
	friend
	std::ostream& operator<<(std::ostream& os, const TypesArrayND<D,T> &a)
	{
		if (a._size == 0)
		{
			os << "(empty)";
			return os;
		}

		if (D == 1)
		{
			os << "[";
			for (int i0 = 0; i0 < a._shape[0]; i0++)
			{
				os << a._getConstNoOffset(i0);

				if (i0 != a._shape[0]-1)
					os << ",\t";
			}
			os << "] ";
			a._printOffset(os);
			os << std::endl;
		}
		else if (D == 2)
		{
			os << "[" << std::endl;
			for (int i0 = 0; i0 < a._shape[0]; i0++)
			{
				os << "\t[";
				for (int i1 = 0; i1 < a._shape[1]; i1++)
				{
					os << a._getConstNoOffset(i0, i1);

					if (i1 != a._shape[1]-1)
						os << ",\t";
				}
				os << "]";
				os << std::endl;
			}
			os << "] ";
			a._printOffset(os);
		}
		else if (D == 3)
		{
			os << "[" << std::endl;
			for (int i0 = 0; i0 < a._shape[0]; i0++)
			{
				os << "\t[" << std::endl;
				for (int i1 = 0; i1 < a._shape[1]; i1++)
				{
					os << "\t\t[";
					for (int i2 = 0; i2 < a._shape[2]; i2++)
					{
						os << a._getConstNoOffset(i0, i1, i2);

						if (i2 != a._shape[2]-1)
							os << ",\t";
					}
					os << "]";
					os << std::endl;
				}
				os << "\t]";
				os << std::endl;
			}
			os << "] ";
			a._printOffset(os);
		}
		else
		{
			SWEETErrorFatal("Not supported!");
		}
		return os;
	}
};

}}

#endif
