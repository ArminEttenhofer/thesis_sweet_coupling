/*
 * Dict.hpp
 *
 *  Created on: Feb 13, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DICT__DICTELEMENT_HPP
#define INCLUDE_SWEET_DICT__DICTELEMENT_HPP

#include <sweet/Dict/_DictFileRead.hpp>
#include <sweet/Dict/_DictFileWrite.hpp>
#include <sweet/Dict/_TypesBasic.hpp>
#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/Error/Fatal.hpp>
#include <string>
#include <vector>
#include <complex>
#include <fstream>


namespace sweet {
namespace Dict {


class _DictElementTypeEnum
{
public:
	enum Enum
	{
		SWEET_FILE_DICT_NONE = 0,

		SWEET_FILE_DICT_STRING = 100,
		SWEET_FILE_DICT_INT64 = 230,
		SWEET_FILE_DICT_FLOAT64 = 240,
		SWEET_FILE_DICT_COMPLEX128 = 250,
		SWEET_FILE_DICT_ARRAY_1D_FLOAT64 = 401,
		SWEET_FILE_DICT_ARRAY_2D_FLOAT64 = 402,
		SWEET_FILE_DICT_ARRAY_3D_FLOAT64 = 403,
		SWEET_FILE_DICT_ARRAY_1D_COMPLEX128 = 501,
		SWEET_FILE_DICT_ARRAY_2D_COMPLEX128 = 502,
		SWEET_FILE_DICT_ARRAY_3D_COMPLEX128 = 503,

		SWEET_FILE_DICT_DEFAULT_TYPE,	// special default type
	};
};



class _DictElement : private _TypesBasic
{
private:
	std::string _key;

public:
	const std::string &getKey() const
	{
		return _key;
	}

private:
	_DictElementTypeEnum::Enum _typeId;

public:
	_DictElementTypeEnum::Enum getTypeID()	const
	{
		return _typeId;
	}

public:
	const std::string getTypeIDString() const
	{

		switch(_typeId)
		{
			default:
				SWEETErrorFatal("Unknown type");
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_STRING:
				return "string";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_INT64:
				return "int64";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64:
				return "float64";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128:
				return "complex128";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
				return "array1dFloat64";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
				return "array2dFloat64";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
				return "array3dFloat64";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
				return "array1dComplex128";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
				return "array2dComplex128";
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
				return "array3dComplex128";
				break;
		}

		return "";
	}

public:
	_DictElement()	:
		_typeId(_DictElementTypeEnum::SWEET_FILE_DICT_NONE)
	{
	}

	template <typename T>
	_DictElement(
			const std::string &i_key,
			const T &i_value
	)
	{
		set(i_key, i_value);
	}

	/*
	 * GETTERS and SETTERS
	 */
	template <typename T>
	void get(T &o_value)	const;

	template <typename T>
	void set(const std::string &i_key, const T &i_value);


	/*!
	 * Element of type std::string
	 */
	std::string _valueStr;

	void get(std::string &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_STRING)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueStr;
	}

	void set(const std::string &i_key, const std::string &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_STRING;
		_key = i_key;
		_valueStr = i_value;
	}

	/*
	 * Special handlers which automatically convert values
	 */
	void set(const std::string &i_key, const char *i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_STRING;
		_key = i_key;
		_valueStr = i_value;
	}



	/*
	 * Element of type int64
	 */

	int64 _valueScalarInt64;

	void get(int64 &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_INT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueScalarInt64;
	}

	void set(const std::string &i_key, const int64 &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_INT64;
		_key = i_key;
		_valueScalarInt64 = i_value;
	}

	void get(bool &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_INT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = (bool)_valueScalarInt64;
	}

	void get(int &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_INT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = (int)_valueScalarInt64;
	}

	void set(const std::string &i_key, const int &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_INT64;
		_key = i_key;
		_valueScalarInt64 = (int64)i_value;
	}



	/*
	 * Element of type float64
	 */

	double _valueScalarFloat64;

	void get(float64 &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueScalarFloat64;
	}

	void set(const std::string &i_key, const float64 &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64;
		_key = i_key;
		_valueScalarFloat64 = i_value;
	}


	/*
	 * Element of type complex128
	 */

	complex128 _valueScalarComplex128;

	void get(complex128 &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueScalarComplex128;
	}

	void set(const std::string &i_key, const complex128 &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128;
		_key = i_key;
		_valueScalarComplex128 = i_value;
	}



	/*
	 * Element of type TypesArrayND<1,float64>
	 */

	TypesArrayND<1,double> _valueArray1dFloat64;

	void get(TypesArrayND<1,float64> &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueArray1dFloat64;
	}

	void set(const std::string &i_key, const TypesArrayND<1,float64> &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64;
		_key = i_key;
		_valueArray1dFloat64 = i_value;
	}



	/*
	 * Element of type TypesArrayND<2,float64>
	 */

	TypesArrayND<2,double> _valueArray2dFloat64;

	void get(TypesArrayND<2,float64> &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueArray2dFloat64;
	}

	void set(const std::string &i_key, const TypesArrayND<2,float64> &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64;
		_key = i_key;
		_valueArray2dFloat64 = i_value;
	}



	/*
	 * Element of type TypesArrayND<3,float64>
	 */

	TypesArrayND<3,double> _valueArray3dFloat64;

	void get(TypesArrayND<3,float64> &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueArray3dFloat64;
	}

	void set(const std::string &i_key, const TypesArrayND<3,float64> &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64;
		_key = i_key;
		_valueArray3dFloat64 = i_value;
	}



	/*
	 * Element of type TypesArrayND<1,complex128>
	 */

	TypesArrayND<1,complex128> _valueArray1dComplex128;

	void get(TypesArrayND<1,complex128> &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueArray1dComplex128;
	}

	void set(const std::string &i_key, const TypesArrayND<1,complex128> &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128;
		_key = i_key;
		_valueArray1dComplex128 = i_value;
	}



	/*
	 * Element of type TypesArrayND<2,complex128>
	 */

	TypesArrayND<2,complex128> _valueArray2dComplex128;

	void get(TypesArrayND<2,complex128> &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueArray2dComplex128;
	}

	void set(const std::string &i_key, const TypesArrayND<2,complex128> &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128;
		_key = i_key;
		_valueArray2dComplex128 = i_value;
	}



	/*
	 * Element of type TypesArrayND<3,complex128>
	 */

	TypesArrayND<3,complex128> _valueArray3dComplex128;

	void get(TypesArrayND<3,complex128> &o_value)	const
	{
		if (_typeId != _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128)
			SWEETErrorFatal("Type mismatch for key '"+_key+"'");

		o_value = _valueArray3dComplex128;
	}

	void set(const std::string &i_key, const TypesArrayND<3,complex128> &i_value)
	{
		_typeId = _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128;
		_key = i_key;
		_valueArray3dComplex128 = i_value;
	}


	/*
	 * Print value
	 */
	std::string getValueAsString()	const
	{
		std::ostringstream ss;

		switch(_typeId)
		{
			default:
				SWEETErrorFatal("Unknown type");
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_STRING:
				ss << _valueStr;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_INT64:
				ss << _valueScalarInt64;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64:
				ss << _valueScalarFloat64;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128:
				ss << _valueScalarComplex128;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
				ss << _valueArray1dFloat64;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
				ss << _valueArray2dFloat64;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
				ss << _valueArray3dFloat64;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
				ss << _valueArray1dComplex128;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
				ss << _valueArray2dComplex128;
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
				ss << _valueArray3dComplex128;
				break;
		}

		return ss.str();
	}

	/*
	 * Print value
	 */
	friend
	std::ostream& operator<<(
			std::ostream& io_os,
			const _DictElement &i_element
	)
	{
		io_os << "'" << i_element._key << "' (" << i_element.getTypeIDString() << ") => '" <<  i_element.getValueAsString() << "'";

		return io_os;
	}


public:
	void fileLoadKeyTypeValue(
			_DictFileRead &io_f
	)
	{
		std::string key = io_f.loadStr0();
		_DictElementTypeEnum::Enum type_id = (_DictElementTypeEnum::Enum)io_f.loadData<int64>();

		switch(type_id)
		{
			case _DictElementTypeEnum::SWEET_FILE_DICT_STRING:
				set(key, io_f.loadStr0());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_INT64:
				set(key, io_f.loadData<int64>());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64:
				set(key, io_f.loadData<float64>());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128:
				set(key, io_f.loadData<complex128>());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
			{
				std::array<int,1> shape;
				shape[0] = io_f.loadData<int64>();
				TypesArrayND<1, float64> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
			{
				std::array<int,2> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				TypesArrayND<2, float64> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
			{
				std::array<int,3> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				shape[2] = io_f.loadData<int64>();
				TypesArrayND<3, float64> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
			{
				std::array<int,1> shape;
				shape[0] = io_f.loadData<int64>();
				TypesArrayND<1, complex128> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
			{
				std::array<int,2> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				TypesArrayND<2, complex128> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
			{
				std::array<int,3> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				shape[2] = io_f.loadData<int64>();
				TypesArrayND<3, complex128> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			default:
			{
				std::ostringstream ss;
				ss << "Unknown type id '" << type_id << "'";
				SWEETErrorFatal(ss.str());
			}
		}

	}


public:
	void fileReadKeyTypeValue(
			_DictFileRead &io_f
	)
	{
		std::string key = io_f.loadStr0();
		_DictElementTypeEnum::Enum type_id = (_DictElementTypeEnum::Enum)io_f.loadData<int64>();

		switch(type_id)
		{
			case _DictElementTypeEnum::SWEET_FILE_DICT_STRING:
				set(key, io_f.loadStr0());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_INT64:
				set(key, io_f.loadData<int64>());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64:
				set(key, io_f.loadData<float64>());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128:
				set(key, io_f.loadData<complex128>());
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
			{
				std::array<int,1> shape;
				shape[0] = io_f.loadData<int64>();
				TypesArrayND<1, float64> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
			{
				std::array<int,2> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				TypesArrayND<2, float64> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
			{
				std::array<int,3> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				shape[2] = io_f.loadData<int64>();
				TypesArrayND<3, float64> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
			{
				std::array<int,1> shape;
				shape[0] = io_f.loadData<int64>();
				TypesArrayND<1, complex128> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
			{
				std::array<int,2> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				TypesArrayND<2, complex128> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
			{
				std::array<int,3> shape;
				shape[0] = io_f.loadData<int64>();
				shape[1] = io_f.loadData<int64>();
				shape[2] = io_f.loadData<int64>();
				TypesArrayND<3, complex128> array(shape);

				io_f.loadArray(array);
				set(key, array);
			}
				break;

			default:
			{
				std::ostringstream ss;
				ss << "Unknown type id '" << type_id << "'";
				SWEETErrorFatal(ss.str());
			}
		}

	}


public:
	void fileWriteKeyTypeValue(
			_DictFileWrite &io_f
	)
	{
		io_f.writeStr0(_key);
		io_f.writeData<int64>(_typeId);

		switch(_typeId)
		{
			case _DictElementTypeEnum::SWEET_FILE_DICT_STRING:
				io_f.writeStr0(_valueStr);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_INT64:
				io_f.writeData(_valueScalarInt64);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64:
				io_f.writeData(_valueScalarFloat64);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128:
				io_f.writeData(_valueScalarComplex128);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
				io_f.writeDictArray(_valueArray1dFloat64);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
				io_f.writeDictArray(_valueArray2dFloat64);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
				io_f.writeDictArray(_valueArray3dFloat64);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
				io_f.writeDictArray(_valueArray1dComplex128);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
				io_f.writeDictArray(_valueArray2dComplex128);
				break;

			case _DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
				io_f.writeDictArray(_valueArray3dComplex128);
				break;

			default:
			{
				std::ostringstream ss;
				ss << "Unknown type id '" << _typeId << "'";
				SWEETErrorFatal(ss.str());
			}
		}

	}
};

}}

#endif
