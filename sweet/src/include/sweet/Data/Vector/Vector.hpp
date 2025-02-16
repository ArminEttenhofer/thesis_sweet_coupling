/*
 * ScalarDataArray.hpp
 *
 *  Created on: 28 Jun 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef INCLUDE_SWEET_DATA_VECTOR_VECTOR_HPP
#define INCLUDE_SWEET_DATA_VECTOR_VECTOR_HPP

#include <string>
#include <ostream>
#include <fstream>
#include <limits>
#include <functional>
#include <cmath>
#include <iomanip>

#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>

/*
 * Precompiler helper functions to handle loops in spectral and physical space
 */

namespace sweet {
namespace Data {
namespace Vector {


/*!
 * Storage space for a simple vector of a generic datatype
 */
template <typename T = double>
class Vector
{

public:
	/*!
	 * Number of elements stored in this vector
	 */
	std::size_t numberOfElements;

	/*!
	 * physical space data
	 */
	T *data;

	/*!
	 * allow empty initialization
	 */
public:
	Vector()	:
		numberOfElements(0),
		data(nullptr)
	{
	}



private:
	void _allocate_buffers()
	{
		data = sweet::Memory::MemBlockAlloc::alloc<T>(
				numberOfElements*sizeof(T)
		);
	}

public:
	std::size_t size()	const
	{
		return static_cast<std::size_t>(numberOfElements);
	}

public:
	/*!
	 * copy constructor, used e.g. in
	 * 	ScalarDataArray tmp_h = h;
	 * 	ScalarDataArray tmp_h2(h);
	 *
	 * Duplicate all data
	 */
	Vector(
			const Vector &i_dataArray
	)
	{
		numberOfElements = i_dataArray.numberOfElements;

		_allocate_buffers();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
				data[idx] = i_dataArray.data[idx];
	}



	/*!
	 * default constructor
	 */
public:
	Vector(
			std::size_t i_number_of_elements
	)	: numberOfElements(i_number_of_elements)
	{
		numberOfElements = i_number_of_elements;

		if (numberOfElements == 0)
			return;

		_allocate_buffers();
	}


public:
	/*!
	 * setup the ScalarDataArray in case that the special
	 * empty constructor with int as a parameter was used.
	 *
	 * Calling this setup function should be in general avoided.
	 */
public:
	void setup(
			std::size_t i_numberOfElements
	)
	{
		clear();

		numberOfElements = i_numberOfElements;

		_allocate_buffers();
	}

public:
	void swap(
			Vector &i_vector
	)
	{
		SWEET_ASSERT(numberOfElements == i_vector.numberOfElements);
		SWEET_ASSERT(data != nullptr);
		SWEET_ASSERT(i_vector.data != nullptr);

		std::swap(data, i_vector.data);
	}

public:
	/*!
	 * setup the ScalarDataArray in case that the special
	 * empty constructor with int as a parameter was used.
	 *
	 * Calling this setup function should be in general avoided.
	 */
public:
	void setup_if_required(std::size_t i_number_of_elements)
	{
		if (data != nullptr)
			return;

		numberOfElements = i_number_of_elements;

		_allocate_buffers();
	}


public:
	/*!
	 * setup the ScalarDataArray in case that the special
	 * empty constructor with int as a parameter was used.
	 *
	 * Calling this setup function should be in general avoided.
	 */
public:
	void setup_if_required(const Vector &i_data)
	{
		if (data != nullptr)
			return;

		numberOfElements = i_data.numberOfElements;

		_allocate_buffers();
	}



	void clear()
	{
		if (data == nullptr)
			return;

		sweet::Memory::MemBlockAlloc::free(data, numberOfElements*sizeof(T));
		data = nullptr;
	}


public:
	~Vector()
	{
		clear();
	}


	inline
	void set(
			std::size_t i,
			T i_value
	)
	{
		data[i] = i_value;
	}



	inline
	void set_all(
			T i_value
	)
	{
		for (std::size_t i = 0; i < numberOfElements; i++)
			data[i] = i_value;
	}


	inline
	void set_zero()
	{
		for (std::size_t i = 0; i < numberOfElements; i++)
			data[i] = 0;
	}


	inline
	T get(
			std::size_t i
	)
	{
		return data[i];
	}


	inline
	void update_lambda_array_indices(
			std::function<void(int,T&)> i_lambda	//!< lambda function to return value for lat/mu
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			i_lambda(idx, data[idx]);
	}


	inline
	T grid_get(
			std::size_t i
	)	const
	{
		return data[i];
	}


	inline
	void grid_set_all(
			T i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] = i_value;
	}


	inline
	void grid_setZero()
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] = 0;
	}

	bool reduce_isAnyNaNorInf()	const
	{
		for (std::size_t i = 0; i < numberOfElements; i++)
		{
			if (
					std::isnan(data[i]) ||
					std::isinf(data[i])
			)
				return true;
		}

		return false;
	}



	/*!
	 * return true, if any value is infinity
	 */
	bool reduce_boolean_all_finite() const
	{
		bool isallfinite = true;

#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(&&:isallfinite)
#endif
		for (std::size_t i = 0; i < numberOfElements; i++)
			isallfinite = isallfinite && std::isfinite(data[i]);

		return isallfinite;
	}



	/*!
	 * return the maximum of all absolute values
	 */
	T reduce_maxAbs()	const
	{
		T maxabs = -1.0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(maxabs)
		for (std::size_t i = 0; i < numberOfElements; i++)
			maxabs = std::max(maxabs, std::abs(data[i]));

		return maxabs;
	}



	/*!
	 * reduce to root mean square
	 */
	T reduce_rms()
	{
		T sum = 0;
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
		for (std::size_t i = 0; i < numberOfElements; i++)
			sum += data[i]*data[i];

		sum = std::sqrt(sum/(T)(numberOfElements));

		return sum;
	}


	/*!
	 * reduce to root mean square
	 */
	T reduce_rms_quad()
	{
		T sum = 0;
		T c = 0;
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
		for (std::size_t i = 0; i < numberOfElements; i++)
		{
			T value = data[i]*data[i];

			// Use Kahan summation
			T y = value - c;
			T t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		sum = std::sqrt(sum/(T)(numberOfElements));

		return sum;
	}



	/*!
	 * return the maximum of all absolute values
	 */
	T reduce_max()	const
	{
		T maxvalue = -std::numeric_limits<T>::max();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_MAX(maxvalue)
		for (std::size_t i = 0; i < numberOfElements; i++)
			maxvalue = std::max(maxvalue, data[i]);

		return maxvalue;
	}


	/*!
	 * return the maximum of all absolute values
	 */
	T reduce_min()	const
	{
		T minvalue = std::numeric_limits<T>::max();
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(min:minvalue)
#endif
		for (std::size_t i = 0; i < numberOfElements; i++)
			minvalue = std::min(minvalue, data[i]);

		return minvalue;
	}


	/*!
	 * return the maximum of all absolute values
	 */
	T reduce_sum()	const
	{
		T sum = 0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
		for (std::size_t i = 0; i < numberOfElements; i++)
			sum += data[i];

		return sum;
	}


	/*!
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	T reduce_sum_quad()	const
	{
		T sum = 0;
		T c = 0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
		for (std::size_t i = 0; i < numberOfElements; i++)
		{
			T value = data[i];

			// Use Kahan summation
			T y = value - c;
			T t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return sum;
	}

	/*!
	 * return the maximum of all absolute values
	 */
	T reduce_norm1()	const
	{
		T sum = 0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
		for (std::size_t i = 0; i < numberOfElements; i++)
			sum += std::abs(data[i]);


		return sum;
	}

	/*!
	 * return the sum of the absolute values.
	 */
	T reduce_norm1_quad()	const
	{
		T sum = 0;
		T c = 0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
		for (std::size_t i = 0; i < numberOfElements; i++)
		{

			T value = std::abs(data[i]);
			// Use Kahan summation
			T y = value - c;
			T t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return sum;
	}


	/*!
	 * return the sqrt of the sum of the squared values
	 */
	T reduce_norm2()	const
	{
		T sum = 0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM(sum)
		for (std::size_t i = 0; i < numberOfElements; i++)
			sum += data[i]*data[i];


		return std::sqrt(sum);
	}


	/*!
	 * return the sqrt of the sum of the squared values, use quad precision for reduction
	 */
	T reduce_norm2_quad()	const
	{
		T sum = 0.0;
		T c = 0.0;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_REDUCE_SUM2(sum,c)
		for (std::size_t i = 0; i < numberOfElements; i++)
		{
			T value = data[i]*data[i];

			// Use Kahan summation
			T y = value - c;
			T t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return std::sqrt(sum);
	}


public:
	/*!
	 * array operator
	 */
	T operator[](std::size_t i_index)	const
	{
		return data[i_index];
	}



public:
	/*!
	 * array operator
	 */
	T& operator[](std::size_t i_index)
	{
		return data[i_index];
	}



public:
	/*!
	 * assignment operator
	 */
	Vector &operator=(T i_value)
	{
		grid_set_all(i_value);

		return *this;
	}


public:
	/*!
	 * assignment operator
	 */
	Vector &operator=(int i_value)
	{
		grid_set_all(i_value);

		return *this;
	}


public:
	/*!
	 * assignment operator
	 */
	Vector &operator=(
			const Vector &i_data
	)
	{
		if (numberOfElements != i_data.numberOfElements)
			setup(i_data.numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] = i_data.data[idx];

		return *this;
	}

public:
	/*!
	 * assignment reference operator
	 */
	Vector &operator=(
			Vector &&i_data
	)
	{
		if (numberOfElements != i_data.numberOfElements)
			setup(i_data.numberOfElements);

		std::swap(data, i_data.data);

		return *this;
	}


	/*!
	 * Compute element-wise addition
	 */
	inline
	Vector operator+(
			const Vector &i_array_data
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx] + i_array_data.data[idx];

		return out;
	}

	/*!
	 * Compute element-wise square root
	 */
	inline
	Vector sqrt()	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = std::sqrt(data[idx]);

		return out;
	}


	/*!
	 * Compute element-wise inverse square root
	 */
	inline
	Vector inv_sqrt()	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = 1.0/std::sqrt(data[idx]);

		return out;
	}



	/*!
	 * Compute element-wise addition
	 */
	inline
	Vector operator+(
			const T i_value
	)	const
	{
		Vector out = *this;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]+i_value;

		return out;
	}



	/*!
	 * Compute element-wise addition
	 */
	inline
	Vector& operator+=(
			const Vector &i_array_data
	)
	{

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] += i_array_data.data[idx];

		return *this;
	}



	/*!
	 * Compute element-wise addition
	 */
	inline
	Vector& operator+=(
			const T i_value
	)
	{

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] += i_value;

		return *this;
	}



	/*!
	 * Compute multiplication with scalar
	 */
	inline
	Vector& operator*=(
			const T i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] *= i_value;

		return *this;
	}



	/*!
	 * Compute division with scalar
	 */
	inline
	Vector& operator/=(
			const T i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] /= i_value;

		return *this;
	}



	/*!
	 * Compute element-wise subtraction
	 */
	inline
	Vector& operator-=(
			const Vector &i_array_data
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
				data[idx] -= i_array_data.data[idx];

		return *this;
	}


	/*!
	 * Compute element-wise subtraction
	 */
	inline
	Vector operator-(
			const Vector &i_array_data
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]-i_array_data.data[idx];

		return out;
	}



	/*!
	 * Compute element-wise subtraction
	 */
	inline
	Vector operator-(
			const T i_value
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]-i_value;

		return out;
	}



	/*!
	 * Compute element-wise subtraction
	 */
	inline
	Vector valueMinusThis(
			const T i_value
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = i_value - data[idx];

		return out;
	}


	/*!
	 * Compute element-wise subtraction
	 */
	inline
	Vector valueDivThis(
			const T i_value
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = i_value/data[idx];

		return out;
	}



	/*!
	 * Compute sine
	 */
	inline
	Vector sin()	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = std::sin(data[idx]);

		return out;
	}



	/*!
	 * Compute cosine
	 */
	inline
	Vector cos()	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = std::cos(data[idx]);

		return out;
	}


	/*!
	 * Compute power of two
	 */
	inline
	Vector pow(T i_pow)		const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = std::pow(data[idx], i_pow);

		return out;
	}


	/*!
	 * Compute power of two
	 */
	inline
	Vector pow2()	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]*data[idx];

		return out;
	}

	/*!
	 * Compute power of three
	 */
	inline
	Vector pow3()	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]*data[idx]*data[idx];

		return out;
	}



	/*!
	 * Invert sign
	 */
	inline
	Vector operator-()
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = -data[idx];

		return out;
	}



	/*!
	 * Compute element-wise multiplication
	 */
	inline
	Vector operator*(
			const Vector &i_array_data	//!< this class times i_array_data
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
				out.data[idx] = data[idx]*i_array_data.data[idx];

		return out;
	}


	/*!
	 * Compute element-wise multiplication
	 */
	inline
	Vector& operator*=(
			const Vector &i_array_data	//!< this class times i_array_data
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			data[idx] *= i_array_data.data[idx];

		return *this;
	}



	/*!
	 * Compute element-wise division
	 */
	inline
	Vector operator/(
			const Vector &i_array_data	//!< this class times i_array_data
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]/i_array_data.data[idx];

		return out;
	}



	/*!
	 * Compute multiplication with a scalar
	 */
	inline
	Vector operator*(
			const T i_value
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx]*i_value;

		return out;
	}


	/*!
	 * Compute element-wise division
	 */
	inline
	Vector operator/(
			const T &i_value
	)	const
	{
		Vector out(numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < numberOfElements; idx++)
			out.data[idx] = data[idx] / i_value;

		return out;
	}



	/*!
	 * Print data
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool print(
			int i_precision = 8		//!< number of floating point digits
			)	const
	{
		std::ostream &o_ostream = std::cout;

		o_ostream << std::setprecision(i_precision);

		for (std::size_t i = 0; i < numberOfElements; i++)
			o_ostream << data[i] << "\t";
		std::cout << std::endl;

		return true;
	}



	void file_write_raw(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const
	{
		std::fstream file(i_filename, std::ios::out | std::ios::binary);
		file.write((const char*)data, sizeof(T)*numberOfElements);
	}



	void file_read_raw(
			const std::string &i_filename		//!< Name of file to load data from
	)
	{
		std::fstream file(i_filename, std::ios::in | std::ios::binary);
		file.read((char*)data, sizeof(T)*numberOfElements);
	}


	friend
	inline
	std::ostream& operator<<(
			std::ostream &o_ostream,
			const Vector &i_dataArray
	)
	{
		for (std::size_t i = 0; i < i_dataArray.numberOfElements; i++)
			std::cout << i_dataArray.data[i] << "\t";

		return o_ostream;
	}

};


/*!
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
Vector<double> operator*(
		const double i_value,
		const Vector<double> &i_array_data
)
{
	return ((Vector<double>&)i_array_data)*i_value;
}

/*!
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 *
 * Otherwise, we'd have to write it as arrayData-1.5
 *
 */
inline
static
Vector<double> operator-(
		const double i_value,
		const Vector<double> &i_array_data
)
{
	return i_array_data.valueMinusThis(i_value);
}

/*!
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */
inline
static
Vector<double> operator+(
		const double i_value,
		const Vector<double> &i_array_data
)
{
	return i_array_data+i_value;
}



inline
static
Vector<double> operator/(
		const double i_value,
		const Vector<double> &i_array_data
)
{
	return i_array_data.valueDivThis(i_value);
}


/*
 * Namespace to use for convenient sin/cos/pow/... calls
 */
namespace DataVector_ops
{
	inline
	static
	double pow2(double i_value)
	{
		return i_value*i_value;
	}

	inline
	static
	double pow3(double i_value)
	{
		return i_value*i_value*i_value;
	}

	inline
	static
	Vector<double> sin(
			const Vector<double> &i_array_data
	)
	{
		return i_array_data.sin();
	}

	inline
	static
	Vector<double> cos(
			const Vector<double> &i_array_data
	)
	{
		return i_array_data.cos();
	}

	inline
	static
	Vector<double> pow2(
			const Vector<double> &i_array_data
	)
	{
		return i_array_data.pow2();
	}

	inline
	static
	Vector<double> pow3(
			const Vector<double> &i_array_data
	)
	{
		return i_array_data.pow3();
	}

	inline
	static
	Vector<double> pow(
			const Vector<double> &i_array_data,
			double i_value
	)
	{
		return i_array_data.pow(i_value);
	}

};

}}}

#endif
