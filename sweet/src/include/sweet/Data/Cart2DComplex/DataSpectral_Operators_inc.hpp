#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATASPECTRAL_OPERATORS_INC_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATASPECTRAL_OPERATORS_INC_HPP


inline
static
Cart2DComplex::DataSpectral operator*(
		const double i_value,
		const Cart2DComplex::DataSpectral &i_array_data
)
{
	return i_array_data*i_value;
}


inline
static
Cart2DComplex::DataSpectral operator*(
		const std::complex<double> &i_value,
		const Cart2DComplex::DataSpectral &i_array_data
)
{
	return i_array_data*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */

inline
static
Cart2DComplex::DataSpectral operator+(
		const double i_value,
		const Cart2DComplex::DataSpectral &i_array_data
)
{
	return ((Cart2DComplex::DataSpectral&)i_array_data)+i_value;
}

inline
static
Cart2DComplex::DataSpectral operator+(
		const std::complex<double> &i_value,
		const Cart2DComplex::DataSpectral &i_array_data
)
{
	return i_array_data+i_value;
}

inline
static
Cart2DComplex::DataSpectral operator-(
		const std::complex<double> &i_value,
		const Cart2DComplex::DataSpectral &i_array_data
)
{
	Cart2DComplex::DataSpectral out_cart2d_data(i_array_data.cart2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < i_array_data.cart2DDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_cart2d_data.spectral_space_data[idx] = -i_array_data.spectral_space_data[idx];

	out_cart2d_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_cart2d_data;

}

#endif
