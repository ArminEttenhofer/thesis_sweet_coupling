#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATAGRID_OPERATORS_INC_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_DATAGRID_OPERATORS_INC_HPP

/**
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
sweet::Data::Cart2DComplex::DataGrid operator*(
		const double i_value,
		const sweet::Data::Cart2DComplex::DataGrid &i_array_data
)
{
	return ((sweet::Data::Cart2DComplex::DataGrid&)i_array_data)*i_value;
}


inline
static
sweet::Data::Cart2DComplex::DataGrid operator*(
		const std::complex<double> &i_value,
		const sweet::Data::Cart2DComplex::DataGrid &i_array_data
)
{
	return ((sweet::Data::Cart2DComplex::DataGrid&)i_array_data)*i_value;
}


#endif
