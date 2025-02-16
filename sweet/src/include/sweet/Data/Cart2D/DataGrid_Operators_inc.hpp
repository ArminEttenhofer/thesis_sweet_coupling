#ifndef INCLUDE_SWEET_DATA_CART2D_DATAGRID_OPERATORS_INC_HPP
#define INCLUDE_SWEET_DATA_CART2D_DATAGRID_OPERATORS_INC_HPP

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
Cart2D::DataGrid operator*(
		const double i_value,
		const Cart2D::DataGrid &i_array_data
)
{
	return i_array_data*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 */
inline
static
Cart2D::DataGrid operator-(
		const double i_value,
		const Cart2D::DataGrid &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}


/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 */
inline
static
Cart2D::DataGrid operator+(
		const double i_value,
		const Cart2D::DataGrid &i_array_data
)
{
	return i_array_data+i_value;
}

#endif
