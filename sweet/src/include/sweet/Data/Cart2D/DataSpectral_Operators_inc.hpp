#ifndef INCLUDE_SWEET_DATA_CART2D_DATASPECTRAL_OPERATORS_INC_HPP
#define INCLUDE_SWEET_DATA_CART2D_DATASPECTRAL_OPERATORS_INC_HPP

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
Cart2D::DataSpectral operator*(
		double i_value,
		const Cart2D::DataSpectral &i_array_data
)
{
	return ((Cart2D::DataSpectral&)i_array_data)*i_value;
}



/*!
 * operator to support operations such as:
 *
 * 1.5 + arrayData
 *
 */
inline
static
Cart2D::DataSpectral operator+(
		double i_value,
		const Cart2D::DataSpectral &i_array_data
)
{
	return ((Cart2D::DataSpectral&)i_array_data)+i_value;
}



/*!
 * operator to support operations such as:
 *
 * 1.5 - arrayData
 *
 */
inline
static
Cart2D::DataSpectral operator-(
		double i_value,
		const Cart2D::DataSpectral &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}


/*!
 * operator to support operations such as:
 *
 * array_data_physical - arrayData
 */
inline
static
Cart2D::DataSpectral operator-(
		const Cart2D::DataGrid &i_cart2d_data_physical,
		const Cart2D::DataSpectral &i_array_data
)
{
	return - (i_array_data - i_cart2d_data_physical);
}

#endif
