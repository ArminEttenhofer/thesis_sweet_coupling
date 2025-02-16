#ifndef INCLUDE_SWEET_DATA_GENERICCONTAINER_CASTHELPER_HPP
#define INCLUDE_SWEET_DATA_GENERICCONTAINER_CASTHELPER_HPP


#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/GenericContainer/ConfigBase.hpp>

namespace sweet {
namespace Data {
namespace GenericContainer {


template <typename TData, typename TConfig>
class CastHelper
{
public:
	/*
	 * Casting functions which are used for convenience
	 */
	static inline
	TData& cast(GenericContainer::Base &i_U)
	{
		return static_cast<TData&>(i_U);
	}

	static inline
	const TData& cast(const GenericContainer::Base &i_U)
	{
		return static_cast<const TData&>(i_U);
	}

	const TConfig& cast(
			const GenericContainer::ConfigBase& i_deTermConfig
	)
	{
		return static_cast<const TConfig&>(i_deTermConfig);
	}
};

}}}

#endif
