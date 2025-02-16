#ifndef INCLUDE_SWEET_DATA_GENERICCONTAINER_CONFIGBASE_HPP
#define INCLUDE_SWEET_DATA_GENERICCONTAINER_CONFIGBASE_HPP


#include <sweet/Data/GenericContainer/Base.hpp>


namespace sweet {
namespace Data {
namespace GenericContainer {



/*!
 * A special class which is forwarded to all
 *
 *  - time stepper instances and
 *  - DE term instances
 *
 * to set up data buffers and other things which are required.
 */
class ConfigBase
{
public:
	ConfigBase()
	{
	}

	virtual ~ConfigBase()
	{
	}

	/**
	 * Return a new instance of a data object related to this time stepper
	 */
	virtual
	Base* getNewDataContainerInstance(int i_id = Base::DATA_SIMULATION) const = 0;
};

}}}

#endif
