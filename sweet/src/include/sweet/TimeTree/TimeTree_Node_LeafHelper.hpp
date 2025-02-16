#ifndef INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_LEAFHELPER_HPP
#define INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_LEAFHELPER_HPP

#include <sweet/TimeTree/TimeTree_Node_Base.hpp>


namespace sweet {
namespace TimeTree {


/*!
 * Helper class for interior tree node
 *
 * It provides default member variables which are often used
 */
template <typename MainInteriorNodeClass>
class TimeTree_Node_LeafHelper	:
		public TimeTree_Node_Base
{
protected:
	double _timestepSize;
	double &_dt = _timestepSize;

	// Number of stages to allocate buffers
	std::vector<sweet::Data::GenericContainer::Base*> _tmpDataContainer;


public:
	TimeTree_Node_LeafHelper():
		_timestepSize(-1)
	{
	}


	/*!
	 * Copy constructor
	 *
	 * Simply copy the raw data over here
	 */
	TimeTree_Node_LeafHelper(
		const TimeTree_Node_LeafHelper &i_src	//!<! Source node to copy from
	)	:
		TimeTree_Node_Base(i_src)
	{
		_timestepSize = i_src._timestepSize;

		_tmpDataContainer.resize(_tmpDataContainer.size());
		for (std::size_t i = 0; i < i_src._tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_src._tmpDataContainer[i]->getNewDataContainer();
	}

	virtual
	~TimeTree_Node_LeafHelper()
	{
		clear();
	}


	/*!
	 * Setup the temporary containers
	 */
	void setupTmpContainer(
			const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,	//!< Config for allocating new containers
			std::size_t i_size,		//!< Size of the temporary arrays
			int i_instanceTypeID = -1	//!< Particular type

	)
	{
		// Assume that no tmp container was allocated before
		SWEET_ASSERT(_tmpDataContainer.size() == 0);

		_tmpDataContainer.resize(i_size);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance(i_instanceTypeID);
	}

	/*!
	 * Simply set the time step size
	 */
	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;
		return true;
	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];

		_tmpDataContainer.clear();

		TimeTree_Node_Base::clear();
	}

	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		MainInteriorNodeClass *m = new MainInteriorNodeClass(static_cast<MainInteriorNodeClass&>(*this));
		return std::shared_ptr<TimeTree_Node_Base>(m);
	}

};

}}

#endif
