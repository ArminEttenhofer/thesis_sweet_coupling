/*
 * PDETermRegistry.hpp
 */

#ifndef INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_REGISTRY_HPP
#define INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_REGISTRY_HPP

#include <string>
#include <map>
#include <memory>
#include <list>
#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>


namespace sweet {
namespace TimeTree {


class TimeTree_Node_Registry
{
public:
	sweet::Error::Base error;

	/*!
	 * This is to lookup for a certain implementation.
	 *
	 * There might be duplicates in it since some node can be registered under different names
	 */
private:
	std::map<std::string, std::shared_ptr<TimeTree_Node_Base>> _lookupRegistry;

	/*!
	 * Unique list where each node exists only once
	 */
private:
	std::list<std::shared_ptr<TimeTree_Node_Base>> uniqueList;

public:
	TimeTree_Node_Registry()
	{
	}

public:
	~TimeTree_Node_Registry()
	{
		clear();
	}


	/*
	 * Register shacks
	 */
public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		// insert as keys to registry
		for (auto iter : _lookupRegistry)
		{
			TimeTree_Node_Base *a = iter.second.get();
			a->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
		}

		return true;
	}

	/*!
	 * Register a new time tree node given by the template parameter
	 */
public:
	template <typename T>
	bool registerTimeTreeNode()
	{
		std::shared_ptr<TimeTree_Node_Base> b = std::make_shared<T>();

		return registerTimeTreeNode(b);
	}


	/*!
	 * Register a new time tree node which is allocated
	 */
public:
	bool registerTimeTreeNode(std::shared_ptr<TimeTree_Node_Base> i_b)
	{
		// get names of all time steppers
		const std::vector<std::string> ts = i_b->getNodeNames();

		// insert as keys to registry
		for (auto &iter : ts)
		{
			auto i = _lookupRegistry.find(iter);

			if (i != _lookupRegistry.end())
			{
				return error.set("Time tree node handler for '"+iter+"' already registered!");
			}

			_lookupRegistry[iter] = i_b;
		}

		uniqueList.push_back(i_b);

		return true;
	}

	/*!
	 * Find a tree node according to its string
	 *
	 * \return copy of tree node
	 */
	bool getTimeTreeNodeNewInstance(
			const std::string &i_ts_string,
			std::shared_ptr<TimeTree_Node_Base> &o_timestepper_instance,
			bool i_triggerError = true
	)
	{
		auto iter = _lookupRegistry.find(i_ts_string);

		if (iter == _lookupRegistry.end())
		{
			if (i_triggerError)
				return error.set("Timestepper with string '"+i_ts_string+"' not found");
			else
				return false;
		}

		o_timestepper_instance = iter->second->getInstanceCopy();
		return true;
	}

	void clear()
	{
		_lookupRegistry.clear();
		uniqueList.clear();
	}


	/*!
	 * Print out help for all registered nodes
	 */
	virtual
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = "",
			int i_verbosity = 0
	)
	{
		for (auto iter : uniqueList)
		{
			TimeTree_Node_Base *a = iter.get();
			a->outputHelp(o_ostream, i_prefix);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
			o_ostream << i_prefix << std::endl;


			if (i_verbosity > 0)
			{
				std::cout << i_prefix << "    PROVIDES: " << a->getSupportedEvalsAsString() << std::endl;
			}

			o_ostream << " -------------------------------------------------------------------------------" << std::endl;
		}
		return true;
	}
};

}}

#endif
