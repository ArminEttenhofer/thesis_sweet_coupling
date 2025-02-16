/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_SHACKS_DICTIONARY_HPP
#define INCLUDE_SWEET_SHACKS_DICTIONARY_HPP


#include <list>
#include <memory>
#include <typeinfo>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Base.hpp>

namespace sweet {
namespace Shacks {

/*!
 * \brief A dictionary with keys given by shack classes (inheriting 'Base').
 */
class Dictionary
{
public:
	Error::Base error;

private:
	std::list<Base*> _list;

	bool _registerationOfClassInstanceFinished;
	bool _getFinished;



public:
	Dictionary()	:
		_registerationOfClassInstanceFinished(false),
		_getFinished(false)
	{
	}


public:
	~Dictionary()
	{
		clear();
	}


	/*!
	 * Clear all data and make ready for reutilization
	 */
public:
	void clear()
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
			delete *i;

		_list.clear();
		error.clear();

		_registerationOfClassInstanceFinished = false;
		_getFinished = false;
	}


	/*!
	 * Close registration: After this, registering a new Shack triggers an error
	 */
public:
	void closeRegistration()
	{
		_registerationOfClassInstanceFinished = true;
	}


	/*!
	 * Close getting the shacks: After this, getting a Shack triggers an error
	 */
public:
	void closeGet()
	{
		_getFinished = true;
	}


	/*!
	 * Register a Shack for the first time.
	 *
	 * Note, that the Shack is handed over as a Type in the template parameter
	 */
public:
	template<typename TShack>
	bool registerFirstTime()
	{
		if (_registerationOfClassInstanceFinished)
		{
			const std::string& tname = typeid(TShack).name();
			error.set("Registration already finished (type '"+tname+"')");
			return false;
		}

		if (exists<TShack>())
		{
			const std::string& tname = typeid(TShack).name();
			error.set("Class of type '"+tname+"' already exists");
			return false;
		}

		TShack* newClass = new TShack();
		_list.push_back(newClass);
		return true;
	}

	/*!
	 * Return true if a shack with this type already exists
	 */
public:
	template<typename TShack>
	bool exists()
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			//! check whether generic interface can be casted to type T
			TShack* derived = dynamic_cast<TShack*>(*i);

			if (derived != nullptr)
				return true;
		}

		return false;
	}


	/*!
	 * Return a Shack if it exists.
	 *
	 * In case of an error, an error is set and a nullptr is returned
	 */
public:
	template<typename TShack>
	TShack* get(bool i_auto_registration = false)
	{
		if (!i_auto_registration)
		{
			if (!_registerationOfClassInstanceFinished)
			{
				const std::string& tname = typeid(TShack).name();
				error.set("Registration of class instances needs to be finished first (type '"+tname+"')");
				return nullptr;
			}
		}

		if (_getFinished)
		{
			const std::string& tname = typeid(TShack).name();
			error.set("Getting a dictionary element class already finished (type '"+tname+"')");
			return nullptr;
		}

		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			//! check whether generic interface can be casted to type T
			TShack* derived = dynamic_cast<TShack*>(*i);

			if (derived != nullptr)
				return derived;
		}

		const std::string& tname = typeid(TShack).name();
		error.set("Type '"+tname+"' not found in dictionary");
		return nullptr;
	}


	/*!
	 * Get a registered Shack
	 */
	template<typename TShack>
	bool get(TShack **o_shack, bool i_autoRegistration = false)
	{
		*o_shack = get<TShack>(i_autoRegistration);
		return o_shack != nullptr;
	}


	/*!
	 * Auto registrate this particular class if it doesn't exist and return an instance
	 */
public:
	template<typename TShack>
	TShack* getAutoRegistration()
	{
		if (!exists<TShack>())
			if (!registerFirstTime<TShack>())
				return nullptr;

		return get<TShack>(true);
	}

	/*!
	 * Auto registrate this particular class if it doesn't exist and return an instance
	 *
	 * That's an laternative version which avoids specifying the TShack again
	 */
public:
	template<typename TShack>
	TShack* getAutoRegistration(TShack** o_shack)
	{
		*o_shack = getAutoRegistration<TShack>();
		return *o_shack;
	}


	/*!
	 * Print out program arguments which are supported by each Shack
	 */
public:
	void printProgramArguments(
			const std::string& i_prefix = ""
	)
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			(*i)->printProgramArguments(i_prefix);
		}
	}

public:
	bool processProgramArguments(
			sweet::Tools::ProgramArguments &i_pa,
			bool i_skipProcessedShacks = true
	)
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			if (i_skipProcessedShacks)
				if ((*i)->argumentsProcessed)
					continue;

			if (!((*i)->processProgramArguments(i_pa)))
			{
				error.forward((*i)->error);
				return false;
			}
			(*i)->argumentsProcessed = true;
		}
		return true;
	}

public:
	void printShackData(const std::string& i_prefix = "")
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			(*i)->printShack(i_prefix);
		}
	}

};

}}

#endif
