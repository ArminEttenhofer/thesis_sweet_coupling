#ifndef INCLUDE_SWEET_TIMETREE_TIMETREEIR_HPP
#define INCLUDE_SWEET_TIMETREE_TIMETREEIR_HPP

#include <vector>
#include <memory>

namespace sweet {

/*!
 * \brief Time integration based on hierarchical tree composition
 */
namespace TimeTree {


/*!
 * \brief Intermediate representation (based on TimeString)
 */
class TimeTreeIR
{
public:
	TimeTreeIR()
	{
	}

	~TimeTreeIR()
	{
		clear();
	}


public:
	class Argument;
	class Function;


public:
	std::shared_ptr<TimeTreeIR::Function> mainFunction;


	void clear()
	{
		mainFunction.reset();
	}

	/*!
	 * \brief Function in time tree
	 *
	 * This represents a function of the form where parameters might be recursive again
	 *
	 * "FunctionNameAnd123SomeNumbers(these, are=parameters, and, can, also, be, a, recursive, function(yadda))"
	 */
public:
	class Function
	{
	public:
		//! Function name itself
		std::string function_name;

	public:
		//! Arguments of function
		std::vector<std::shared_ptr<Argument>> arguments;

	private:
		/*!
		 * String which contains information about the function to debug
		 * (in case of an error in parsing the arguments in the time stepper)
		 */
		std::string _debugMessage;

	public:
		Function(
				const std::string &i_function_name
		)	:
			function_name(i_function_name)
		{
		}

		~Function()
		{
			clear();
		}

	public:
		void clear()
		{
			function_name = "";
			arguments.clear();
		}

	public:
		/*
		 * Output information about this function and its arguments including
		 * recursive calls to other functions
		 */
		void print(const std::string &i_prefix_str = "")
		{
			std::cout << i_prefix_str << function_name << "(" << std::endl;
			std::string new_prefix_str = i_prefix_str + "  ";

			for (std::size_t i = 0; i < arguments.size(); i++)
				arguments[i]->print(new_prefix_str);


			std::cout << i_prefix_str << ")" << std::endl;
		}

public:
		void setDebugMessage(const std::string &i_debugMessage)
		{
			_debugMessage = i_debugMessage;
		}
		std::string getDebugMessage()
		{
			return _debugMessage;
		}
		std::string getNewLineDebugMessage()
		{
			return "\n"+_debugMessage;
		}
	};

public:
	/*!
	 * \brief An argument of a function
	 *
	 * Examples for different arguments are given here:
	 *
	 * function(function_argument(yadda),key=value_argument,fun=function_argument(yadda),value_argument)
	 */
	class Argument
	{
	public:
		Error::Base error;

		/*!
		 * Type of the argument
		 */
		enum ArgumentType
		{
			ARG_INVALID,
			ARG_TYPE_FUNCTION,	//!< just a function
			ARG_TYPE_KEY_VALUE,	//!< a key and a value (2 strings)
			ARG_TYPE_KEY_FUNCTION,	//!< a key and a function as a value
			ARG_TYPE_VALUE,		//!< only a value
		};
		ArgumentType argType;


		/*!
		 * Function in case this argument is a function
		 */
		std::shared_ptr<Function> function;

		//! Key of this argument
		std::string key;

		//! Value of this argument
		std::string value;

		/*!
		 * Argument debugging message which can be printed in case there's something
		 * wrong with this argument
		 */
private:
		std::string _debugMessage;


public:
		Argument()	:
			argType(ARG_INVALID)
		{
		}

		void reset()
		{
			key = "";
			value = "";
		}

		/*!
		 * Return the value assuming it's boolean
		 */
		bool getValue(bool &o_value)
		{
			if (value == "true")
			{
				o_value = true;
				return true;
			}

			if (value == "false")
			{
				o_value = false;
				return true;
			}

			/*
			 * Interpret 0 as false and everything else as true
			 */
			int value_int;
			try
			{
				value_int = std::stoi(value);
			}
			catch (const std::exception &e)
			{
				return error.set("Exception caught during conversion of value '"+value+"' to integer: "+e.what());
			}

			o_value = (value_int != 0);
			return true;
		}

		/*!
		 * Return the value assuming it's an integer
		 */
		bool getValue(int &o_value)
		{
			try
			{
				o_value = std::stoi(value);
			}
			catch (const std::exception &e)
			{
				error.set("Exception caught during conversion of value '"+value+"' to integer: "+e.what());
				return false;
			}

			return true;
		}

		/*!
		 * Return the value assuming it's a floating point value
		 */
		bool getValue(double &o_value)
		{
			try
			{
				o_value = std::stod(value);
			}
			catch (const std::exception &e)
			{
				error.set("Exception caught during conversion of value '"+value+"' to integer: "+e.what());
				return false;
			}

			return true;
		}


		/*!
		 * Return the value assuming it's an integer
		 */
		bool getValue(std::string &o_value)
		{
			o_value = value;

			return true;
		}

		/*!
		 * Print out the argument
		 */
	public:
		void print(const std::string &i_prefix_str = "")
		{
			std::string new_prefix_str = i_prefix_str + "  ";

			switch(argType)
			{
			case ARG_TYPE_FUNCTION:
				function->print(new_prefix_str);
				break;

			case ARG_TYPE_KEY_VALUE:
				std::cout << i_prefix_str << "'" << key << "' => '" << value << "'" << std::endl;
				break;

			case ARG_TYPE_KEY_FUNCTION:
				std::cout << i_prefix_str << "'" << key << "' => FUNCTION" << std::endl;
				function->print(new_prefix_str);
				break;

			case ARG_TYPE_VALUE:
				std::cout << i_prefix_str << "'" << value << "'" << std::endl;
				break;

			case ARG_INVALID:
				break;
			}
		}

public:
		void setDebugMessage(const std::string &i_debugMessage)
		{
			_debugMessage = i_debugMessage;
		}
		std::string getDebugMessage()
		{
			return _debugMessage;
		}
		std::string getNewLineDebugMessage()
		{
			return "\n"+_debugMessage;
		}
	};


	/**
	 * Proxy printing
	 */
	void print(const std::string i_prefix_str = "")
	{
		mainFunction->print(i_prefix_str);
	}

};

}}

#endif
