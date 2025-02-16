/*
 * Parareal_CoutPrefix.hpp
 *
 *  Created on: 18 Apr 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_CONSOLEPREFIX_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PARAREAL_CONSOLEPREFIX_HPP


#include <iostream>

namespace sweet {
namespace DEPRECATED_pint {

/**
 * Prefix the std::cout output with a given string
 */
class Parareal_ConsolePrefix
{
#if 1
private:
	/*
	 * Prefix std::cout with string
	 *
	 * Source: http://stackoverflow.com/questions/27336335/c-cout-with-prefix
	 */
	class PrefixBuffer
		: public std::streambuf
	{
		std::string     prefix;
		std::streambuf* sbuf;
		bool            need_prefix;

		int sync()
		{
			return sbuf->pubsync();
		}

		int overflow(int c)
		{
			if (c != std::char_traits<char>::eof()) {
				if (need_prefix
					&& !prefix.empty()
					&& (int)prefix.size() != sbuf->sputn(&prefix[0], prefix.size())) {
					return std::char_traits<char>::eof();
				}
				need_prefix = c == '\n';
			}
			return sbuf->sputc(c);
		}

	public:
		PrefixBuffer()	:
			sbuf(0),
			need_prefix(true)
		{
		}

	public:
		void setPrefix(const std::string& i_prefix)
		{
			prefix = i_prefix;
		}

	public:
		void setPrefix(const char *i_prefix)
		{
			prefix = i_prefix;
		}


	public:
		void setStreamBuf(std::streambuf* i_sbuf)
		{
			sbuf = i_sbuf;
		}
	};

	PrefixBuffer coutPrefixBuffer, cerrPrefixBuffer;
#endif

	std::streambuf *coutbuf, *cerrbuf;

public:
	void start(
			int i_number		//!< prefix number will be extended to "[NUMBER] "
	)
	{
		std::ostringstream ss;
		ss << "[" << i_number << "] ";

		coutPrefixBuffer.setPrefix(ss.str());
		cerrPrefixBuffer.setPrefix(ss.str());

		std::cout.rdbuf(&coutPrefixBuffer);
		std::cerr.rdbuf(&cerrPrefixBuffer);
	}


	void start(
			const char *i_prefix	//!< prefix string
	)
	{
		coutPrefixBuffer.setPrefix(i_prefix);
		std::cout.rdbuf(&coutPrefixBuffer);

		cerrPrefixBuffer.setPrefix(i_prefix);
		std::cerr.rdbuf(&cerrPrefixBuffer);
	}



	void end()
	{
		std::cout.rdbuf(coutbuf);
		std::cerr.rdbuf(cerrbuf);
	}


	Parareal_ConsolePrefix()
	{
		coutbuf = std::cout.rdbuf();
		cerrbuf = std::cerr.rdbuf();

		coutPrefixBuffer.setStreamBuf(coutbuf);
		cerrPrefixBuffer.setStreamBuf(cerrbuf);
	}


	~Parareal_ConsolePrefix()
	{
		end();
	}

};


}}

#endif
