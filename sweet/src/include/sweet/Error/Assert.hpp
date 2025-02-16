#ifndef INCLUDE_SWEET_ERROR_ASSERT_HPP
#define INCLUDE_SWEET_ERROR_ASSERT_HPP

#include <sweet/Error/Fatal.hpp>
#include <string>

/*!
 * SWEET Debug assertions active during debug mode
 */

#ifdef NDEBUG

	#define SWEET_ASSERT_MSG(assertion, msg)
	#define SWEET_ASSERT(assertion)

#else

	namespace sweet {
	namespace Error {

	class _AssertInternal
	{
	public:
		_AssertInternal(
				bool i_assertion,
				const std::string &i_error_message,
				const char* i_filename,
				int i_line_no,
				const char* i_func
		)
		{
			if (i_assertion)
				return;

			sweet::Error::_Fatal("ASSERTION ERROR", i_error_message, i_filename, i_line_no, i_func);
		}
	};

	}}

	#define SWEET_ASSERT_MSG(assertion, msg)	::sweet::Error::_AssertInternal(assertion, msg, __FILE__, __LINE__, __func__)

	#define SWEET_ASSERT(assertion)			::sweet::Error::_AssertInternal(assertion, "", __FILE__, __LINE__, __func__)

#endif


#endif
