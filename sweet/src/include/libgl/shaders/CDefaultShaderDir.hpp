/*
 * CDefaultShaderDir.hpp
 *
 *  Created on: Aug 4, 2011
 * Author: schreibm
 */

#ifndef INCLUDE_SWEET_LIBGL_SHADERS_CDEFAULTSHADERDIR_HPP
#define INCLUDE_SWEET_LIBGL_SHADERS_CDEFAULTSHADERDIR_HPP

#include <string>
#include <stdlib.h>

#ifndef SHADER_GLSL_DEFAULT_DIR
#define SHADER_GLSL_DEFAULT_DIR "data/shaders_glsl/"
#endif

namespace libgl {

class ShaderDir
{
public:
	static std::string getDirectory()
	{
		std::string sweet_root;

		char *a = getenv("MULE_SOFTWARE_ROOT");
		if (a == nullptr)
		{
			sweet_root = "./";
		}
		else
		{
			sweet_root = a;
			sweet_root = sweet_root + "/";
		}

		sweet_root += SHADER_GLSL_DEFAULT_DIR;

		return sweet_root;
	}
};

}

#endif
