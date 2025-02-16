#ifndef INCLUDE_SWEET_LIBGL_SHADERS_SHADER_TEXTURIZE_CSHADERTEXTURIZE_HPP
#define INCLUDE_SWEET_LIBGL_SHADERS_SHADER_TEXTURIZE_CSHADERTEXTURIZE_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlTexture.hpp>

namespace libgl {

class GlShaderTexturize	: public GlProgram
{
public:
	GlUniform pvm_matrix_uniform;

	GlShaderTexturize()
	{
		if (!initVertFragShadersFromDirectory("shader_texturize"))
			return;

		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: linking: " << infoLog << std::endl;
			return;
		}

		setupUniform(pvm_matrix_uniform, "pvm_matrix");
	}

	~GlShaderTexturize()
	{
	}
};

}

#endif
