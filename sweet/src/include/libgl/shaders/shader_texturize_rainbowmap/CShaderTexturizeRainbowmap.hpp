#ifndef INCLUDE_SWEET_LIBGL_SHADERS_SHADER_TEXTURIZE_RAINBOWMAP_CSHADERTEXTURIZERAINBOWMAP_HPP
#define INCLUDE_SWEET_LIBGL_SHADERS_SHADER_TEXTURIZE_RAINBOWMAP_CSHADERTEXTURIZERAINBOWMAP_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlTexture.hpp>

namespace libgl {

class GlShaderTexturizeRainbowmap	: public GlProgram
{
public:
	GlUniform pvm_matrix_uniform;

	GlShaderTexturizeRainbowmap()
	{
		if (!initVertFragShadersFromDirectory("shader_texturize_rainbowmap"))
			return;

		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: linking: " << infoLog << std::endl;
			return;
		}

		use();

			setupUniform(pvm_matrix_uniform, "pvm_matrix");

			bindAttribLocation(0, "vertex_position4");
			bindAttribLocation(1, "vertex_normal3");
			bindAttribLocation(2, "vertex_texture_coord2");

		disable();
	}

	~GlShaderTexturizeRainbowmap()
	{
	}
};

}

#endif
