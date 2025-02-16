#ifndef INCLUDE_SWEET_LIBGL_SHADERS_SHADER_BATHYMETRY_CSHADERBATHYMETRY_HPP
#define INCLUDE_SWEET_LIBGL_SHADERS_SHADER_BATHYMETRY_CSHADERBATHYMETRY_HPP

#include <libgl/shaders/CDefaultShaderDir.hpp>
#include <libgl/core/GlTexture.hpp>
#include <libgl/core/CGlError.hpp>


/**
 * general blinn shader to use for rendering vertices
 */
#include <libgl/core/GlProgram.hpp>
#include <libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp>

namespace libgl {

class CShaderBathymetry	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	//!< uniform to enable and disable texturing

	CShaderBathymetry()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_bathymetry");
		attachFragShader(ShaderDir::getDirectory()+"shader_blinn/fragment_shader_skeleton.glsl");

		// link programs
		link();
		if (error())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: during linking: " << infoLog << std::endl;
			return;
		}

		initBlinnSkeleton(*this);
	}

	~CShaderBathymetry()
	{
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniforms(
			GlMaterial	&material,
			CGlLights &lights,
			const vec3 &light_view_pos3
	)
	{
		CShaderBlinnSkeleton::setupUniforms(material, lights, light_view_pos3);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniformsMaterial(
			GlMaterial	&material
	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};

}

#endif
