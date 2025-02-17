#ifndef INCLUDE_SWEET_LIBGL_SHADERS_SHADER_BLINN_SHADOW_MAP_CSHADERBLINNSHADOWMAP_HPP
#define INCLUDE_SWEET_LIBGL_SHADERS_SHADER_BLINN_SHADOW_MAP_CSHADERBLINNSHADOWMAP_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlProgram.hpp>
#include <libgl/core/GlTexture.hpp>
#include <libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp>

namespace libgl {

class GlShaderBlinnShadowMap	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	//!< uniform to enable and disable texturing
	GlUniform shadow_map_matrix_uniform;	//!< shadow map matrix


	GlShaderBlinnShadowMap()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_blinn_shadow_map");
		attachFragShader((ShaderDir::getDirectory()+"shader_blinn/fragment_shader_skeleton.glsl").c_str());

		// link programs
		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: during linking: " << infoLog << std::endl;
			return;
		}

		setupUniform(texture0_enabled, "texture0_enabled");
		setupUniform(shadow_map_matrix_uniform, "shadow_map_matrix");

		use();
		setUniform1i("texture0", 0);
		setUniform1i("texture_shadow_map", 1);
		disable();

		initBlinnSkeleton(*this);
	}

	~GlShaderBlinnShadowMap()
	{
	}

	/**
	 * setup the light for rendering
	 */
	void setupUniforms(	GlMaterial	&material,
						Lights &lights,
						const vec3 &light_view_pos3,
						const mat4 &shadow_map_matrix
	)
	{
		CShaderBlinnSkeleton::setupUniformsLights(lights, light_view_pos3);
		setupUniformsMaterial(material);
		setupUniformsShadowMapping(shadow_map_matrix);
	}

	/**
	 * setup the uniforms for rendering
	 */
	void setupUniformsShadowMapping(	const mat4 &shadow_map_matrix	)
	{

		shadow_map_matrix_uniform.set(shadow_map_matrix);
	}

	/**
	 * setup the uniforms for rendering
	 */
	void setupUniformsMaterial(	GlMaterial	&material	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};

}

#endif
