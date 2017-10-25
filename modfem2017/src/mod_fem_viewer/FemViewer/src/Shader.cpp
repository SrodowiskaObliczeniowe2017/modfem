/*
l182
 * Shader.cpp
 *
 *  Created on: 5 sie 2014
 *      Author: dwg
 */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include"defs.h"
#include"Log.h"
#include "fv_txt_utls.h"
#include "fv_dictstr.h"
#include"Light.h"
#include"Shader.h"
#include"Legend.h"
#include"ViewManager.h"


namespace FemViewer {


typedef struct {
	const char * vert_src;
	const char * geom_src;
	const char * frag_src;
} shader_srcs;

shader_srcs g_shaders[SH_ALL] = {
	// Wireframe rendering
	{
	"shEdge.vert",
	nullptr,
	"shEdge.frag"
	},
	// Linear rendering
	{
	"shTriVGF.vert",
	"shTriVGF.geom",
	"shTriVGF.frag"
	},
	// High order rendering
	{
	"shTriVGF.vert",
	"shTriStripVGF.geom",
	"shTriVGF.frag"
	}
};

static const char* shader_files[] = {
		"shEdge.vert", "shEdge.frag",
		"shTriVGF.vert","shTriVGF.geom","shTriVGF.frag",
		"shTriStripVGF.geom"
};

StrDict ShaderDictionary[] = {
		{EDGE_VERT, shader_files[0]},
		{EDGE_FRAG, shader_files[1]},
		{TRI_VERT, shader_files[2]},
		{TRI_GEOM, shader_files[3]},
		{TRI_FRAG, shader_files[4]},
		{TRISTRIP_GEOM, shader_files[5]}
};

const int Shader::ProjectionBlkIdx = 0;
const int Shader::IsoValuesBlkIdx = 1;

const int GLShaderProgram::MatricesBlokIndex = 0;
const int GLShaderProgram::ParametersBlockIndex = 1;

Shader::Shader(eRenderType type)
: _type(type)
, _programId(0)
, _UnifBlockIndexlocation(-1)
, _IsoValuesUnifBlk(-1)
{
	//mfp_log_debug("Shader ctr. type: %d\n",type);
}

Shader::~Shader()
{
	//mfp_log_debug("Shader::Dtr\n");
	for (vShaders::iterator it = _shaders.begin(); it != _shaders.end(); ++it) {
		glDeleteShader(*it);
	}


	if (_programId) {
		glDeleteProgram(_programId);
		_programId = 0;
	}
}

bool Shader::Init()
{
	//mfp_log_debug("Shader::Init\n");
	bool result(initGLEW());

	if (result) _programId = glCreateProgram();
	result = (_programId != 0) && (g_shaders[this->_type].vert_src != nullptr);

	if (result) {
		result = Attach(GL_VERTEX_SHADER,g_shaders[_type].vert_src);
	}
	if (result && g_shaders[this->_type].geom_src != NULL) {
		//mfp_debug("geometry shader creatin\n");
		result = Attach(GL_GEOMETRY_SHADER,g_shaders[_type].geom_src);
	}
	if (result && g_shaders[this->_type].frag_src != NULL) {
		result = Attach(GL_FRAGMENT_SHADER,g_shaders[_type].frag_src);
	}
	if (result) result = Complete();

	return result;
}

void Shader::Enable()
{
	assert(_programId > 0);
	glUseProgram(_programId);
}

//void Shader::SetMatrixMP(const Matrixf& MP)
//{
//	glUniformMatrix4fv(_MPlocation, 1, GL_FALSE, MP.matrix.data());
//}
//
//void Shader::SetMatrixMV(const Matrixf& MV)
//{
//	glUniformMatrix4fv(_MVlocation, 1, GL_FALSE, MV.matrix.data());
//}

bool Shader::Attach(GLenum Type,const char* FilePath)
{
	assert(_programId != 0);
	//mfp_debug("Loading file: %s\n",FilePath);
	std::string source;
	if (! loadFile(FilePath,source)) {
		mfp_log_err("Error while reading shader file %s\n",FilePath);
		return false;
	}
	//mfp_debug("Before create shader\n");
	GLuint shaderId = glCreateShader(Type);
	if (shaderId == 0) {
		mfp_log_err("Error while creating shader type %d\n",Type);
		return false;
	}
	_shaders.push_back(shaderId);

	//mfp_debug("before compiling shader\n");
	char const * srcPtr = source.c_str();
	glShaderSource(shaderId, 1, &srcPtr, NULL);
	glCompileShader(shaderId);

	GLint Result = GL_FALSE;

	// Check shader
	glGetShaderiv(shaderId, GL_COMPILE_STATUS, &Result);

	if (Result == GL_FALSE) {
		GLint length = 0;
		glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &length);
		std::vector<char> ShaderErrorMessage(length+1);
		glGetShaderInfoLog(shaderId, length, NULL, &ShaderErrorMessage[0]);
		fprintf(stderr,"%s\n", &ShaderErrorMessage[0]);
	}

	glAttachShader(_programId, shaderId);
	return (glGetError() == GL_NO_ERROR);
}

bool Shader::Complete()
{
	assert(_programId != 0);
	GLint Result = GL_FALSE;

	// Link program
	glLinkProgram(_programId);

	// Check the program
	glGetProgramiv(_programId, GL_LINK_STATUS, &Result);

	if (Result == GL_FALSE) {
		GLint length = 0;
		glGetProgramiv(_programId, GL_INFO_LOG_LENGTH, &length);
		std::vector<char> LinkerErrorMessage(length+1);
		glGetProgramInfoLog(_programId, length, NULL, &LinkerErrorMessage[0]);
		mfp_log_err("%s\n", &LinkerErrorMessage[0]);
		return false;
	}

	//mfp_debug("Albo tu program id = %u\n",_programId);

	for (vShaders::iterator it = _shaders.begin(); it != _shaders.end(); ++it) {
		glDeleteShader(*it);
	}

	_UnifBlockIndexlocation = GetUniformBlockIndex("Projection");
	if (_UnifBlockIndexlocation == INVALID_LOCATION) {
		return false;
	}
	glUniformBlockBinding(_programId, _UnifBlockIndexlocation, ProjectionBlkIdx);

	_IsoValuesUnifBlk = GetUniformBlockIndex("Parameters");
	if (_IsoValuesUnifBlk == INVALID_LOCATION) {
		return false;
	}
	glUniformBlockBinding(_programId, _IsoValuesUnifBlk, IsoValuesBlkIdx);
	FV_CHECK_ERROR_GL();

	//mfp_debug("After binding uniform block pf Params --- complete shader %u %u\n",_IsoValuesUnifBlk,IsoValuesBlkIdx);
	_shaders.clear();
	return (glGetError() == GL_NO_ERROR);
}

//bool Shader::SetProjection(const Matrixf& pPorj, const Matrixf& pModelView)
//{
//	return true;
//}

GLuint Shader::GetUniformParam(const char* NameOfParam)
{
	GLuint location = glGetUniformLocation(_programId, NameOfParam);
	if (location == INVALID_LOCATION) {
		mfp_log_warn("Can't get location for parameter %s\n",NameOfParam);
		//exit(-1);
	}

	return location;
}

GLuint Shader::GetUniformBlockIndex(const char* UnifBlockName)
{
	GLuint location = glGetUniformBlockIndex(_programId, UnifBlockName);
	if (location == INVALID_LOCATION) {
		mfp_log_warn("Can't get location for uniform block index of a name: %s\n",UnifBlockName);
	}

	return location;
}

EdgeShader::EdgeShader()
: Shader(SH_EDGE)
{
	//mfp_log_debug("EdgeShader ctr\n");
}


TriangleShader::TriangleShader(eRenderType type)
: Shader(type)
{
	memset(_shParams,0xFF,sizeof(_shParams));
}

bool TriangleShader::Init()
{

	bool result = Shader::Init();
	if (!result) return result;

	// Read parameters location
	//mfp_log_debug("Init Triangle Shader\n");
	_shParams[NUM_OF_TRIFACES] = GetUniformParam("nTriangles");
	_shParams[LIGHT_POSITION]  = GetUniformParam("posLight");
	//_shParams[LIGHT_INTENISTY] = GetUniformParam("lightIntensity");
	//_shParams[AMBIENT_INTENSITY] = GetUniformParam("ambientIntensity");
	//_shParams[DRAW_EDGES] = GetUniformParam("bDrawEdges");
	//_shParams[DRAW_ISOLINES] = GetUniformParam("bDrawIsoLines");

	if (_shParams[NUM_OF_TRIFACES] == INVALID_LOCATION ||
		_shParams[LIGHT_POSITION] == INVALID_LOCATION //||
		//_shParams[LIGHT_INTENISTY] == INVALID_LOCATION ||
		//_shParams[AMBIENT_INTENSITY] == INVALID_LOCATION ||
		//_shParams[DRAW_EDGES] == INVALID_LOCATION ||
		//_shParams[DRAW_ISOLINES] == INVALID_LOCATION
		) {
		result = false;
	}

	return result;
}

//void TriangleShader::SetMatrixNM(const float normMatrix[])
//{
//	//assert(_NormalModelViewMatrixUnif != INVALID_LOCATION);
//	//glUniformMatrix3fv(	_NormalModelViewMatrixUnif, 1, GL_FALSE, normMatrix);
//}

void TriangleShader::SetNumOfTriangles(int nTriangles)
{
	//mfp_debug("Number of trianglefaces: %u\n",nTriangles);
	assert(_shParams[NUM_OF_TRIFACES] != INVALID_LOCATION);
	glUniform1i(_shParams[NUM_OF_TRIFACES], nTriangles);
}

void TriangleShader::SetLight(const Light& light)
{
	assert(_shParams[LIGHT_POSITION] != INVALID_LOCATION);
	//assert(_shParams[LIGHT_INTENISTY] != INVALID_LOCATION);
	//assert(_shParams[AMBIENT_INTENSITY] != INVALID_LOCATION);
	//glUniform2f(_WinScale, WinScale[0], WinScale[1]);
	//FV_CHECK_ERROR_GL();
	//mfp_debug("Setting light position: %f %f %f\n",light.Position().x,light.Position().y,light.Position().z);
	glUniform3fv(_shParams[LIGHT_POSITION], 1, light.Position().v);
	//glUniform3fv(_shParams[LIGHT_INTENISTY], 1, light.DiffuseIntensity().v);
	//glUniform3fv(_shParams[AMBIENT_INTENSITY], 1, light.AmbientIntensity().v);
	//glUniform3fv(_LDirlocation, 1, light.Direction().v);
}

void TriangleShader::EdgesOn(bool flag)
{
	//mfp_debug("Setting bDrawEdges %u to : %d\n",_DrawEdgesUnif,int(flag));
	assert(_shParams[DRAW_EDGES] != INVALID_LOCATION);
	glUniform1i(_shParams[DRAW_EDGES], flag ? 1 : 0);
}

void TriangleShader::IsoLinesOn(bool flag)
{
	//mfp_debug("Setting bDrawIsovalues %u to : %d\n",_DrawIsoLinesUnif,int(flag));
	assert(_shParams[DRAW_ISOLINES] != INVALID_LOCATION);
	glUniform1i(_shParams[DRAW_ISOLINES], flag ? 1 : 0);
}

void TriangleShader::ColoredIsoLines(bool flag)
{
	//mfp_debug("Setting bColoredIsoLines %u to : %d\n",_ColoredIsoLinesUnif,int(flag));
	assert(_shParams[ISOLINES_COLORED] != INVALID_LOCATION);
	glUniform1i(_shParams[ISOLINES_COLORED], flag ? 1 : 0);
}


const char* GLShader::GetPath(int type)
{
	const char* path(0);
	if (type >= EDGE_VERT && type <= TRISTRIP_GEOM) {
		const unsigned int size = FV_SIZEOF_ARRAY(ShaderDictionary);
		path = getString(ShaderDictionary,size,type);
	}
	return path;
}

GLShader::GLShader(const char* path,int type)
{
	//mfp_log_debug("GLShader ctr");
	std::string source;
	if (! loadFile(path, source)) {
		mfp_log_err("Error while reading shader file %s\n",path);
		throw;
	}

	//mfp_debug("Before create shader\n");
	m_id = glCreateShader(type);
	if (m_id == 0) {
		mfp_log_err("Error while creating shader type %d\n",type);
		throw;
	}

	char const * srcPtr = source.c_str();
	glShaderSource(m_id, 1, &srcPtr, NULL);
	glCompileShader(m_id);

	GLint Result = GL_FALSE;

	// Check shader
	glGetShaderiv(m_id, GL_COMPILE_STATUS, &Result);

	if (Result == GL_FALSE) {
		GLint length = 0;
		glGetShaderiv(m_id, GL_INFO_LOG_LENGTH, &length);
		std::vector<char> ShaderErrorMessage(length+1);
		glGetShaderInfoLog(m_id, length, NULL, &ShaderErrorMessage[0]);
		mfp_log_err("%s\n", &ShaderErrorMessage[0]);
		throw;
	}
	m_type = type;
	m_loaded = true;
}

GLShader::GLShader(const GLShader& rh)
: m_id(rh.m_id)
, m_type(rh.m_type)
, m_loaded(rh.m_loaded)
{

}

void GLShader::Delete()
{
	if (!m_loaded) return;
	m_loaded = false;
	glDeleteShader(m_id);
}

GLShaderProgram::GLShaderProgram()
: m_shaders()
, m_linked(false)
, m_id(0)
, m_blkProjID(INVALID_LOCATION)
, m_blkParamID(INVALID_LOCATION)
{
	//mfp_log_debug("GLShaderProgram ctr\n");
	m_id = glCreateProgram();
	FV_CHECK_ERROR_GL();
}

GLShaderProgram::~GLShaderProgram()
{
	//mfp_log_debug("GLShaderProgram dtr\n");
	Delete();
}


void GLShaderProgram::Delete()
{
	for (vShaders::iterator it = m_shaders.begin(); it != m_shaders.end(); ++it) {
		it->Delete();
	}
	m_shaders.clear();
	glDeleteProgram(m_id);
}

bool GLShaderProgram::AddShader(const char* path,int type)
{
	bool ret = false;
	auto it = std::find_if(m_shaders.begin(),m_shaders.end(),
			[type](GLShader const& s)
			{
				return s.GetType() == type;
			});
	if (it == m_shaders.end())  {
		//mfp_log_debug("Adding shader\n");
		GLShader shader(path,type);
		glAttachShader(m_id, shader.GetID());

		ret = (glGetError() == GL_NO_ERROR);
		if (ret) m_shaders.push_back(shader);
	}
	return ret;
}

bool GLShaderProgram::Complete()
{
	assert(m_id != 0);
	GLint Result = GL_FALSE;

	// Link program
	glLinkProgram(m_id);

	// Check the program
	glGetProgramiv(m_id, GL_LINK_STATUS, &Result);

	if (Result == GL_FALSE) {
		GLint length = 0;
		glGetProgramiv(m_id, GL_INFO_LOG_LENGTH, &length);
		std::vector<char> LinkerErrorMessage(length+1);
		glGetProgramInfoLog(m_id, length, NULL, &LinkerErrorMessage[0]);
		//mfp_log_err("%s\n", &LinkerErrorMessage[0]);
		return false;
	}

	//mfp_debug("Albo tu program id = %u\n",_programId);

	for (vShaders::iterator it = m_shaders.begin(); it != m_shaders.end(); ++it) {
		it->Delete();
	}

	m_blkProjID = GetUniformBlockIndex("Projection");
	if (m_blkProjID == INVALID_LOCATION) {
		return false;
	}

	glUniformBlockBinding(m_id, m_blkProjID, MatricesBlokIndex);

	m_blkParamID = GetUniformBlockIndex("Parameters");
	if (m_blkParamID == INVALID_LOCATION) {
		return false;
	}

	glUniformBlockBinding(m_id, m_blkParamID, ParametersBlockIndex);
	FV_CHECK_ERROR_GL();

	//mfp_debug("After binding uniform block pf Params --- complete shader %u %u\n",_IsoValuesUnifBlk,IsoValuesBlkIdx);
	m_shaders.clear();
	m_linked = (glGetError() == GL_NO_ERROR);
	return m_linked;

}

void GLShaderProgram::Enable()
{
	assert(m_id != 0);
	if (m_linked) glUseProgram(m_id);
}

GLuint GLShaderProgram::GetID()
{
	return m_id;
}

GLuint GLShaderProgram::GetUniformBlockIndex(const char* blkName)
{
	GLuint location = glGetUniformBlockIndex(m_id, blkName);
	if (location == INVALID_LOCATION) {
		mfp_log_warn("Can't get location for uniform block index of a name: %s\n",blkName);
	}

	return location;
}

void GLShaderProgram::SetUniform(const char* name, int value)
{
	assert(m_id!=0);
	int iLoc = glGetUniformLocation(m_id, name);
	//FV_CHECK_ERROR_GL();
	assert(iLoc!=-1);
	glUniform1i(iLoc, value);
}


void GLShaderProgram::SetUniform(const char* name,const float* values)
{
	assert(m_id!=0);
	int iLoc = glGetUniformLocation(m_id, name);
	//FV_CHECK_ERROR_GL();
	assert(iLoc!=-1);
	glUniform3fv(iLoc, 1, values);
}



}// end namespace


