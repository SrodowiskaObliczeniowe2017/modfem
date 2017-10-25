/*
 * Shader.h
 *
 *  Created on: 5 sie 2014
 *      Author: dwg
 */

#ifndef SHADER_H_
#define SHADER_H_

#include<vector>
#include"fv_inc.h"
#include"Log.h"


struct tagStrDict;
namespace FemViewer {

extern tagStrDict ShaderDictionary[];

class Light;

typedef enum eRenderType {
	SH_UNKNOWN	= -1,
	SH_EDGE		= 0,
	SH_TRIANGLE,
	SH_TRIANGLE_STRIP,
	SH_ALL
} glsl_program_type;

typedef enum {
	EDGE_VERT =0,
	EDGE_FRAG,
	TRI_VERT,
	TRI_GEOM,
	TRI_FRAG,
	TRISTRIP_GEOM
} glslshader_type;

class Shader
{
  public:
	static const int ProjectionBlkIdx;
	static const int IsoValuesBlkIdx;

  public:
	Shader(eRenderType type = SH_UNKNOWN);
	virtual ~Shader();
	virtual bool Init();
	void Enable();
	//static Matrixf MP;
	//static Matrixf MV;
	//bool SetProjection(const Matrixf& pPorj, const Matrixf& pModelView);
	//void SetMatrixMP(const Matrixf& MP);
	//void SetMatrixMV(const Matrixf& MV);
  protected:
	bool 	Attach(GLenum Type,const char* FilePath);
	bool 	Complete();
	GLuint 	GetUniformParam(const char* NameOfParam);
	GLuint 	GetUniformBlockIndex(const char* UnifBlockName);

	eRenderType _type;
	GLuint  _programId;
	GLuint  _UnifBlockIndexlocation;
	GLuint  _IsoValuesUnifBlk;
  private:
	typedef std::vector<GLuint> vShaders;
	vShaders _shaders;
};

class EdgeShader : public Shader
{
  public:
	EdgeShader();
	~EdgeShader() {mfp_log_debug("~EdgeShader");}

  private:
	EdgeShader(const EdgeShader&);
	EdgeShader& operator=(const EdgeShader&);
};

class TriangleShader : public Shader
{
public:
	TriangleShader(eRenderType type = SH_TRIANGLE);
	virtual ~TriangleShader() {mfp_log_debug("~TriangleShader");}
	virtual bool Init();

	void SetMatrixNM(const float martix[]);
	void SetNumOfTriangles(int nTraingles);
	void SetLight(const Light& light);
	void EdgesOn(bool flag = true);
	void IsoLinesOn(bool flag = false);
	void ColoredIsoLines(bool flag = false);

protected:
	enum {
		NUM_OF_TRIFACES = 0,
		LIGHT_POSITION,
		LIGHT_INTENISTY,
		AMBIENT_INTENSITY,
		DRAW_EDGES,
		DRAW_ISOLINES,
		ISOLINES_COLORED,
		ALL_PARAMS
	};

	GLuint _shParams[ALL_PARAMS];
private:
	TriangleShader(const TriangleShader&);
	TriangleShader& operator=(const TriangleShader&);
};


class TriStripsVGFShader : public TriangleShader
{
public:
	TriStripsVGFShader() : TriangleShader(SH_TRIANGLE_STRIP)
	{}
	~TriStripsVGFShader() { mfp_log_debug("~TriStripShader"); }
private:
	TriStripsVGFShader(const TriStripsVGFShader&);
	TriStripsVGFShader& operator=(const TriStripsVGFShader&);
};

class GLShader {
public:
	static const char* GetPath(int type);
public:
	GLShader(const char* path,int type);
	GLShader(const GLShader& rh);
	void Delete();
	bool IsLoaded() const { return m_loaded; }
	GLuint GetID() const { return m_id; }
	GLuint GetType() const { return m_type; }
private:
	GLuint m_id;
	int m_type;
	bool m_loaded;
};


class GLShaderProgram
{
public:
	static const int MatricesBlokIndex;
	static const int ParametersBlockIndex;
public:
	GLShaderProgram();
	~GLShaderProgram();

	//void Create();
	void Delete();
	bool AddShader(const char* path,int type);
	bool Complete();
	void Enable();

	GLuint GetID();

	// Services
	GLuint GetUniformBlockIndex(const char* blkName);

	void SetUniform(const char* name,int value);
	void SetUniform(const char* name,const float* values);


private:
	GLShaderProgram(const GLShaderProgram&);
	GLShaderProgram& operator=(const GLShaderProgram&);

	typedef std::vector<GLShader> vShaders;
	vShaders m_shaders;
	bool m_linked;
	GLuint m_id;

	GLuint m_blkProjID;
	GLuint m_blkParamID;

};

}// end namespace
#endif /* SHADER_H_ */
