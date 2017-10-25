#ifndef _FV_TXT_UTLS_H_
#define _FV_TXT_UTLS_H_

#include "types.h"
#include "fv_config.h"
#include "fv_inc.h"
#include "RenderParams.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <vector>

namespace FemViewer {

enum ePixelFormat {
  PF_UNKNOWN = 0,
  PF_R8,
  PF_RG8,
  PF_RGB8,
  PF_RGBA8,

  // Signed and noramlized
  PF_RGBA8_SRGB,
  PF_BGRA8_SRGB,

  PF_ALL
};


enum eShaderType
{
	SH_VERTEX = 0x01,
	SH_GEOMETRY = 0x02,
	SH_FRAGMENT = 0x04,
	SH_TESELATOR = 0x08,
	SH_COMPUTE = 0x10,
};



/// pixel format descriptor
typedef struct {
  int _internal; 	//< OpenGL internal format (GL_RGBA8)
  int _format; 		//< OpenGL format (GL_RGBA)
  int _type; 		//< OpenGL component type (GL_UNSIGNED_BYTE)
  unsigned int _size; //< byte size of one pixel (4)
  int _components;	//< number of components (4)
  bool _rt; 		//< true if it can be used as render target
  int _sRGB; 		//< sRGB pixel format alternative
  const char *_txt; //< readable description
  bool _compressed; //< true if it is compressed format
} pixel_descr_t, *pixel_descrptr_t;

extern const pixel_descr_t pixel_descriptors[];

inline const pixel_descr_t* getPixelFormatDescriptor(const int id) {
	assert(id > PF_UNKNOWN && id < PF_ALL);
	return (pixel_descriptors + id);
}

// struct variable to store OpenGL info
class GLCore
{
	static GLboolean wasInit;
    std::string vendor;
    std::string renderer;
    std::string version;
    std::string glslVersion;
    std::vector <std::string> extensions;
    int redBits;
    int greenBits;
    int blueBits;
    int alphaBits;
    int depthBits;
    int stencilBits;
    int maxTextureSize;
    int maxLights;
    int maxAttribStacks;
    int maxModelViewStacks;
    int maxProjectionStacks;
    int maxClipPlanes;
    int maxTextureStacks;
    int maxUniformBufferSize;

    GLuint compileList;
    GLsizei range;
    //std::vector<GLuint> bufferIDs;
    explicit GLCore();
public:

    ~GLCore();

    static GLboolean init(int argc,char *argv[]);
    static GLCore& instance();

    GLuint createBuffer(const void* data,unsigned int dataSize, GLenum target, GLenum usage);

    bool initGLLists(const GLsizei size);
    void renderGLList(const GLuint type);
    void deleteGLLists(const GLuint listID, const GLsizei size);
    void deleteBuffers(GLuint* buffer,const GLsizei size);

    void startGLList(GLuint type);
    void endGLList();

    bool getInfo(unsigned int param=0);         // extract info
    void printSelf();                           // print itself
    bool isExtensionSupported(const std::string& ext); // check if a extension is supported

    template< typename TVertex >
    GLuint createGLText(const TVertex* buffer,
    		          const GLsizei stride,
    		          const GLuint num,
					  const GLfloat color[3],
    		          const GLuint type,
    		          const GLboolean drawAll =  GL_TRUE) {}
    //void registerBuffer(const GLuint buffID);
private:
    GLCore(const GLCore&);
    GLCore& operator=(const GLCore&);

    void Check() const;// { if (!wasInit) throw fv_exception("GL Core wasn't initialized!\n"); }

};

}

bool initGLEW(const bool bFlag = false);
bool initGLLists(const GLsizei size);
GLuint initShaders(const char *path_to_vertex_shader, const char *payh_to_fragment_shader);
bool loadFile(const char* path,std::string& source);
GLuint createBuffer(const void* data,unsigned int dataSize, GLenum target, GLenum usage);

GLuint createPixelBuffer(const int width,const int height,GLuint* pbId);


// Returns width of given text in current screen coordinates
float getTextWidthOnScreen(const std::string& text, void* pCurrentFont);
void  drawText(const std::string& pText,
              float pXMin, float pYMin, float pXMax, float pYMax,
              int pFontSize, void* pCurrentFont, bool pAutoLineBreak, bool pFlagCentered);

void drawText3D(const std::string& pText, float pX, float pY, float pZ, void* pCurrentFont);
void drawString(const char *str, int x, int y, float color[4], void *font);
void drawPlane(const float* points, const float* scale, unsigned int size);

inline
bool __CheckErrorGL(const char *file, const int line) {
    bool ret_val = true;

    // check for error
    GLenum gl_error = glGetError();

    if (gl_error != GL_NO_ERROR) {
#ifdef _WIN32
        char tmpStr[512];
        // NOTE: "%s(%i) : " allows Visual Studio to directly jump to the file at the right line
        // when the user double clicks on the error line in the Output pane. Like any compile error.
        sprintf_s(tmpStr, 255, "\n%s(%i) : GL Error : %s\n\n", file, line, gluErrorString(gl_error));
        fprintf(stderr, "%s", tmpStr);
#endif
        fprintf(stderr, "GL Error in file '%s' in line %d :\n", file, line);
        fprintf(stderr, "%s\n", gluErrorString(gl_error));
        fflush(stderr);
        ret_val = false;
    }

    return ret_val;
}
#define FV_DEBUG_GL 1
#ifdef FV_DEBUG_GL
#define FV_CHECK_ERROR_GL()                                             \
    if( false == __CheckErrorGL( __FILE__, __LINE__)) {                 \
        exit(EXIT_FAILURE);												\
    }

inline
void print_array(const float array[],size_t size_of_array) {
	printf("\n");
	for (size_t i(0);i<size_of_array;++i)
		printf("tab[%u] = %f\n",i,array[i]);
}

inline
void print_array_vec4f(const float array[],
		                      const size_t size_of_array,
		                      const size_t offset) {
	printf("\n");
	size_t i(0), j(0);
	for (; i < size_of_array; i++,j+=offset) {
		printf("vec[%u]: ",i);
		print_array(&array[j],4);
	}
}

GLuint createQuadGrid(const double orig[], // minimal point of the grid
		              double dims[], // dimensions of the grid - should be 2 dims, 1 - is 0
		              double density[], // a vector of grid density
		              GLuint* nvertces // handle to number of vertices - out
		              );
GLuint createGrid3D(const float minb[],const float maxb[],const grid_t* drid, GLuint outBuff[]);
void drawGrid(const GLuint params[2],const float color[3],const float linew = 0.0);

struct colorbar_config_t;
GLuint createColorBar(const FemViewer::colorbar_config_t *cfg,GLuint* buffID);

#else
#define FV_CHECK_ERROR_GL()
#endif



#endif /* _FV_TXT_UTLS_H_
		*/
