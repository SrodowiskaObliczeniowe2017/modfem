#version 330

layout(std140) uniform;

layout (location = 0) in vec3 inPos;

out vec4 color;

uniform Projection
{
  mat4 proj;
  mat4 view; 
};

const int numberOfBreaks = 32;
uniform Parameters
{
	mat4 projm;
	mat4 viewm;
	vec4 wireframe_col;
	vec4 border_col;
	vec4 iso_col;
	vec4 iso_values[numberOfBreaks];
	int num_breaks;
	int num_triangles;
	bool edges_on;
	bool isolines_on;
	vec4 light_intensity;
	vec4 light_ambient;
};

void main(void) {                              
	vec4 temp = viewm * vec4(inPos, 1.0);
	gl_Position = projm * temp;
	color = wireframe_col;
}
 

