#version 330

layout (std140) uniform; 
// xyz - position; w - flag
layout (location = 0) in vec4 inPosition;
// x - value; yzw - color
layout (location = 1) in vec4 inDiffuseColor;

out VertexData {
  vec4 color;
  vec4 cs_position;
  float flag;
} vert;

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
    // Calculate camera space position and store
	vert.cs_position = viewm * vec4(inPosition.xyz, 1.0f);
	// Calculate clip position
	gl_Position = projm * vert.cs_position;                           
	// Store color with assigned scalar value
	vert.color = inDiffuseColor;	
	// Store edge flag
	vert.flag = inPosition.w;
}
