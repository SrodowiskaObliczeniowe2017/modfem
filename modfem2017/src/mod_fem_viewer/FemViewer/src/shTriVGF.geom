#version 330
#extension GL_EXT_gpu_shader4 : enable
#extension GL_EXT_geometry_shader4 : enable
 
precision highp float;
 
layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

uniform int nTriangles;

in VertexData {
  vec4 color;
  vec4 cs_position;
  float flag;
} vert[];
 
out FragData {
  noperspective vec3 dist;
  noperspective vec4 color;
  smooth vec3 normal;
  smooth vec4 cs_position;
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

 void main()
 {
   	float Qval = (gl_PrimitiveIDIn < nTriangles) ? 0.0f : 100.0f;
    vec3 n = cross(vert[1].cs_position.xyz - vert[0].cs_position.xyz, 
                   vert[2].cs_position.xyz - vert[0].cs_position.xyz);
    n = normalize(n);
    
    gl_Position = gl_in[0].gl_Position;
    dist = vec3(1.0, Qval, 0.0);
    color = vert[0].color;
    normal = n;
    cs_position = vert[0].cs_position;
    EmitVertex();
    
    gl_Position = gl_in[1].gl_Position;
    dist = vec3(0.0, 1.0, 0.0);
    color = vert[1].color;
    normal = n;
    cs_position = vert[1].cs_position;
    EmitVertex();

    gl_Position = gl_in[2].gl_Position;
    dist = vec3(0.0, Qval, 1.0);
    color = vert[2].color;
    normal = n;
    cs_position = vert[2].cs_position;
    EmitVertex();
   
	EndPrimitive();
}
