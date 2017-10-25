#version 330
layout(std140) uniform;

in FragData {
  noperspective vec3 dist;
  noperspective vec4 color;
  smooth vec3 normal;
  smooth vec4 cs_position;
};

out vec4 outputColor;

//uniform vec3 posCamera;
uniform vec3 posLight;
uniform vec3 lightIntensity;
uniform vec3 ambientIntensity;

uniform bool bDrawEdges;
uniform bool bDrawIsoLines;
//uniform bool bColoredIsoLines;

// legend data for iso-values
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

// Calcluate edge factor
float edgeFactor(){
    vec3 d = fwidth(dist);
    vec3 a3 = smoothstep(vec3(0.0), d*1.5, dist);
    return min(min(a3.x, a3.y), a3.z);
}

// get overlap ratio of range a with range b
float compute_overlap(in float a1,in float a2,in float b1,in float b2)
{
   if (a1>=a2) return(0.0); // malformed range a
   if (b1>=b2) return(0.0); // malformed range b
   if (a1>=b2 || a2<=b1) return(0.0); // no overlap of a and b
   if (a1>=b1 && a2<=b2) return(1.0); // full overlap of a with b
   if (b1>=a1 && b2<=a2) return((b2-b1)/(a2-a1)); // full overlap of b with a
   if (a1<b1) return((a2-b1)/(a2-a1)); // partial left-side overlap of a with b
   if (a2>b2) return((b2-a1)/(a2-a1)); // partial right-side overlap of a with b
   return(0.0);
}


const float ds = 2.0f;
const float dt = 2.0f;
// compute contour overlap
float compute_contour( in float value,            // scalar function value
                       in float isoValue)         // iso value to be extracted
{
   	float slope;
   	float range2;
   	float overlap;
  	float dvds = dFdx(value);
   	float dvdt = dFdy(value);
   	// compute slope
   	slope=sqrt(dvds*dvds+dvdt*dvdt);

   	// map foot prints
   	range2=0.25*(ds+dt)*slope;

   	// compute overlap of cell with contour
   	overlap=compute_overlap(value-range2,value+range2,        // foot print range of pixel,
                           isoValue-range2,isoValue+range2); // foot print range of iso line

   	return(overlap);
}

/*vec3 unpackColor(float f) {
    vec3 color;
    color.b = floor(f / 256.0 / 256.0);
    color.g = floor((f - color.b * 256.0 * 256.0) / 256.0);
    color.r = floor(f - color.b * 256.0 * 256.0 - color.g * 256.0);
    // now we have a vec3 with the 3 components in range [0..256]. Let's normalize it!
    return color / 256.0;
}*/

void main()
{
	int nBreaks = num_breaks - 1;
	float contourOverlap;
	vec3 tColor = color.yzw;
	//vec3 tColor = unpackColor(color.y);
	if (isolines_on) {
		int iv;
		for (iv = 1; iv < nBreaks; iv++) {
		    vec4 iso_data = iso_values[iv];
			contourOverlap = compute_contour(color.x, iso_data.w);
			if (contourOverlap > 0.5f) {
				//tColor = bColoredIsoLines ? iso_data.xyz : border_col.xyz;
				tColor = iso_col.xyz;
				//tColor = border_col.xyx;
				break;
			}
		}
	}
	
	vec3 lightDir = normalize(posLight.xyz - cs_position.xyz);
	float cosAngIncidence = dot(normal, lightDir);
	cosAngIncidence = clamp(cosAngIncidence, 0, 1);
	
	vec3 tmpColor = (tColor.xyz * light_intensity.xyz * cosAngIncidence) +
		(tColor.xyz * light_ambient.xyz);
	if (edges_on) {
		outputColor = mix(border_col, vec4(tmpColor,1.0), edgeFactor());
	} else {    	
		outputColor = vec4(tmpColor,1.0);
	}
}
