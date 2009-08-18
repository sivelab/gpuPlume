// The version numbers are required and since Apple supports version
// 1.20 of the GLSL spec, this is the version we currently support.
// By not supplying a version, the compiler assumes version 1.10.

#version 120

#extension GL_ARB_texture_rectangle : enable

uniform sampler2DRect wind_texunit;

uniform int nxdx;
uniform int nydy;
uniform int numInRow;
uniform vec3 pos;

void main(void) {
     
   // The floor of the position in 3D space is needed to find the index into
   // the 2D Texture.
   int i = int(floor(pos.y));
   int j = int(floor(pos.x));
   int k = int(floor(pos.z));
   
   // This is the initial lookup into the 2D texture that holds the wind field.
   vec2 index;
   index.s = j + mod(float(k),float(numInRow))*float(nxdx);
   index.t = i + floor(float(k)/float(numInRow))*float(nydy);
   vec3 wind = vec3(texture2DRect(wind_texunit, index));
   
   gl_FragColor = vec4(wind.x, wind.y, wind.z, 1.0);
   // gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
}