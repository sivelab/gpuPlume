// The version numbers are required and since Apple supports version
// 1.20 of the GLSL spec, this is the version we currently support.
// By not supplying a version, the compiler assumes version 1.10.

#version 120

#extension GL_ARB_texture_rectangle : enable

// These matricies are for indexing into the depth map taken from
// the sun's position.
uniform mat4x4 sunModelviewMatrix;
uniform mat4x4 sunProjectionMatrix;
uniform mat4x4 sunScaleAndBiasMatrix;

// This is the z position, the shader will get run multiple times for
// each z slice.
uniform float zPos;

// This is the depth map taken from the sun's position.
uniform sampler2DShadow shadowMap;

// This is the input positions.
uniform sampler2DRect posTex;

void main(void){
     float average = 0.0;
     int count = 0;

     float i;
     float j;
     float k;

     for(i = 0; i < 1; i += 0.1) {
       for(j = 0; j < 1; j += 0.1) {
         for(k = 0; k < 1; k += 0.1) {
	   vec4 pos = vec4(gl_FragCoord.x + i , gl_FragCoord.y + j, zPos + k, 1.0);
	   vec4 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
	   average += shadow2DProj(shadowMap, index).x;
	   count++;
	 }
       }
     }
     
     average = average / count;

     gl_FragColor = vec4(average, average, average, 1.0);
}