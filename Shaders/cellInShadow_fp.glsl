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

// This tells the shader which points within the cell to find the shadow values
// for.
uniform int cellPoints;

uniform float dx;

uniform float dy;

uniform float dz;

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

// Sample 100 times within a cell (roughly)...		 
//     for(i = 0.05; i < 1; i += 0.1) {
//       for(j = 0.05; j < 1; j += 0.1) {
//         for(k = 0.05; k < 1; k += 0.1) {
//	   				vec4 pos = vec4(gl_FragCoord.x + i , gl_FragCoord.y + j, zPos + k, 1.0);
//	   				vec4 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
//	   				average += shadow2DProj(shadowMap, index).x;
//	   				count++;
//	 				}
//       }
//     }
     
//     average = average / count;

vec4 pos;
vec4 index;
float pos1, pos2, pos3;

  if(cellPoints == 0) {
  
		 vec4 pos = vec4(gl_FragCoord.x, gl_FragCoord.y, zPos - 0.49 * dz, 1.0);
		 vec4 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
		 pos1 = shadow2DProj(shadowMap, index).x;
		 
		 pos = vec4(gl_FragCoord.x - 0.49 * dx, gl_FragCoord.y, zPos, 1.0);
		 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
		 pos2 = shadow2DProj(shadowMap, index).x;
		 
		 pos = vec4(gl_FragCoord.x, gl_FragCoord.y - 0.49 * dy, zPos, 1.0);
		 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
		 pos3 = shadow2DProj(shadowMap, index).x;
  
  } else {
  
  	 pos = vec4(gl_FragCoord.x, gl_FragCoord.y, zPos + 0.49 * dz, 1.0);
		 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
		 pos1 = shadow2DProj(shadowMap, index).x;
		 
		 pos = vec4(gl_FragCoord.x + 0.49 * dx, gl_FragCoord.y, zPos, 1.0);
		 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
		 pos2 = shadow2DProj(shadowMap, index).x;
		 
		 pos = vec4(gl_FragCoord.x, gl_FragCoord.y + 0.49 * dy, zPos, 1.0);
		 index = (sunScaleAndBiasMatrix * sunProjectionMatrix * sunModelviewMatrix) * pos;
		 pos3 = shadow2DProj(shadowMap, index).x;
  
  }

     gl_FragColor = vec4(pos1, pos2, pos3, 1.0);

}
