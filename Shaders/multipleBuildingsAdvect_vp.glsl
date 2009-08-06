// The version numbers are required and since Apple supports version
// 1.20 of the GLSL spec, this is the version we currently support.
// By not supplying a version, the compiler assumes version 1.10.

#version 120


void main(void){	

	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position = ftransform();
}