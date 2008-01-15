void main(void){	

	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
	
	//vec4 pos = vec4(gl_Vertex);
	//gl_Position = gl_ModelViewProjectionMatrix*pos;
}