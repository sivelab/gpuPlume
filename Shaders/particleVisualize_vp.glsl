//uniform samplerRect vel;
varying vec4 pcolor;

void main(void)
{
	vec4 pos = vec4(gl_Vertex);
	pos.w = 1.0;


	pcolor = gl_Color;
	//gl_Color = vec4(1.0, 1.0, 1.0, 1.0);
	//pcolor = vec3(textureRect(vel,gl_TexCoord[0].xy));
	//pcolor = normalize(pcolor);

	gl_Position = gl_ModelViewProjectionMatrix * pos;
	//gl_TexCoord[0] = gl_MultiTexCoord0;
	//gl_Position = ftransform();
}