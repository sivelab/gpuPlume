//uniform samplerRect render_texunit;
varying vec4 pcolor;

void main(void)
{
	/*	
	vec4 pos = vec4(textureRect(render_texunit, gl_TexCoord[0].st));
	
	pos.x = 1.0;
	pos.y = 1.0;
	pos.z = 1.0;
	pos.a = 1.0;*/
        
	//Draw particles white
	//gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
	//Draw particles black
	//gl_FragColor = vec4(0.0,0.0,0.0,1.0);
	gl_FragColor = pcolor;

}