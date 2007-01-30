uniform samplerRect texture;
uniform float nx;
uniform float ny;
uniform float nz;

void main(void){
		
	vec4 pos = vec4(textureRect(texture, gl_TexCoord[0].st));

	pos.x = pos.x * (nx-1);// + 100;
        pos.z = pos.z * (nz-1);// + 100;
        pos.y = pos.y * ny;// + 100;
	   
	gl_FragColor = pos;
}