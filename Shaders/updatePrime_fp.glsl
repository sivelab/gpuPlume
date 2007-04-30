uniform samplerRect primePrev;
uniform samplerRect pos;
uniform samplerRect wind;
uniform samplerRect random;
uniform samplerRect lambda;
uniform float time_step;
uniform int nx;
uniform int ny;
uniform int nz;
uniform float numInRow;

//
// This variable contains a uniformally generated random 2D vector in
// the range [0,1].  The value is generated on the CPU each iteration
// and is used to offset the random texture (used in this shader).  By
// offsetting the random texture lookup, the random numbers will
// reflect a more uniform random value rather than the same value for
// each particle.
//
uniform int random_texWidth;
uniform int random_texHeight;
uniform vec2 random_texCoordOffset;

void main(void)
{
	vec2 texCoord = gl_TexCoord[0].xy;
	
	// Add the random texture coordinate offset to the
	// texture coordinate.  The texture is set to perform wrapping
	// so texture coordinates outside the range will be valid.
	//
	// GL_TEXTURE_RECTANGLE_ARB does not support GL_REPEAT (at least on windows)
	// so we will do the normalization here given the width and length of the 
	// random texture.
	vec2 rTexCoord = texCoord + random_texCoordOffset;

	// bring the texture coordinate back within the (0,W)x(0,H) range
	if (rTexCoord.s > random_texWidth)
	   rTexCoord.s = rTexCoord.s - random_texWidth;
	if (rTexCoord.t > random_texHeight)
	   rTexCoord.t = rTexCoord.t - random_texHeight;

	// lookup the random value to be used for this particle in
	// this timestep
	vec3 randn = vec3(textureRect(random, rTexCoord));

	float xRandom = randn.x;
	float yRandom = randn.y;
	float zRandom = randn.z;
	
	vec3 PrmPrev = vec3(textureRect(primePrev, texCoord));
	vec3 PrmCurr=PrmPrev;
	
	float upPrev=PrmPrev.x;
	float vpPrev=PrmPrev.y;
	float wpPrev=PrmPrev.z;	
		
	vec4 position = vec4(textureRect(pos, texCoord));
           
    int i = floor(position.y);
    int j = floor(position.x);
    int k = floor(position.z);
    
    if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >=0)){
   

    //This is the initial lookup into the 2D texture that holds the wind field.
 	vec2 index;
   	index.s = j + mod(k,numInRow)*nx;
        index.t = i + floor(k/numInRow)*ny;
	
	vec4 windTex = vec4(textureRect(wind, index));
	vec4 lam = vec4(textureRect(lambda, index));
	
	float Lam11=lam.x;
	float Lam22=lam.y;
	float Lam33=lam.z;
	float Lam13=lam.w;
	
	float CoEps_D2=windTex.w; // grabing 4th vector,--Co*Eps/2 in the wind texture
      
        float du= -CoEps_D2*(Lam11*upPrev+Lam13*wpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* xRandom;
        float dv= -CoEps_D2*(Lam22*vpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* yRandom;
        float dw= -CoEps_D2*(Lam13*upPrev+Lam33*wpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* zRandom;
           
        PrmCurr=vec3(upPrev+du,vpPrev+dv,wpPrev+dw);
	
    }	
	//Currently this code just keeps passing the prime values,
	//keeping them the same.  
	//We need to add the equations to calculate the new values.
    gl_FragColor = vec4(PrmCurr, 1.0);


}