uniform samplerRect primePrev;
uniform samplerRect pos;
uniform samplerRect wind;
uniform samplerRect random;
uniform samplerRect lambda;
uniform float time_step;
uniform float nx;
uniform float ny;
uniform float nz;
uniform float numInRow;

//
// This variable contains a uniformally generated random 2D vector in
// the range [0,1].  The value is generated on the CPU each iteration
// and is used to offset the random texture (used in this shader).  By
// offsetting the random texture lookup, the random numbers will
// reflect a more uniform random value rather than the same value for
// each particle.
//
uniform vec2 random_texCoordOffset;

void main(void)
{
	vec2 texCoord = gl_TexCoord[0].xy;
	
	// Add the random texture coordinate offset to the
	// texture coordinate.  The texture is set to perform wrapping
	// so texture coordinates outside the range will be valid.
	vec2 rTexCoord = texCoord + random_texCoordOffset;
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
           
    float i = floor(position.x);
    float j = floor(position.z);
    float k = floor(position.y);
    
    if((i < nx) && (j < nz) && (k < ny) && (i >= 0) && (j >= 0) && (k >=0)){
   

    //This is the initial lookup into the 2D texture that holds the wind field.
 	vec2 index;
   	index.s = j + mod((int)k,numInRow)*nz;
        index.t = i + floor((int)k/numInRow)*nx;
	
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