uniform sampler2DRect posTexSampler;

uniform int nx;
uniform int ny;
uniform int nz;

uniform int doNorm;

void main(void)
{
  vec2 texCoord = gl_TexCoord[0].st;

  vec4 posColor = texture2DRect(posTexSampler, texCoord); 

  vec3 color;


  if (doNorm)
  {
  float mag = length(posColor);
  posColor = normalize(posColor);
  posColor = abs(posColor);

    color = posColor;
} 
  else 
{
  posColor /= vec4(nx, ny, nz, 1.0); // normalize with domain
  color = posColor.xyz;  
}

  gl_FragColor = vec4(color, 1.0);
}
