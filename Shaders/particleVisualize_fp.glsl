uniform sampler2D pointsprite_texunit;
uniform sampler2D pointspritenormal_texunit;
uniform samplerRect visualization_texunit;

uniform int point_visuals;
uniform int nx;
uniform int ny;
uniform int nz;
uniform int numInRow;

varying vec4 pcolor;
varying vec4 particle_pos;

void main(void)
{
  if (point_visuals == 1)
    {
      //vec4 d = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
      vec4 normal = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
      if (normal.w == 0.0)
	discard;

      else if (normal.x == 0.0 && normal.y == 0.0 && normal.z == 1.0)
	{
	  //The floor of the position in 3D space is needed to find the index into
	  //the 2D Texture.
	  int i = int(floor(particle_pos.y));
	  int j = int(floor(particle_pos.x));
	  int k = int(floor(particle_pos.z));
	
	  vec4 cvis = vec4(1.0, 0.0, 0.0, 1.0);

	  gl_FragColor = cvis;
	} 
      else
	{
	  vec4 l = normalize(vec4(1.0, -1.0, 1.0, 0.0));
	  vec4 n = (normal * 2.0) - vec4(1.0);
	  n.w = 0.0;
	  vec4 h = normalize(l + vec4(0.0,0.0,1.0,0.0));

	  float p = 64.0;
	  float cp = 1.0;
	  vec4 cr = pcolor;
	  vec4 cl = vec4(1.0,1.0,1.0,1.0);
	  vec4 ca = vec4(0.4,0.4,0.4,1.0);

	  vec4 color = cr * (ca + cl * max(0.0,dot(n,l))) + cp * cl * pow(max(0.0,dot(h,n)), p);  
	  gl_FragColor = color;
	}
    }
  else 
    {
      // 
      // If using simple point based particles, this is all that needs to be done
      // 
      gl_FragColor = pcolor;
    }
}
