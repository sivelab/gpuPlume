uniform sampler2D pointsprite_texunit;
uniform sampler2D pointspritenormal_texunit;
uniform samplerRect visualization_texunit;

uniform int nx;
uniform int ny;
uniform int nz;
uniform int numInRow;

varying vec4 pcolor;
varying vec4 particle_pos;

void main(void)
{
  //vec4 d = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
  vec4 normal = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
  if (normal.w == 0.0)
    {
      discard;
    }
  else if (normal.x == 0.0 && normal.y == 0.0 && normal.z == 1.0)
    {
      discard;
    } 
  else
    {
      vec4 l = normalize(vec4(1.0, -1.0, 1.0, 0.0));
      vec4 n = (normal * 2.0) - vec4(1.0);
      n.w = 0.0;
      vec4 h = normalize(l + vec4(0.0,0.0,1.0,0.0));

      float p = 128.0;
      float cp = 0.75;
      vec4 cr = pcolor;
      vec4 cl = vec4(1.0,1.0,1.0,1.0);
      vec4 ca = vec4(0.2,0.2,0.2,1.0);

      vec4 color = cr * (ca + cl * max(0.0,dot(n,l))) + cp * cl * pow(max(0.0,dot(h,n)), p);  
      color.w = 1.0;

      gl_FragColor = color;
    }
}
