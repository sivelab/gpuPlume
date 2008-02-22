uniform int nx;
uniform int ny;
uniform int nz;
uniform int numInRow;

varying vec4 pcolor;
varying vec4 particle_pos;

varying vec4 particle_epos;
varying vec4 lpos;

void main(void)
{
  vec4 cr = vec4(1.0, 1.0, 1.0, 1.0);
  vec4 cl = vec4(1.0,1.0,1.0,1.0);
  vec4 ca = vec4(0.2,0.2,0.2,1.0);

  float g = 0.67;
  float pi = 3.141592654;

  vec3 w = normalize(vec3(-particle_epos));
  vec3 wp = normalize(vec3(particle_epos.xyz - lpos.xyz));
  float costheta = dot(w, wp);
  float phase = 1.0 / (4.0*pi) * (1.0 - g*g) / pow(1.0 + g*g - 2.0 * g * costheta, 1.5);

  vec4 color = phase * cl; // + ca;

  // simple radial expansion for alpha
  // float radius = sqrt( ((0.5-gl_TexCoord[0].s) * (0.5-gl_TexCoord[0].s)) + ((0.5-gl_TexCoord[0].t) * (0.5-gl_TexCoord[0].t)) );
  // float alpha = clamp(0.5 - radius, 0.0, 1.0);

  // Gaussian based expansion for alpha
  // recenter x and y (aka s, t) about origin
  float x = gl_TexCoord[0].s - 0.5;
  float y = gl_TexCoord[0].t - 0.5;

  float a = 20.0;
  float xWidth = 0.5;
  float yWidth = 0.5;
  float expX = exp(-a * xWidth*xWidth);
  float expY = exp(-a * yWidth*yWidth);

  float alpha = max(0.0, exp(-a * x * x) - expX) * max(0.0, exp(-a * y * y) - expY);

  gl_FragColor = vec4(color.rgb, alpha);
}
