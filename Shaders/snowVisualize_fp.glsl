uniform int nx;
uniform int ny;
uniform int nz;
uniform int numInRow;

varying vec4 pcolor;
varying vec4 particle_pos;

varying vec4 particle_epos;
varying vec4 lpos1, lpos2, lpos3;

void main(void)
{
  vec4 cl1 = vec4(1.0,0.0,0.0,1.0);
  vec4 cl2 = vec4(0.0,1.0,0.0,1.0);
  vec4 cl3 = vec4(0.0,0.0,1.0,1.0);
  // vec4 ca = vec4(0.5,0.5,0.5,1.0);
  vec4 ca = vec4(0.0,0.0,0.0,1.0);

  vec3 lnormal = vec3(0.0, 0.0, -1.0);

  // Setting g to 0.0 produces an isotropic scattering, integrated over
  // the sphere (1/4pi).  This can be used as a test case to get
  // transparency and other parameters for ambient light correct.
  // float g = 0.0;
  float g = 0.67;
  float pi = 3.141592654;
  vec3 w = normalize(vec3(-particle_epos));

  vec3 wp = normalize(vec3(particle_epos.xyz - lpos1.xyz));

  // that should be the costheta for the theta between the light's normal and where the particle is relative to light
  // so if this value is greaterthan (45 degrees), so cos(45) or 0.7071 then zero out light...
  float vis = max(0.0, dot(lnormal, wp));
  if (vis > 0.7071) vis = 0.0; else vis = 1.0;

  float costheta = dot(w, wp);
  float phase = 1.0 / (4.0*pi) * (1.0 - g*g) / pow(1.0 + g*g - 2.0 * g * costheta, 1.5);
  vec4 color = phase*cl1 + ca;

  wp = normalize(vec3(particle_epos.xyz - lpos2.xyz));
  vis = max(0.0, dot(lnormal, wp));
  if (vis > 0.7071) vis = 0.0; else vis = 1.0;
  costheta = dot(w, wp);
  phase = 1.0 / (4.0*pi) * (1.0 - g*g) / pow(1.0 + g*g - 2.0 * g * costheta, 1.5);
  color = color + (phase*cl2 + ca);

  wp = normalize(vec3(particle_epos.xyz - lpos3.xyz));
  vis = max(0.0, dot(lnormal, wp));
  if (vis > 0.7071) vis = 0.0; else vis = 1.0;
  costheta = dot(w, wp);
  phase = 1.0 / (4.0*pi) * (1.0 - g*g) / pow(1.0 + g*g - 2.0 * g * costheta, 1.5);
  color = color + (phase*cl3 + ca);

  // simple radial expansion for alpha
  // float radius = sqrt( ((0.5-gl_TexCoord[0].s) * (0.5-gl_TexCoord[0].s)) + ((0.5-gl_TexCoord[0].t) * (0.5-gl_TexCoord[0].t)) );
  // float alpha = clamp(0.5 - radius, 0.0, 1.0);

  // Gaussian based expansion for alpha
  float x = gl_TexCoord[0].s - 0.5;   // recenter x and y (aka s, t) about origin
  float y = gl_TexCoord[0].t - 0.5;
  float a = 0.5;  // height, controls the transparency, currently not greater than 50%
  float c = 0.2;  // width of the Gaussian
  float alpha = a * exp( - ((x*x)/(2.0*c*c)) - ((y*y)/(2.0*c*c)));

  gl_FragColor = vec4(color.rgb, alpha);
}
