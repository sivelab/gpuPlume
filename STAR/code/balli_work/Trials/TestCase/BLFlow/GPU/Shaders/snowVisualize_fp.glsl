// Snow Visualization Shader
// Pete Willemsen <willemsn@d.umn.edu>

uniform int nx;
uniform int ny;
uniform int nz;
uniform int numInRow;

varying vec4 pcolor;
varying vec4 particle_pos;

// Eye space values of the particle and light positions and the
// lights' directions.
varying vec4 particle_epos;
varying vec4 lpos1, lpos2, lpos3;
varying vec3 lnorm1;

void main(void)
{
  // Light intensities
  vec4 cl1 = vec4(248.0/255.0, 180.0/255.0, 0.0, 1.0);
  vec4 cl2 = vec4(1.0,1.0,1.0,1.0);
  vec4 cl3 = vec4(248.0/255.0, 151.0/255.0, 0.0, 1.0);

  // Ambient light parameter
  // vec4 ca = vec4(0.5,0.5,0.5,1.0);
  vec4 ca = vec4(0.0,0.0,0.0,1.0);

  // Setting g to 0.0 produces an isotropic scattering, integrated over
  // the sphere (1/4pi).  This can be used as a test case to get
  // transparency and other parameters for ambient light correct.
  // float g = 0.0;
  float g = 0.67;
  float gsquared = g*g;
  float pi = 3.141592654;
  float oneover4pi = 0.07957747;

  // In eye space, so get the vector from particle to eye
  vec3 w = normalize(vec3(-particle_epos));
  // Get the vector from light to particle
  vec3 wp = normalize(vec3(particle_epos.xyz - lpos1.xyz));

  // Using these vectors, compute a visibility value with respect to
  // the light, compute the phase function, and assign a color value
  // to the particle.  Do this for each of the lights in the scene.
  // 
  // Ways to get visibility of 0 or 1 are the "step" and the "sign"
  // function.  The step function returns 0 if 2ndArg < 1stArg,
  // otherwise 1.  A step function with cone of 45 degrees:
  // float vis = step(0.7071, dot(lnorm1, wp));
  float vis = sign(dot(lnorm1, wp));   // use sign of dot prod to determine visibility
  float costheta = dot(w, wp);
  float phase = oneover4pi * ((1.0 - gsquared) / pow(1.0 + gsquared - 2.0 * g * costheta, 1.5));
  vec4 color = vis*phase*cl1 + ca;
  
  wp = normalize(vec3(particle_epos.xyz - lpos2.xyz));
  vis = sign(dot(lnorm1, wp));
  costheta = dot(w, wp);
  phase = oneover4pi * (1.0 - gsquared) / pow(1.0 + gsquared - 2.0 * g * costheta, 1.5);
  color = color + (vis*phase*cl2 + ca);

  wp = normalize(vec3(particle_epos.xyz - lpos3.xyz));
  vis = sign(dot(lnorm1, wp));
  costheta = dot(w, wp);
  phase = oneover4pi * (1.0 - gsquared) / pow(1.0 + gsquared - 2.0 * g * costheta, 1.5);
  color = color + (vis*phase*cl3 + ca);

  // Gaussian based expansion for alpha
  float x = gl_TexCoord[0].s - 0.5;   // recenter x and y (aka s, t) about origin
  float y = gl_TexCoord[0].t - 0.5;
  float a = 0.5;  // height, controls the transparency, currently not greater than 50%
  float c = 0.2;  // width of the Gaussian
  float alpha = a * exp( - ((x*x)/(2.0*c*c)) - ((y*y)/(2.0*c*c)));

  gl_FragColor = vec4(color.rgb, alpha);
}
