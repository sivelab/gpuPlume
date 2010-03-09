#version 120

// Incoming color value from the vertex shader.
varying vec4 gl_Color;

// Incoming light intensity from the vertex shader.
varying float lightIntensity;

// The currently bound texture.
uniform sampler2D boundTexture;

// The currently bound texture.
uniform sampler2DShadow shadowMap;

varying vec4 shadowCoord;

void main(void)
{
  float shadow = shadow2DProj(shadowMap, shadowCoord).x;
  
//  vec4 shadowCoordinateWdivide = shadowCoord / shadowCoord.w ;
  
//  Used to lower moir pattern and self-shadowing
//  shadowCoordinateWdivide.z += 0.0005;
  
//  float distanceFromLight = texture2D(shadowMap, shadowCoordinateWdivide.st).z;
  
//  float shadow = 1.0;
  
//  if (shadowCoord.w > 0.0)
//    shadow = distanceFromLight < shadowCoordinateWdivide.z ? 0.5 : 1.0 ;
  
  vec4 textureColor = vec4(texture2D(boundTexture, gl_TexCoord[0].st));
  
  // Test Output
  gl_FragColor = textureColor * shadow * lightIntensity;
  
}
