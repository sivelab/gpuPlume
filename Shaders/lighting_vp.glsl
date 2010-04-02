#version 120

// Incoming color value from OpenGL
attribute vec4 gl_Color;

// Outgoing intensity value for the fragment shader
varying float lightIntensity;

// Outgoing color values for the fragment shader
varying vec4 gl_FrontColor;
varying vec4 gl_BackColor;

// lightPosition is the position of the light which is passed in from OpenGL
uniform vec3 lightPosition;

// light transform matricies
uniform mat4x4 lightModelviewMatrix;
uniform mat4x4 lightProjectionMatrix;
uniform mat4x4 lightScaleAndBiasMatrix;
uniform mat4x4 cameraModelviewMatrix;

varying vec4 shadowCoord;

// Light attributes that are specified by constants
const float specularContribution = 0.1;
const float diffuseContribution = 1.0 - specularContribution;

void main()
{
  
  // Light the scene
  vec3 ecPosition = vec3(gl_ModelViewMatrix * gl_Vertex);
  vec3 tnorm = normalize(gl_NormalMatrix * gl_Normal);
  vec3 lightVec = normalize(lightPosition - ecPosition);
  vec3 reflectVec = reflect(-lightVec, tnorm);
  vec3 viewVec = normalize(-ecPosition);
  float spec = clamp(dot(reflectVec, viewVec), 0.0, 1.0);
  spec = pow(spec, 16.0);
  lightIntensity = diffuseContribution * max(dot(lightVec, tnorm), 0.0) + specularContribution * spec;
  
  lightIntensity *= 1.5;
  
  // Pass along the gl_Color value
  gl_FrontColor = gl_Color * lightIntensity;
  gl_BackColor = gl_Color * lightIntensity;
  
  // Pass along the bound texture
  gl_TexCoord[0] = gl_MultiTexCoord0;
  
  // Pass along the bound texture (shadowMap).
  gl_TexCoord[1] = gl_MultiTexCoord1;
  
  // Do some setup for the shadowing in the fragment shader
  shadowCoord = (lightScaleAndBiasMatrix * lightProjectionMatrix * lightModelviewMatrix) * gl_Vertex;
  
  // Perform the standard transformation
  gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * gl_Vertex;
  
}

