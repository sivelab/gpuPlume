uniform sampler2DRect positions;
uniform float x;
uniform float y;
varying vec4 pos;

void main(void)
{
  //vec2 texCoord = (pos.x,pos.y);
  vec2 texCoord;
  texCoord.s = pos.x;
  texCoord.t = pos.y;
  vec3 point = vec3(texture2DRect(positions, texCoord));
  gl_FragColor = vec4(point,1.0);

  //gl_FragColor = vec4(x,y,0.0,1.0);
}
