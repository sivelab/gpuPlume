varying vec4 pcolor;

void main(void)
{
	vec4 pos = vec4(gl_Vertex);
	pos.w = 1.0;

	// We can control the point size so in the future, if we want individual particle sizes to 
	// be different, then we will need to pass that information to the vertex shader.
	//
	// In the meantime, we will use this at the moment to scale the particle
	// size with distance so particles farther away are smaller.
	//	vec4 modv_vertex = gl_ModelViewMatrix * pos;
	//	gl_PointSize = 62.0 / -modv_vertex.z;

	pcolor = gl_Color;

	gl_Position = gl_ModelViewProjectionMatrix * pos;
}