varying vec4 pcolor;
varying vec4 particle_pos;

void main(void)
{
	vec4 pos = vec4(gl_Vertex);
	pos.w = 1.0;

	// pass the particle position to the fragment shader to perform a
	// lookup into the turbulence field.
	particle_pos = pos;

	pcolor = gl_Color;

	// We can control the point size so in the future, if we want individual particle sizes to 
	// be different, then we will need to pass that information to the vertex shader.
	//
	// In the meantime, we will use this at the moment to scale the particle
	// size with distance so particles farther away are smaller.
	vec4 modv_vertex = gl_ModelViewMatrix * pos;
	gl_PointSize = 100.0 / -modv_vertex.z;

	gl_Position = gl_ModelViewProjectionMatrix * pos;
}