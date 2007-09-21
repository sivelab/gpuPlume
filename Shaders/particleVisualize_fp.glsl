uniform sampler2D pointsprite_texunit;
uniform sampler2D pointspritenormal_texunit;
uniform int point_visuals;

varying vec4 pcolor;

void main(void)
{
	if (point_visuals == 1)
	{
		vec4 d = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
		vec4 normal = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
		if (normal.w == 0.0)
			discard;
		else if (d.x == 0.0 && d.y == 0.0 && d.z == 1.0)
		{
		    vec4 color = vec4(0.0, 0.0, 1.0, 1.0);
		    gl_FragColor = color;
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