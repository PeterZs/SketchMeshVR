#version 330
#extension GL_ARB_explicit_uniform_location : enable

layout(location = 0) out vec4 fragcolor;   

layout(location = 1) uniform sampler2D tex;
layout(location = 2) uniform vec4 color = vec4(1.0);

in vec2 tex_coord;

void main(void)
{   
   fragcolor = texture(tex, tex_coord)*color;
}




















