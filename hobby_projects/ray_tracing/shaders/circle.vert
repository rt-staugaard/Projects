#version 330 core
layout (location = 0) in vec4 aPos;

void main(){
    vec3 position = aPos.yzw;
    gl_Position = vec4(position, 1.0); 
    gl_PointSize = 60.0;
}