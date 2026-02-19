#version 330 core
layout (location = 0) in vec4 aPos;

uniform mat4 model;

void main() {
    vec3 spatial_postion = vec3(aPos.y, aPos.z, aPos.w);
    gl_Position = model * vec4(spatial_postion, 1.0);
}