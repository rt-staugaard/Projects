#version 330 core

uniform float ray_time;
uniform float real_time;

out vec4 FragColor;

void main() {
    float time_delay =  ray_time/(real_time + 0.001);
    FragColor = vec4(1.0, 1.0/time_delay, 1.0/time_delay, 1.0);
}