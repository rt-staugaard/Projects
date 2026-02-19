#version 330 core
uniform vec3 objectColor;
out vec4 FragColor;

void main(){
    if (distance(gl_PointCoord, vec2(0.5, 0.5)) > 0.5) {
        discard;
    }
    FragColor = vec4(objectColor, 1.0);
}