#version 330 core
out vec4 fragColor;
in vec2 TexCoord;
uniform sampler2D tempTexture;

void main(){
    float temp = texture(tempTexture, TexCoord).r;
    
    float normTemp = (temp - 10.0) / 20.0; 
    vec3 coldColor = vec3(0.0, 0.0, 1.0);
    vec3 hotColor = vec3(1.0, 0.0, 0.0);
    
    fragColor = vec4(mix(coldColor, hotColor, clamp(normTemp, 0.0, 1.0)), 1.0);
}