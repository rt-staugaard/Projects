#version 330 core
out vec4 fragColor;
in vec2 TexCoord;
uniform sampler2D Texture;

void main(){
    float pixel = texture(Texture, TexCoord).r;
    
    float intensity = pixel > 0.0 ? 1.0 : 0.0;
    
    fragColor = vec4(vec3(intensity), 1.0);
}