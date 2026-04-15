#version 330 core


uniform vec2 u_resolution;
uniform sampler2D u_posTex;
uniform sampler2D u_velTex;
uniform vec3 u_sourcePos;
uniform float sourceRadius;
uniform float mass_BH;

out vec4 fragColor;

void main() {
    vec2 uv = gl_FragCoord.xy / u_resolution;
    vec4 pos = texture(u_posTex, uv);
    vec4 vel = texture(u_velTex, uv);

    float r = pos.y;

    if (r < 2.3 * mass_BH) {
        fragColor = vec4(0.0);
        return;
    }
    
    if (r > 55.0) {
        float theta1 = mod(vel.z, 3.14159265); 
        float phi1 = mod(vel.w, 6.28318530);

        float theta2 = mod(u_sourcePos.y, 3.14159265); 
        float phi2 = mod(u_sourcePos.z, 6.28318530);

        float cosGamma = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
        float alpha = atan(sourceRadius / u_sourcePos.x);
        if (cosGamma > cos(alpha)) {    
            fragColor = vec4(1.0); 
            return;
        }
        fragColor = vec4(0.5); 
    }
    else {
        fragColor = vec4(0.2, 0.0, 0.0, 1.0); 
    }
}
