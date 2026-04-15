#version 330 core
out vec4 outPos;
out vec4 outVel;

uniform vec3 u_cameraPos;
uniform mat3 u_viewMatrix;
uniform vec2 u_resolution;

void Spherical(inout vec4 pos, inout vec4 vel){
    float r = length(pos.yzw) + 1e-4;
    float theta = acos(clamp(pos.w / r, -1.0, 1.0)); 
    float phi = atan(pos.z, pos.y);
    
    float dr = (pos.y * vel.y + pos.z * vel.z + pos.w * vel.w)/(r);
    float dtheta = cos(theta) * cos(phi) / r * vel.y + cos(theta) * sin(phi) / r * vel.z - sin(theta)/r * vel.w;
    float dphi = -pos.z / (pos.y * pos.y + pos.z * pos.z) * vel.y + pos.y /(pos.y * pos.y + pos.z * pos.z) * vel.z;

    pos = vec4(0.0, r, theta, phi);
    vel = vec4(1.0, dr, dtheta, dphi);
}

void main(){
    vec2 uv = gl_FragCoord.xy / u_resolution;
    vec2 scoords = uv * 2.0 - 1.0;
    scoords.x *= (u_resolution.x / u_resolution.y);

    vec3 rayDir = normalize(vec3(scoords.x, scoords.y, -1.0)); 
    vec3 worldDir = u_viewMatrix * rayDir;

    outPos = vec4(0.0, u_cameraPos);
    outVel = vec4(1.0, worldDir);

    Spherical(outPos, outVel);
}
