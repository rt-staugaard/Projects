#version 330 core
out vec4 nextPos;
out vec4 nextVel;

uniform sampler2D u_lastPos;
uniform sampler2D u_lastVel;
uniform float h;
uniform float mass_BH;

void get_acceleration(vec4 pos, vec4 vel, inout vec4 acc) {
    float r = max(pos.y, 2.05 * mass_BH); 
    float theta = pos.z;
    float r_s = 2.0 * mass_BH;
    
    float sinT = sin(theta);
    float cosT = cos(theta);
    float cotT = cosT / (sign(sinT) * max(abs(sinT), 1e-4));
    
    float inv_r = 1.0 / r;
    float omrs = 1.0 - r_s * inv_r; 
    float inv_omrs = 1.0 / max(omrs, 1e-5);
    
    float dt = vel.x;
    float dr = vel.y;
    float dth = vel.z;
    float dph = vel.w;

    acc.x = -(r_s * inv_r * inv_omrs * inv_r) * dt * dr;

    acc.y = -(r_s * omrs / (2.0 * r * r)) * dt * dt 
                   + (r_s * inv_r * inv_omrs / 2.0) * dr * dr 
                   + (r * omrs) * dth * dth 
                   + (r * sinT * sinT * omrs) * dph * dph;

    acc.z = -(2.0 * inv_r) * dr * dth + (sinT * cosT) * dph * dph;

    acc.w = -(2.0 * inv_r) * dr * dph - (2.0 * cotT) * dth * dph;
}

void update(inout vec4 pos, inout vec4 vel, float h){
    vec4 startPos = pos;
    vec4 startVel = vel;
    vec4 acc = vec4(0.0);
    get_acceleration(pos, vel, acc);

    vec4 midPos = startPos +  startVel * (h * 0.5);
    vec4 midVel = startVel + acc * (h * 0.5);
    vec4 new_acc = vec4(0.0);
    get_acceleration(midPos, midVel, new_acc);

    pos = startPos + midVel * h;
    vel = startVel + new_acc * h;
}


void main() {
    ivec2 texSize = textureSize(u_lastPos, 0);
    vec2 uv = gl_FragCoord.xy / vec2(texSize);

    vec4 pos = texture(u_lastPos, uv);
    vec4 vel = texture(u_lastVel, uv);
    
    if (pos.y > 2.3 * mass_BH && pos.y < 60.0) {
        update(pos, vel, h);
    }

    nextPos = pos;
    nextVel = vel;
}
