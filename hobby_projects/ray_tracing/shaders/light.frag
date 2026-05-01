#version 330 core

const float G = 1.0;

struct State{
    vec4 pos;
    vec4 vel;
};

uniform vec2 u_resolution;
uniform vec3 u_cameraPos;
uniform vec3 u_sourcePos;
uniform float mass_BH;
uniform float sourceRadius;
uniform mat3 u_viewMatrix;

out vec4 fragColor;

State Spherical(vec4 pos, vec4 vel){
    float r = length(pos.yzw) + 1e-4;
    float theta = acos(clamp(pos.w / r, -1.0, 1.0)); 
    float phi = atan(pos.z, pos.y);
    
    float dr = (pos.y * vel.y + pos.z * vel.z + pos.w * vel.w)/(r);
    float dtheta = cos(theta) * cos(phi) / r * vel.y + cos(theta) * sin(phi) / r * vel.z - sin(theta)/r * vel.w;
    float dphi = -pos.z/(pos.y * pos.y + pos.z * pos.z) * vel.y + pos.y /(pos.y * pos.y + pos.z * pos.z) * vel.z;

    return State(vec4(pos.x, r, theta, phi), vec4(vel.x, dr, dtheta, dphi));
}

State get_acceleration(vec4 pos, vec4 vel) {
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

    State change;
    change.pos = vel; 

    change.vel.x = -(r_s * inv_r * inv_omrs * inv_r) * dt * dr;

    change.vel.y = -(r_s * omrs / (2.0 * r * r)) * dt * dt 
                   + (r_s * inv_r * inv_omrs / 2.0) * dr * dr 
                   + (r * omrs) * dth * dth 
                   + (r * sinT * sinT * omrs) * dph * dph;

    change.vel.z = -(2.0 * inv_r) * dr * dth + (sinT * cosT) * dph * dph;

    change.vel.w = -(2.0 * inv_r) * dr * dph - (2.0 * cotT) * dth * dph;

    return change;
}

void update(inout State s, float h){
    State k = get_acceleration(s.pos, s.vel);
    State total = k;

    k = get_acceleration(s.pos + k.pos * (h * 0.5), s.vel + k.vel * (h * 0.5));
    total.pos += 2.0 * k.pos;
    total.vel += 2.0 * k.vel;

    s.pos += total.pos * h;
    s.vel += total.vel * h;
}

float hash(vec2 p) {
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

float starNoise(vec2 uv) {
    vec2 i = floor(uv);
    vec2 f = fract(uv);
    float a = hash(i);
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));
    vec2 u = f * f * (3.0 - 2.0 * f);
    return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}


void main() {
    vec2 st = gl_FragCoord.xy / u_resolution.xy;
    vec2 uv = st * 2.0 - 1.0;
    uv.x *= (u_resolution.x / u_resolution.y);

    vec3 rayDir = normalize(vec3(uv.x, uv.y, -1.0)); 
    vec3 dir = u_viewMatrix * rayDir;

    vec4 startPos = vec4(0.0, u_cameraPos);
    vec4 startVel = vec4(1, dir);
    State ray = State(startPos, startVel);

    float theta2 = clamp(u_sourcePos.y,0.001, 3.141);
    float phi2 = clamp(u_sourcePos.z, 0.001, 2 * 3.141);

    ray = Spherical(ray.pos, ray.vel);

    vec3 camToBH = normalize(-u_cameraPos);
    float cosAngle = dot(dir, camToBH);

    float R = length(u_cameraPos);
    float shadowRadius = 2.4 * 2.0 * mass_BH; // 2.4 effective boundary
    float sinAlpha = shadowRadius / R;
    float cosAlpha = sqrt(1.0 - sinAlpha * sinAlpha);

    if (cosAngle > cosAlpha) {
        fragColor = vec4(0.0, 0.0, 0.0, 1.0);
        return;
    }

    vec3 pixelColor = vec3(1.0, 0.0, 0.0);

    for(int i = 0; i < 1024; i++) {
        float h_radial = min(0.35, (ray.pos.y - 2.0 * mass_BH) * 0.3);
        float max_angular_vel = max(abs(ray.vel.z), abs(ray.vel.w)) + 1e-4;
        float h_angular = 0.0035 / max_angular_vel;
        float h = min(h_radial, h_angular);
        h = clamp(h, 0.0001, 0.15);

        update(ray, h);
    
        float r = ray.pos.y;

        // Check Event Horizon
        if (r < 2.1 * mass_BH) {
            pixelColor = vec3(0.0);
            break;
        }
    
        if (r > 50) {
            ray.pos.z = mod(ray.pos.z, 3.14159265); 
            ray.pos.w = mod(ray.pos.w, 6.28318530);
            float theta1 = clamp(ray.pos.z, 0.001, 3.141);
            float phi1 = clamp(ray.pos.w, 0.001, 2 * 3.141);

            float cosGamma = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
            float alpha = atan(sourceRadius / u_sourcePos.x);
            if (cosGamma > cos(alpha)) {
                float dist = acos(cosGamma) / alpha; 
                float centerFactor = 1.0 - smoothstep(0.0, 1.0, dist);

                vec2 starUV = vec2(ray.pos.z * 5.0, ray.pos.w * 10.0);
                float noise = starNoise(starUV) * 0.5 + 0.5;

                vec3 baseColor = vec3(0.9, 0.8, 0.8); 
                vec3 hotColor = vec3(1.0, 1.0, 1.0);  
                vec3 finalStar = mix(baseColor, hotColor, noise);

                float limbDarkening = pow(centerFactor, 0.2);
                pixelColor = finalStar * limbDarkening * 1.5;    
                break;
            }
            pixelColor = vec3(0.5); 
            break;
        }
    }
    fragColor = vec4(pixelColor, 1.0);
}
