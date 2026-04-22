#pragma once
#include <glm.hpp>

glm::mat3 get_spherical_jacobian(glm::vec3 p) {
    float x = p.x;
    float y = p.y;
    float z = p.z;
    
    float r2 = x * x + y * y + z * z + 1e-4f;
    float r  = sqrt(r2);
    float rho2 = x * x + y * y + 1e-4f;
    float rho  = sqrt(rho2); 
    
    glm::mat3 J;

    J[0][0] = x / r;
    J[0][1] = y / r;
    J[0][2] = z / r;

    J[1][0] = (x * z) / (r2 * rho);
    J[1][1] = (y * z) / (r2 * rho);
    J[1][2] = -rho / r2;

    J[2][0] = -y / rho2;
    J[2][1] =  x / rho2;
    J[2][2] = 0.0f;

    return glm::transpose(J);
}

glm::vec3 convert_to_spherical(glm::vec3 cartestian_vec){
    float x = cartestian_vec.x;
    float y = cartestian_vec.y;
    float z = cartestian_vec.z;

    float r = std::sqrt(x * x + y * y + z * z);
    float theta = std::acos(cartestian_vec.z / (r + 1e-7f));
    float phi = std::atan2(cartestian_vec.y, cartestian_vec.x);

    return glm::vec3(r, theta, phi);
}