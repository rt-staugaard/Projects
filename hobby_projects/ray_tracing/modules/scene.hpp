#pragma once
#include "physics.cpp"
#include "rays and stars.hpp"
#include "shader.hpp"
#include "camera.hpp"
#include "renderer.hpp"

void make_ray_bundle(int number_of_rays){
    std::vector<std::unique_ptr<Ray>> rayBundle;
    rayBundle.reserve(number_of_rays);

    std::filesystem::path rayVertDir = exe_path / "shaders" / "line.vert";
    std::filesystem::path rayFragDir = exe_path / "shaders" / "line.frag";
    auto ray_shader = std::make_shared<Shader>(rayVertDir.string(), rayFragDir.string());

    for (int i = 0; i < number_of_rays; ++i){
        for(int j = 0; j < number_of_rays; ++j){
            State initial_state(glm::vec4{0.0f, -10.0f, i-4, j-4}, glm::vec4{1.0f, 1, 0.0f, 0.0f});
            rayBundle.push_back(std::make_unique<Ray>(central_object, initial_state));
            rayBundle.back()->set_Shader(ray_shader); 
        }
    }
}

void set_up_scene(std::filesystem::path exe_path){
    int number_of_rays = 10;


}