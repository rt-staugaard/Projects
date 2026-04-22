#include "../../modules/camera/camera.hpp"
#include "../../modules/graphics/GLFW_setup.hpp"
#include "../../modules/graphics/DisplayDrawer.hpp"
#include "../../modules/physics/math_functions.hpp"
#include "../../modules/physics/physics.hpp"
#include "../../modules/physics/rays_and_stars.hpp"

float screenHeight = 600;
float screenWidth = 600;
float screenFar = 100;
float screenNear = 0.1;

const char vertDir[] = "../../shaders/ray.vert";
const char fragDir[] = "../../shaders/ray.frag";

// CAMERA
glm::vec3 initPosition{0,0,30};
glm::vec3 initFront{0,0,-1};
glm::vec3 initUp{0,1,0};

// BLACK HOLE
glm::vec4 BH_pos{0.0f, 0.0f, 0.0f, 0.0f};
glm::vec4 BH_vec{0.0f, 0.0f, 0.0f, 0.0f};

// STAR
glm::vec4 Star_pos{0.0f, -5.0f, 10.0f, -50.0f};
glm::vec4 Star_vel{0.0f, 0.0f, 0.0f, 0.0f};

int main(int argc, char* argv[]){

    std::unique_ptr<Window> window = std::make_unique<Window>(screenWidth, screenHeight, "Simulation Window");
    std::unique_ptr<DisplayDrawer> drawer = std::make_unique<DisplayDrawer>(argv, vertDir, fragDir);
    std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();


    // Setting up camera

    Camera camera(initPosition, initFront, initUp);
    window->bindCamera(camera);

    // Making Spherical Central Object 
    State blackhole_state = State(BH_pos, BH_vec);
    std::shared_ptr<Schwarzchild> blackhole = std::make_shared<Schwarzchild>(1, blackhole_state);

    float blackhole_radius = blackhole->radius;
    std::filesystem::path circleVertDir = exe_path / ".." / ".." / "shaders" / "circle.vert";
    std::filesystem::path circleFragDir = exe_path / ".." / ".." / "shaders" / "circle.frag";
    auto circle_shader = std::make_shared<Shader>(circleVertDir.string(), circleFragDir.string());
    blackhole->set_Shader(circle_shader);

    // Making Light Source
    uint sectorCount = 30;
    uint stackCount = 30;
    float mass = 0.01;
    float star_radius = 2;
    State star_state = State(Star_pos, Star_vel);
    Star star = Star(star_state, mass, star_radius);
    star.set_Shader(circle_shader);

    // Making Ray Bundle
    int number_of_rays = 10;
    std::vector<std::unique_ptr<Ray>> rayBundle;
    rayBundle.reserve(number_of_rays);

    std::filesystem::path rayVertDir = exe_path / ".." / ".." / "shaders" / "line.vert";
    std::filesystem::path rayFragDir = exe_path / ".." / ".." / "shaders" / "line.frag";
    auto ray_shader = std::make_shared<Shader>(rayVertDir.string(), rayFragDir.string());

    for (int i = 0; i < number_of_rays; ++i){
        for(int j = 0; j < number_of_rays; ++j){
            State initial_state(glm::vec4{0.0f, -10.0f, 2*i-10, 2*j-10}, glm::vec4{1.0f, 1, 0.0f, 0.0f});
            rayBundle.push_back(std::make_unique<Ray>(blackhole, initial_state));
            rayBundle.back()->set_Shader(ray_shader); 
        }
    }

    glEnable(GL_DEPTH_TEST);
    float lastFrame = 0;
    glm::mat4 ray_start = glm::mat4(1.0f);
    glm::mat4 blackhole_position = glm::mat4(1.0f);
    glm::mat4 star_position = glm::translate(glm::mat4(1.0f), glm::vec3{star_state.pos.y, star_state.pos.z, star_state.pos.w});

    while (!glfwWindowShouldClose(window->instance)){

        float currentFrame = glfwGetTime();
        float deltaTime = currentFrame - lastFrame;
        window->ManageInput(camera, deltaTime);
        lastFrame = currentFrame;  

        glm::mat4 blackhole_model = camera.transform_Model(blackhole_position);
        glm::mat4 ray_model = camera.transform_Model(ray_start);
        glm::mat4 star_model = camera.transform_Model(star_position);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Exportation of graphics to GPU
        //blackhole->drawObject(blackhole_model);
        star.drawObject(star_model);
        for (auto& ri : rayBundle) {
            ri->propagate_ray(ray_model);
        }

        glfwSwapBuffers(window->instance);
        glfwPollEvents();
    }
    return 0;
}