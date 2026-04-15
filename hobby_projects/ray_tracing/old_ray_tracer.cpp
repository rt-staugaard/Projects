#include "modules/physics.cpp"
#include "modules/rays and stars.hpp"
#include "modules/shader.hpp"
#include "modules/camera.hpp"
#include "modules/processInput.hpp"

int main(int argc, char* argv[])
{   
    std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
    GLFWwindow* window;

    if (!glfwInit()){
        return -1;
    }
        
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    float screenHeight = 600;
    float screenWidth = 600;
    float screenFar = 100;
    float screenNear = 0.1;

    window = glfwCreateWindow(screenWidth, screenHeight, "Simulation Window", NULL, NULL);
    if (!window){
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // Setting up camera
    glm::vec3 cameraPosition{0,0,30};
    glm::vec3 cameraFront{0,0,-1};
    glm::vec3 cameraUp{0,1,0};
    Camera camera(cameraPosition, cameraFront, cameraUp);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);  
    glfwSetWindowUserPointer(window, &camera);
    glfwSetCursorPosCallback(window, mouse_callback_bridge);
    glfwSetScrollCallback(window, scroll_callback_bridge); 

    // Making Spherical Central Object 
    State blackhole_state = State(glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f}, glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f});
    std::shared_ptr<Schwarzchild> blackhole = std::make_shared<Schwarzchild>(1, blackhole_state);
    float blackhole_radius = blackhole->radius;
    std::filesystem::path circleVertDir = exe_path / "shaders" / "circle.vert";
    std::filesystem::path circleFragDir = exe_path / "shaders" / "circle.frag";
    auto circle_shader = std::make_shared<Shader>(circleVertDir.string(), circleFragDir.string());
    blackhole->set_Shader(circle_shader);

    // Making Light Source
    uint sectorCount = 30;
    uint stackCount = 30;
    float mass = 0.01;
    float star_radius = 2;
    State star_state = State(glm::vec4 {0.0f, -5.0f, 10.0f, -50.0f}, glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f});
    Star star = Star(star_state, mass, star_radius);
    star.set_Shader(circle_shader);

    // Making Ray Bundle
    int number_of_rays = 10;
    std::vector<std::unique_ptr<Ray>> rayBundle;
    rayBundle.reserve(number_of_rays);

    std::filesystem::path rayVertDir = exe_path / "shaders" / "line.vert";
    std::filesystem::path rayFragDir = exe_path / "shaders" / "line.frag";
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

    while (!glfwWindowShouldClose(window)){

        float currentFrame = glfwGetTime();
        float deltaTime = currentFrame - lastFrame;
        processInput(window, camera, deltaTime);
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

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}