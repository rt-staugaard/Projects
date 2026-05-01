#include <iostream>
#include <vector>
#include <cstdlib>
#include <thread>
#include <gtc/type_ptr.hpp>
#include "modules/shader.hpp"
#include "modules/camera.hpp"
#include "modules/renderer.hpp"

int rand_int(int max){
    int result = max + 1;
    while (result > max){
        result = std::rand() % max;
    }
    return result;
}

struct System{
    std::vector<float> masses, xs, ys, zs, vx, vy, vz;
    
    System(int n) : masses(n), xs(n), ys(n), zs(n), vx(n, 0), vy(n, 0), vz(n, 0) {}
};

void propagate_states(int state_nr, System &old_system, System &new_system, float t = 0.1){
    float m = old_system.masses[state_nr];
    float x = old_system.xs[state_nr];
    float y = old_system.ys[state_nr];
    float z = old_system.zs[state_nr];
    float vx = old_system.vx[state_nr];
    float vy = old_system.vy[state_nr];
    float vz = old_system.vz[state_nr];

    float acc_x = 0;
    float acc_y = 0;
    float acc_z = 0;
    float epsilon = 1e-9;
    for (int i = 0; i < old_system.xs.size(); ++i){
        if (i != state_nr){
            float dx =  x - old_system.xs[i];
            float dy =  y - old_system.ys[i];
            float dz =  z - old_system.zs[i];
            float R = std::sqrt(dx * dx + dy * dy + dz * dz);

            acc_x +=  m * old_system.masses[i] * dx / (R * R * R + epsilon);
            acc_y +=  m * old_system.masses[i] * dy / (R * R * R + epsilon); 
            acc_z +=  m * old_system.masses[i] * dz / (R * R * R + epsilon);  
        }
    }
    new_system.vx[state_nr] += t * acc_x;
    new_system.vy[state_nr] += t * acc_y;
    new_system.vz[state_nr] += t * acc_z;

    new_system.xs[state_nr] += t * new_system.vx[state_nr];
    new_system.ys[state_nr] += t * new_system.vy[state_nr];
    new_system.zs[state_nr] += t * new_system.vz[state_nr];
}

void manual_parallel_propagate(System &old_system, System &new_system) {
    int n = old_system.xs.size();
    int num_threads = std::thread::hardware_concurrency(); 
    int chunk_size = n / num_threads;

    std::vector<std::thread> workers;

    for (int layer = 0; layer < num_threads; ++layer) {
        int start = layer * chunk_size;
        int end = (layer == num_threads - 1) ? n : (num_threads + 1) * chunk_size;

        workers.emplace_back([start, end, &old_system, &new_system](){
            for (int i = start; i < end; ++i){
                propagate_states(i, old_system, new_system);
            }
        });
    }

    for (auto &t : workers){
        if(t.joinable()){
            t.join();
        }
    }
    std::swap(old_system.masses, new_system.masses);
    std::swap(old_system.xs, new_system.xs);
    std::swap(old_system.ys, new_system.ys);
    std::swap(old_system.zs, new_system.zs);
    std::swap(old_system.vx, new_system.vx);
    std::swap(old_system.vy, new_system.vy);
    std::swap(old_system.vz, new_system.vz);
}

void make_sphere(std::vector<float> &sphereVertices, std::vector<unsigned int> &indices, unsigned int sectorCount = 30, unsigned int stackCount = 30){
    float sectorStep = 2 * M_PI / sectorCount;
    float stackStep = M_PI / stackCount;    
    float sectorAngle, stackAngle;
    int radius = 1;
        
    for (int i = 0; i <= stackCount; ++i){
       stackAngle = M_PI / 2 - i * stackStep;

       float xy = radius * cosf(stackAngle);

        for (int j = 0; j <= sectorCount; ++j){
            sectorAngle = j * sectorStep;

            sphereVertices.insert(sphereVertices.end(), { xy * cosf(sectorAngle), xy * sinf(sectorAngle), radius * sinf(stackAngle)});
        }
    }

    for (int i = 0; i < stackCount; ++i){
       int k1 = i * (sectorCount + 1);
       int k2 = k1 + sectorCount + 1;

        for (int j = 0; j < sectorCount; ++j, ++k1, ++k2){
            indices.push_back(k1);
            indices.push_back(k2);
            indices.push_back(k1 + 1);

            indices.push_back(k1 + 1);
            indices.push_back(k2);
            indices.push_back(k2 + 1);
        }
    }
}


std::vector<unsigned int> setup_buffers() {
    unsigned int VBO, VAO, EBO;

    std::vector<float> sphereVertices;
    std::vector<unsigned int> indices;
    make_sphere(sphereVertices, indices);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sphereVertices.size() * sizeof(float), sphereVertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(VAO);

    return indices;
}    

void drawObject(Shader &shader, glm::mat4 &model, unsigned int modelLoc, std::vector<unsigned int> &indices) {
    shader.use();
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
}    


int main(int argc, char* argv[]){
    std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
    int num_objects = 100;

    float radius = 1;
    std::vector<float> masses;
    std::vector<float> x_values;
    std::vector<float> y_values;
    std::vector<float> z_values;

    std::vector<float> vx_values(num_objects, 0);
    std::vector<float> vy_values(num_objects, 0);
    std::vector<float> vz_values(num_objects, 0);

    int grid_size = 100;

    for (int i = 0; i < num_objects; ++i){
        masses.push_back(10);
        x_values.push_back(rand_int(grid_size));
        y_values.push_back(rand_int(grid_size));
        z_values.push_back(rand_int(grid_size));
    }

    System init_system(num_objects);
    init_system.masses = masses;
    init_system.xs = x_values;
    init_system.ys = y_values;
    init_system.zs = z_values;
    //init_system.vx = vx_values;
    //init_system.vy = vy_values;
    //init_system.vz = vz_values;

    System new_system = init_system;

    GLFWwindow *window;
    if (!glfwInit()){
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(600, 600, "Simulation window", NULL, NULL);
    if (!window){
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glm::vec3 cameraPostion{0,0,0};
    glm::vec3 cameraFront{0,0,-4};
    glm::vec3 cameraUp{0,1,0};
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 600.0f / 600.0f, 0.1f, 100.0f);
    Camera camera(cameraPostion, cameraFront, cameraUp);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);  
    glfwSetWindowUserPointer(window, &camera);
    glfwSetCursorPosCallback(window, mouse_callback_bridge);
    glfwSetScrollCallback(window, scroll_callback_bridge); 

    std::filesystem::path circleVertDir = exe_path / "shaders" / "circle.vert";
    std::filesystem::path circleFragDir = exe_path / "shaders" / "circle.frag";
    auto circle_shader = Shader(circleVertDir.string(), circleFragDir.string());
    auto modelLoc = circle_shader.getID();
    std::vector<unsigned int> indices = setup_buffers();

    float lastFrame;     
    while (!glfwWindowShouldClose(window)){
        float currentFrame = glfwGetTime();
        float deltaTime = currentFrame - lastFrame;
        processInput(window, camera, deltaTime);
        lastFrame = currentFrame;  
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        manual_parallel_propagate(init_system, new_system);
        for (int i = 0; i < num_objects; ++i){
            glm::mat4 object_position = glm::translate(glm::mat4(1.0f), glm::vec3(init_system.xs[i], init_system.ys[i], init_system.zs[i]));
            drawObject(circle_shader, object_position, modelLoc, indices);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
};
