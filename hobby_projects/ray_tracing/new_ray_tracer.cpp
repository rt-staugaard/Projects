#include "modules/camera.hpp"
#include "modules/shader.hpp"
#include "modules/rays and stars.hpp"
#include "modules/renderer.hpp"
#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <glad/glad.h>
#include <gtc/type_ptr.hpp>
#include <cmath>

const float screenWidth = 600;
const float screenHeight = 600;
const float screenNear = 0.1;
const float screenFar = 100;

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

void setup_buffers(unsigned int &VAO, unsigned int &VBO) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 18 * sizeof(float), NULL, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
}

void draw(std::unique_ptr<Shader> &shader, const float vertices[18], unsigned int &VAO, unsigned int &VBO, glm::mat3 &model) {
    glUniformMatrix3fv(glGetUniformLocation(shader->getID(), "u_viewMatrix"), 1, GL_FALSE, glm::value_ptr(model));
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    glBufferSubData(GL_ARRAY_BUFFER, 0,  18 * sizeof(float), &vertices[0]);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}


int main(int argc, char* argv[]){

GLFWwindow* window;

    if (!glfwInit()){
        return -1;
    }
        
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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

    std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
    std::filesystem::path rayVertDir = exe_path / "shaders" / "light.vert";
    std::filesystem::path rayFragDir = exe_path / "shaders" / "light.frag";
    auto lightShader = std::make_unique<Shader>(rayVertDir.string(), rayFragDir.string());

    lightShader->use();
    auto shaderID = lightShader->getID();
    

    int pixelWidth, pixelHeight;
    glfwGetFramebufferSize(window, &pixelWidth, &pixelHeight);
    glViewport(0, 0, pixelWidth, pixelHeight);
    glUniform2f(glGetUniformLocation(shaderID, "u_resolution"), (float)pixelWidth, (float)pixelHeight);
   
    glUniform1f(glGetUniformLocation(shaderID, "mass_BH"), 1.0f);

    glm::vec3 cameraPosition{5.0f, 1.01f, 30.0f};
    glm::vec3 cameraFront{0.001f, 0.01f, -1};
    glm::vec3 cameraUp{0.01f, 1, 0.01f};

    Camera camera(cameraPosition, cameraFront, cameraUp, screenWidth, screenHeight, screenNear, screenFar);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);  
    glfwSetWindowUserPointer(window, &camera);
    glfwSetCursorPosCallback(window, mouse_callback_bridge);
    glfwSetScrollCallback(window, scroll_callback_bridge); 

    glm::vec3 sourcePosition{-10.0f, -15.0f, -60.0f};
    float sourceRadius = 10.0;
    glUniform1f(glGetUniformLocation(shaderID, "sourceRadius"), sourceRadius);

    const float vertices[] = {
        -1.0f,  1.0f, 0.0f,  
        -1.0f, -1.0f, 0.0f,  
         1.0f, -1.0f, 0.0f,  

        -1.0f,  1.0f, 0.0f,  
        1.0f, -1.0f, 0.0f,  
        1.0f,  1.0f, 0.0f  
    };

    unsigned int VAO, VBO;
    setup_buffers(VAO, VBO);

    glm::mat4 identity(1.0f);
    GLint cameraLoc =glGetUniformLocation(shaderID, "u_cameraPos");
    GLint sourceLoc =glGetUniformLocation(shaderID, "u_sourcePos");

    double lastTime = glfwGetTime();
    int nbFrames = 0;
    while (!glfwWindowShouldClose(window)){
        double currentTime = glfwGetTime();
        nbFrames++;

        float deltaTime = currentTime - lastTime;
        processInput(window, camera, deltaTime);


        if (deltaTime >= 1.0) { 
            double msPerFrame = 1000.0 / double(nbFrames);
            std::cout << msPerFrame << " ms/frame (" << nbFrames << " FPS)" << std::endl;
            nbFrames = 0;
            lastTime = currentTime; 
        }    
        
        glClear(GL_COLOR_BUFFER_BIT);
        glm::mat4 model = camera.look_at_matrix;
        glm::mat3 viewRot = glm::inverse(glm::mat3(model));
        glm::vec3 sourceSphPos = convert_to_spherical(viewRot * sourcePosition);

        glUniform3fv(cameraLoc, 1, &camera.position[0]);
        glUniform3fv(sourceLoc, 1, &sourceSphPos[0]);
        draw(lightShader, vertices, VAO, VBO, viewRot);

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (currentTime > 60.0) break;
    }

    glfwTerminate();
    return 0;
}