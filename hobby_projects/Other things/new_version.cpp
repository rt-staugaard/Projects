#include "modules/camera.hpp"
#include "modules/shader.hpp"
#include "modules/rays and stars.hpp"
#include "modules/processInput.hpp"
#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <glad/glad.h>
#include <gtc/type_ptr.hpp>
#include <cmath>

const float screenWidth = 600;
const float screenHeight = 600;
const float screenNear = 0.1;
const float screenFar = 100;

glm::vec3 convert_to_spherical(glm::vec3 cartestian_vec){
    float x = cartestian_vec.x;
    float y = cartestian_vec.y;
    float z = cartestian_vec.z;

    float r = std::sqrt(x * x + y * y + z * z);
    float theta = std::acos(cartestian_vec.z / (r + 1e-7f));
    float phi = std::atan2(cartestian_vec.y, cartestian_vec.x);

    return glm::vec3(r, theta, phi);
}

void setup_buffers(unsigned int &VAO, unsigned int &VBO, const float vertices[]) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 18 * sizeof(float), vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
}

unsigned int createPhysicsTexture(int width, int height) {
    unsigned int tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
    return tex;
}

void setupPhysicsFBO(GLuint &fbo, GLuint posTex, GLuint velTex) {
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, posTex, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, velTex, 0);

    GLenum drawBuffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
    glDrawBuffers(2, drawBuffers);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cout << "Error: Framebuffer is not complete!" << std::endl;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void resetSimulation(unsigned int FBO, Shader& initShader, unsigned int posTex, unsigned int velTex, Camera& camera, int w, int h, unsigned int VAO) {
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glViewport(0, 0, w, h);
    
    initShader.use();
    
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, posTex, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, velTex, 0);

    glm::mat3 viewRot = glm::inverse(glm::mat3(camera.look_at_matrix));
    glUniformMatrix3fv(glGetUniformLocation(initShader.getID(), "u_viewMatrix"), 1, GL_FALSE, glm::value_ptr(viewRot));
    glUniform3fv(glGetUniformLocation(initShader.getID(), "u_cameraPos"), 1, &camera.position[0]);
    glUniform2f(glGetUniformLocation(initShader.getID(), "u_resolution"), (float)w, (float)h);

    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
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
    std::filesystem::path rayFragDir = exe_path / "shaders" / "new_light.frag";
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

    glm::vec3 sourcePosition{-10.0f, -25.0f, -100.0f};
    float sourceRadius = 10.0;
    glUniform1f(glGetUniformLocation(shaderID, "sourceRadius"), sourceRadius);
    glm::vec3 SphSourcePosition = convert_to_spherical(sourcePosition);
    glUniform3fv(glGetUniformLocation(shaderID, "u_sourcePos"), 1, &SphSourcePosition[0]);

    const float vertices[] = {
        -1.0f,  1.0f, 0.0f,  
        -1.0f, -1.0f, 0.0f,  
         1.0f, -1.0f, 0.0f,  

        -1.0f,  1.0f, 0.0f,  
        1.0f, -1.0f, 0.0f,  
        1.0f,  1.0f, 0.0f  
    };

    unsigned int VAO, VBO;
    setup_buffers(VAO, VBO, vertices);

    std::filesystem::path rayInitDir = exe_path / "initRays.frag";
    std::filesystem::path rayUpdateDir = exe_path / "updateRays.frag";
    auto initShader = std::make_unique<Shader>(rayVertDir.string(), rayInitDir.string(), true);
    auto updateShader = std::make_unique<Shader>(rayVertDir.string(), rayUpdateDir.string(), true);

    glm::mat4 identity(1.0f);
    GLint cameraLoc = glGetUniformLocation(shaderID, "u_cameraPos");

    double lastTime = glfwGetTime();
    int nbFrames = 0;

    unsigned int posTex_A = createPhysicsTexture(pixelWidth, pixelHeight);
    unsigned int velTex_A = createPhysicsTexture(pixelWidth, pixelHeight);
    unsigned int posTex_B = createPhysicsTexture(pixelWidth, pixelHeight);
    unsigned int velTex_B = createPhysicsTexture(pixelWidth, pixelHeight);

    unsigned int FBO;
    setupPhysicsFBO(FBO, posTex_A, velTex_A);
    resetSimulation(FBO, *initShader, posTex_A, velTex_A, camera, pixelWidth, pixelHeight, VAO);

    glm::vec3 lastCamPos = camera.position;

    while (!glfwWindowShouldClose(window)){
        double currentTime = glfwGetTime();
        nbFrames++;

        float deltaTime = currentTime - lastTime;
        processInput(window, camera, deltaTime);
        if (glm::distance(camera.position, lastCamPos) > 0.01f) {
            resetSimulation(FBO, *initShader, posTex_A, velTex_A, camera, pixelWidth, pixelHeight, VAO);
            lastCamPos = camera.position;
        }

        if (deltaTime >= 1.0) { 
            double msPerFrame = 1000.0 / double(nbFrames);
            std::cout << msPerFrame << " ms/frame (" << nbFrames << " FPS)" << std::endl;
            nbFrames = 0;
            lastTime = currentTime; 
        }    
        
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        glViewport(0, 0, pixelWidth, pixelHeight);
        updateShader->use();
        glUniform1f(glGetUniformLocation(updateShader->getID(), "h"), 0.1f); 
        glUniform1f(glGetUniformLocation(updateShader->getID(), "mass_BH"), 1.0f);

        for (int i = 0; i < 128; i++) {
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, posTex_A);
            glUniform1i(glGetUniformLocation(updateShader->getID(), "u_lastPos"), 0);

            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_2D, velTex_A);
            glUniform1i(glGetUniformLocation(updateShader->getID(), "u_lastVel"), 1);

            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, posTex_B, 0);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, velTex_B, 0);

            glBindVertexArray(VAO);
            glDrawArrays(GL_TRIANGLES, 0, 6);

            std::swap(posTex_A, posTex_B);
            std::swap(velTex_A, velTex_B);
        }

    glBindFramebuffer(GL_FRAMEBUFFER, 0); 
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, pixelWidth, pixelHeight);
    
    lightShader->use();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, posTex_A); 
    glUniform1i(glGetUniformLocation(lightShader->getID(), "u_posTex"), 0);    
    glUniform2f(glGetUniformLocation(shaderID, "u_resolution"), (float)pixelWidth, (float)pixelHeight);
    glUniform1f(glGetUniformLocation(shaderID, "sourceRadius"), sourceRadius);
    glUniform3fv(glGetUniformLocation(shaderID, "u_sourcePos"), 1, &SphSourcePosition[0]);
    glUniform1f(glGetUniformLocation(shaderID, "mass_BH"), 1.0f);
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    glfwSwapBuffers(window);
    glfwPollEvents();

        if (currentTime > 60.0) break;
    }

    glfwTerminate();
    return 0;
}