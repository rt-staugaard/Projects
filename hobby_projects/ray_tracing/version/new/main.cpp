#include "../../modules/camera/camera.hpp"
#include "../../modules/graphics/GLFW_setup.hpp"
#include "../../modules/graphics/DisplayDrawer.hpp"
#include "../../modules/physics/math_functions.hpp"
#include "../../modules/physics/physics.hpp"

const float screenWidth = 600;
const float screenHeight = 600;
const float screenNear = 0.1;
const float screenFar = 100;

const char vertDir[] = "../../shaders/light.vert";
const char fragDir[] = "../../shaders/light.frag";

const float mass_BH = 1.0;

// CAMERA 
glm::vec3 initPosition{5.0f, 1.01f, 30.0f};
glm::vec3 initFront{0.001f, 0.01f, -1};
glm::vec3 initUp{0.01f, 1, 0.01f};

// SOURCE
glm::vec3 sourcePosition{-10.0f, -25.0f, -60.0f};
float sourceRadius = 10.0;


int main(int argc, char* argv[]){

    std::unique_ptr<Window> window = std::make_unique<Window>(screenWidth, screenHeight, "Simulation Window");
    glfwSetCursorPos(window->instance, screenWidth / 2.0, screenHeight / 2.0);

    std::unique_ptr<DisplayDrawer> drawer = std::make_unique<DisplayDrawer>(argv, vertDir, fragDir);

    int pixelWidth, pixelHeight;
    glfwGetFramebufferSize(window->instance, &pixelWidth, &pixelHeight);
    glViewport(0, 0, pixelWidth, pixelHeight);

    GLuint shaderID = drawer->shader->getID();
    glUniform2f(glGetUniformLocation(shaderID, "u_resolution"), (float)pixelWidth, (float)pixelHeight);
    glUniform1f(glGetUniformLocation(shaderID, "mass_BH"), mass_BH);

    Camera camera(initPosition, initFront, initUp, screenWidth, screenHeight, screenNear, screenFar);
    window->bindCamera(camera);

    glUniform1f(glGetUniformLocation(shaderID, "sourceRadius"), sourceRadius);

    glm::mat4 identity(1.0f);
    GLint cameraLoc = glGetUniformLocation(shaderID, "u_cameraPos");
    GLint sourceLoc = glGetUniformLocation(shaderID, "u_sourcePos");

    double lastTime = glfwGetTime();
    int nbFrames = 0;
    while (!glfwWindowShouldClose(window->instance)){
        double currentTime = glfwGetTime();
        nbFrames++;

        float deltaTime = currentTime - lastTime;
        window->ManageInput(camera, deltaTime);

        if (deltaTime >= 1.0) { 
            double msPerFrame = 1000.0 / double(nbFrames);
            std::cout << msPerFrame << " ms/frame (" << nbFrames << " FPS)" << std::endl;
            nbFrames = 0;
            lastTime = currentTime; 
        }    
        
        glClear(GL_COLOR_BUFFER_BIT);
        glm::mat4 view = camera.look_at_matrix;
        glm::mat3 viewRot = glm::inverse(glm::mat3(view));
        glm::vec3 sourceSphPos = convert_to_spherical(viewRot * sourcePosition);

        glUniform3fv(cameraLoc, 1, &camera.position[0]);
        glUniform3fv(sourceLoc, 1, &sourceSphPos[0]);
        drawer->draw(viewRot);

        glfwSwapBuffers(window->instance);
        glfwPollEvents();

        if (currentTime > 60.0) break;
    }

    return 0;
}