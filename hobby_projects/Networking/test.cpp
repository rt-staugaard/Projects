#include "modules/physics.hpp"
#include "modules/graphics/Window.hpp"
#include "modules/graphics/DisplayDrawer.hpp"
#include "modules/graphics/InputHandling.hpp"

const float screenWidth = 600;
const float screenHeight = 600;

const char vertDir[] = "shaders/grid.vert";
const char fragDir[] = "shaders/grid.frag";

int main(int argc, char* argv[]){

    std::unique_ptr<Window> window = std::make_unique<Window>(screenWidth, screenHeight, "Simulation Window");
    std::unique_ptr<DisplayDrawer> drawer = std::make_unique<DisplayDrawer>(argv, vertDir, fragDir);

    int pixelWidth, pixelHeight;
    glfwGetFramebufferSize(window->instance, &pixelWidth, &pixelHeight);
    glViewport(0, 0, pixelWidth, pixelHeight);

    int winW, winH;
    glfwGetWindowSize(window->instance, &winW, &winH);

    double lastTime = glfwGetTime();
    int nbFrames = 0;

    std::unique_ptr<Grid> current_grid = std::make_unique<Grid>();
    std::unique_ptr<Grid> previous_grid = std::make_unique<Grid>();
    std::unique_ptr<Grid> buffer = std::make_unique<Grid>();

    auto forceX = std::make_unique<float[]>(gridLength * gridHeight);
    auto forceY = std::make_unique<float[]>(gridLength * gridHeight);

    Interaction inter;

    while (!glfwWindowShouldClose(window->instance)){

        double currentTime = glfwGetTime();
        nbFrames++;

        float deltaTime = currentTime - lastTime;

        if (deltaTime >= 1.0) { 
            double msPerFrame = 1000.0 / double(nbFrames);
            std::cout << msPerFrame << " ms/frame (" << nbFrames << " FPS)" << std::endl;
            nbFrames = 0;
            lastTime = currentTime; 
        }

        processInput(window->instance, inter, pixelWidth, pixelHeight,  winW, winH);

        if (!inter.isDragging && inter.targetIdx != -1) {
            current_grid->Xs[inter.targetIdx] = inter.targetX;
            current_grid->Ys[inter.targetIdx] = inter.targetY;

            previous_grid->Xs[inter.targetIdx] = inter.targetX;
            previous_grid->Ys[inter.targetIdx] = inter.targetY;
        }

        updateRows(0, gridHeight - 1, current_grid.get(), previous_grid.get(), buffer.get(), forceX.get(), forceY.get());    
        
        std::unique_ptr<Grid> temp = std::move(previous_grid);
        previous_grid = std::move(current_grid);
        current_grid = std::move(buffer);
        buffer = std::move(temp);

        glClear(GL_COLOR_BUFFER_BIT);

        drawer->draw(&current_grid->Xs[0], &current_grid->Ys[0], gridLength * gridHeight);

        glfwSwapBuffers(window->instance);
        glfwPollEvents();
    }

    return 0;
}