#include "modules/graphics/Window.hpp"
#include "modules/graphics/DisplayDrawer.hpp"
#include "modules/graphics/InputHandling.hpp"
#include "modules/network/client.hpp"
#include "modules/config.hpp"
#include <thread>

const float screenWidth = 300;
const float screenHeight = 300;

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

    float* Xs = new float[gridHeight * gridLength];
    float* Ys = new float[gridHeight * gridLength];
    float* XsBuffer = new float[gridHeight * gridLength];
    float* YsBuffer = new float[gridHeight * gridLength];

    Client client;

    client.ConnectToServer("127.0.0.1", 8080);

    std::mutex mtx;
    bool freshData = false;

    std::thread listenerThread(NetworkWorker, 
                           std::ref(client),    
                           std::ref(Xs),        
                           std::ref(Ys),        
                           std::ref(XsBuffer),  
                           std::ref(YsBuffer),  
                           std::ref(mtx),
                           gridHeight * gridLength,       
                           std::ref(freshData)
                           );

    listenerThread.detach();

    Interaction inter;
    while(client.getConnectivity() && !glfwWindowShouldClose(window->instance)){
        glClear(GL_COLOR_BUFFER_BIT);

        if (freshData){
            std::lock_guard<std::mutex> lock(mtx);
            freshData = false;
        }

        drawer->draw(&Xs[0], &Ys[0], gridLength * gridHeight);
        processInput(window->instance, inter, pixelWidth, pixelHeight,  winW, winH);

        if (!inter.isDragging && inter.targetIdx != -1){
            client.sendPacket(inter);
            inter.targetIdx = -1;
        }

        glfwSwapBuffers(window->instance);
        glfwPollEvents();

    }
    delete[] Xs; delete[] Ys; delete[] XsBuffer; delete[] YsBuffer;
    return 0;
}
