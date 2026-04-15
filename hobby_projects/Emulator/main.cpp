#include "modules/GLFW_setup.hpp"
#include "modules/Graphics.hpp"
#include "modules/CHIP-8.hpp"
#include <chrono>
#include <thread>

int main(int argc, char* argv[]){

    const int width = 640;
    const int height = 320;
    const char* title = "CHIP-8 EMULATOR";

    std::unique_ptr<Window> AppWindow = std::make_unique<Window>(width, height, title);
    CHIP8 Emulator = CHIP8();
    std::unique_ptr<DisplayDrawer> drawer = std::make_unique<DisplayDrawer>(argv, &Emulator.Display[0]);

    Emulator.LoadProgram("../roms/spaceinvaders.ch8");

    const auto targetCycleTime = std::chrono::microseconds(1000000 / 500);
    auto lastTimerUpdate = std::chrono::high_resolution_clock::now();

    while((!glfwWindowShouldClose(AppWindow->instance)) && Emulator.running){
        auto startTime = std::chrono::high_resolution_clock::now();

        glfwPollEvents();
        AppWindow->processInput(&Emulator.KeyPad[0]);

        uint16_t nextOpcode = Emulator.Fetch();
        Emulator.Decode(nextOpcode);

        if (Emulator.drawFlag){
            glClear(GL_COLOR_BUFFER_BIT);
            drawer->draw(&Emulator.Display[0]);
            glfwSwapBuffers(AppWindow->instance);
            Emulator.drawFlag = false;
        }

        if (std::chrono::duration<float>(startTime - lastTimerUpdate).count() >= (1.0f / 60.0f)) { 
            Emulator.UpdateTimers();
            lastTimerUpdate = std::chrono::high_resolution_clock::now();
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);


        if (elapsedTime < targetCycleTime) {
            std::this_thread::sleep_for(targetCycleTime - elapsedTime);
        }
    }
    return 1;
}