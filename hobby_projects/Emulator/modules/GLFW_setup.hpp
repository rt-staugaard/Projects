#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>


struct Window{
    GLFWwindow* instance;

    const int HEX_MAP[16] = {
        GLFW_KEY_X, //0
        GLFW_KEY_1, //1
        GLFW_KEY_2, //2
        GLFW_KEY_3, //3
        GLFW_KEY_Q, //4
        GLFW_KEY_W, //5
        GLFW_KEY_E, //6
        GLFW_KEY_A, //7
        GLFW_KEY_S, //8
        GLFW_KEY_D, //9
        GLFW_KEY_Z, //A
        GLFW_KEY_C, //B
        GLFW_KEY_4, //C
        GLFW_KEY_R, //D
        GLFW_KEY_F, //E
        GLFW_KEY_V  //F
    };

    Window(int width, int height, const char* title){
        if (!glfwInit()){
            std::cerr << "[ERROR] Could not initiate window\n";
        }
        
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        GLFWwindow* instance = glfwCreateWindow(width, height, title, NULL, NULL);

        if (!instance){
            std::cerr << "[ERROR] Window not constructed, terminating GLFW\n";
            glfwTerminate();
        }

        glfwMakeContextCurrent(instance);

        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
            std::cerr << "[ERROR] Failed to initialize GLAD\n";
        }

        this->instance = instance;
    }

    void processInput(uint8_t* emulatorKeypad){
        for (int i = 0; i < 16; i++) {
            if (glfwGetKey(instance, HEX_MAP[i]) == GLFW_PRESS) {
                emulatorKeypad[i] = 1;
            } else {
                emulatorKeypad[i] = 0;
            }
        }
    }

    ~Window(){
        if (instance) {
            glfwDestroyWindow(instance);
        }
        glfwTerminate();
    }

};


