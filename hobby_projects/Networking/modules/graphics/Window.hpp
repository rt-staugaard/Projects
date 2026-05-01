#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>

struct Window{
    GLFWwindow* instance;
    double xpos, ypos;
    int idx;

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

    ~Window(){
        if (instance) {
            glfwDestroyWindow(instance);
        }
        glfwTerminate();
    }

};
