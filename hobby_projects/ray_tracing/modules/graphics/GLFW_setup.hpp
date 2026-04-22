#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "../camera/processInput.hpp"
#include "../camera/camera.hpp"


struct Window{
    GLFWwindow* instance;

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

    void bindCamera(Camera &camera){
        glfwSetWindowUserPointer(instance, &camera);
        glfwSetCursorPosCallback(instance, mouse_callback_bridge);
        glfwSetScrollCallback(instance, scroll_callback_bridge);
        glfwSetInputMode(instance, GLFW_CURSOR, GLFW_CURSOR_DISABLED);   
    }

    void ManageInput(Camera &camera, float deltaTime){
        processInput(this->instance, camera, deltaTime);
    }

    ~Window(){
        if (instance) {
            glfwDestroyWindow(instance);
        }
        glfwTerminate();
    }

};


