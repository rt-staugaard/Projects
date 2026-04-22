#pragma once
#include <GLFW/glfw3.h>
#include "camera.hpp"

void processInput(GLFWwindow *window, Camera &camera, float deltaTime){
    
    const float cameraSpeed = 2.5 * deltaTime;
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
        camera.position += cameraSpeed * camera.front;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
        camera.position -= cameraSpeed * camera.front;
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
        camera.position -= cameraSpeed * glm::normalize(glm::cross(camera.front, camera.up));
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
        camera.position += cameraSpeed * glm::normalize(glm::cross(camera.front, camera.up));
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS){
        camera.position += cameraSpeed * camera.up;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
        camera.position -= cameraSpeed * camera.up;
    }

    float distance = glm::length(camera.position);
    if (distance > 40.0f) {
        camera.position = glm::normalize(camera.position) * 40.0f;
    }
    camera.look_at_matrix = glm::lookAt(camera.position, camera.position + camera.front, camera.up);
}

void mouse_callback_bridge(GLFWwindow* window, double xpos, double ypos) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    
    if (camera) {
        camera->mouse_movement(xpos, ypos);
    }
}

void scroll_callback_bridge(GLFWwindow* window, double xoffset, double yoffset) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    
    if (camera) {
        camera->scroll_movement(xoffset, yoffset);
    }
}




