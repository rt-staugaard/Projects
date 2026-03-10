#pragma once
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>

class Camera{
    private:
    const float sensitivity = 0.1f;
    bool firstMouse = true;
    float lastX, lastY;
    float yaw = -90.0;
    float pitch;
    float fov = 45;
    float screenWidth = 600.0f;
    float screenHight = 600.0f;
    float nearDistance = 0.1f;
    float farDistance = 100.0f;

    public:
    glm::vec3 position;
    glm::vec3 front;
    glm::vec3 up;
    glm::mat4 look_at_matrix;
    glm::mat4 projection = glm::perspective(glm::radians(fov),screenWidth / screenHight, nearDistance, farDistance);

    Camera(glm::vec3 &position, glm::vec3 &front, glm::vec3 &up){
        this->position = position;
        this->front = front;
        this->up = up;
        look_at_matrix = glm::lookAt(position, front, up);
    }

    void mouse_movement(double xpos, double ypos){
        if (firstMouse) {
            lastX = xpos;
            lastY = ypos;
            firstMouse = false;
        }

        float xoffset = sensitivity * (xpos - lastX);
        float yoffset = sensitivity * (lastY - ypos);
        lastX = xpos;
        lastY = ypos;

        yaw   += xoffset;
        pitch += yoffset;

        if(pitch > 89.0f){
            pitch = 89.0f;
        }
        if(pitch < -89.0f){
            pitch = -89.0f;
        }

        glm::vec3 view;
        view.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        view.y = sin(glm::radians(pitch));
        view.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        front = glm::normalize(view);
        
        glm::vec3 worldUp = glm::vec3(0.0f, 1.0f, 0.0f);
        glm::vec3 right = glm::normalize(glm::cross(front, worldUp));
        up = glm::normalize(glm::cross(right, front));
        look_at_matrix = glm::lookAt(position, position + front, up);
    }

    void scroll_movement(double xoffset, double yoffset){
        fov -= (float)yoffset;
        if (fov < 1.0f){
            fov = 1.0f;
        }
        if (fov > 45.0f){
            fov = 45.0f; 
        }
        projection = glm::perspective(fov, screenWidth / screenHight, nearDistance, farDistance);
    }

    glm::mat4 transform_Model(glm::mat4 &model){
        glm::mat4 transformedModel = projection * look_at_matrix * model;
        return transformedModel;
    }
};