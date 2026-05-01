#pragma once
#include <GLFW/glfw3.h>
#include <cmath>
#include "../config.hpp"
#include "InteractionPacket.hpp"

void processInput(GLFWwindow *window, Interaction &inter, int fbWidth, int fbHeight, int winW, int winH){

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        float pixelX = (float)fbWidth / (float)winW * xpos;
        float pixelY = (float)fbHeight / (float)winH * ypos;

        float normX = pixelX / (float)fbWidth;
        float normY = 1.0f - (pixelY / (float)fbHeight);

        inter.targetX = 0.5f + (normX * (gridLength - 1));
        inter.targetY = 0.5f + (normY * (gridHeight - 1));

        if(!inter.isDragging){
            int xCoord = std::clamp((int)std::round(inter.targetX - 0.5f), 0, gridLength - 1);
            int yCoord = std::clamp((int)std::round(inter.targetY - 0.5f), 0, gridHeight - 1);

            inter.targetIdx = xCoord + yCoord * gridLength;

            if (xCoord == 0 || xCoord == gridLength - 1 || yCoord == 0 || yCoord == gridHeight - 1){
                inter.targetIdx = -1;
            }
        }
    }

    else{
        inter.isDragging = false;
        inter.targetIdx = -1;
    }
}


