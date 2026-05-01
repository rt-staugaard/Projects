#pragma once
#include <glad/glad.h>
#include <shader.hpp>
#include <gtc/type_ptr.hpp>
#include "../config.hpp"

struct DisplayDrawer{    
    GLuint VAO, VBO, EBO;
    std::unique_ptr<Shader> shader;
    std::vector<unsigned int> indices;

    DisplayDrawer(char* argv[], const char* vert, const char* frag){
        std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
        std::filesystem::path vertDir = exe_path / vert;
        std::filesystem::path fragDir = exe_path / frag;
        shader = std::make_unique<Shader>(vertDir.string(), fragDir.string());
        shader->use();

        setupIndices();

        setup_buffers();
    }

    void setupIndices(){
        for (int y = 0; y < gridHeight; y++){
            for (int x = 0; x < gridLength; x++){
                int idx = x + y * gridLength; 

                if (x < gridLength - 1){
                    int right = (x + 1) + y * gridLength;
                    indices.push_back(idx);
                    indices.push_back(right);
                }

                if (y < gridHeight - 1){
                    int up = x + (y + 1) * gridLength;
                    indices.push_back(idx);
                    indices.push_back(up);
                }
            }
        }
    }

    void setup_buffers() {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, 2 * gridHeight * gridLength * sizeof(float), NULL, GL_DYNAMIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void draw(float *Xs, float *Ys, int NumPoints) {
        std::vector<float> vertexData;
        vertexData.reserve(2 * NumPoints);
        for (int i = 0; i < NumPoints; i++){
            float x_mapped = Xs[i] / (gridLength / 2) - 1.0f;
            float y_mapped = Ys[i] / (gridHeight / 2) - 1.0f;

            vertexData.push_back(x_mapped);
            vertexData.push_back(y_mapped);
        }

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glBufferSubData(GL_ARRAY_BUFFER, 0,  vertexData.size() * sizeof(float), vertexData.data());
        glDrawElements(GL_LINES, indices.size(), GL_UNSIGNED_INT, 0);
    }

    ~DisplayDrawer(){}
};


