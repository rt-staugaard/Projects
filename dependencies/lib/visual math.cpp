#include <vector>
#include <stdlib.h>
#include <glad/glad.h>
#include <glm.hpp>

class Grid{
    std::vector<float> xs;
    std::vector<float> ys;

    Grid(float x_start, float y_start, float x_end, float y_end, float spacing_x = 1, float spacing_y = 1){
        int n = static_cast<int>((x_start - x_end) / spacing_x);
        int m = static_cast<int>((y_start - y_end) / spacing_y);
        
        std::vector<float> xs(n);
        for (int i = 0; i < n; ++i){
            xs[i] = x_start + spacing_x * i;
        }

        std::vector<float> ys(m);
        for (int j = 0; j < n; ++j){
            ys[j] = y_start + spacing_y * j;
        }
    }


    void setup_buffer_for_grid(unsigned int &VAO, unsigned int &VBO) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, 18 * sizeof(float), NULL, GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }
};


void draw(std::unique_ptr<Shader> &shader, const float vertices[18], unsigned int &VAO, unsigned int &VBO) {
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    glBufferSubData(GL_ARRAY_BUFFER, 0,  18 * sizeof(float), &vertices[0]);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}