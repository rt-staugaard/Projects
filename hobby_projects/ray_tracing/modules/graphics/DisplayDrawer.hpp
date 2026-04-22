#pragma once
#include <glad/glad.h>
#include <shader.hpp>
#include <gtc/type_ptr.hpp>

struct DisplayDrawer{    
    const float vertices[18] = {
        -1.0f,  1.0f, 0.0f,  
        -1.0f, -1.0f, 0.0f,  
         1.0f, -1.0f, 0.0f,  

        -1.0f,  1.0f, 0.0f,  
        1.0f, -1.0f, 0.0f,  
        1.0f,  1.0f, 0.0f  
    };

    GLuint VAO, VBO;
    std::unique_ptr<Shader> shader;


    DisplayDrawer(char* argv[], const char* vert, const char* frag){
        std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
        std::filesystem::path vertDir = exe_path / vert;
        std::filesystem::path fragDir = exe_path / frag;
        shader = std::make_unique<Shader>(vertDir.string(), fragDir.string());
        shader->use();

        setup_buffers();
    }

    void setup_buffers() {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, 18 * sizeof(float), NULL, GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void draw(glm::mat3 &model) {
        glUniformMatrix3fv(glGetUniformLocation(shader->getID(), "u_viewMatrix"), 1, GL_FALSE, glm::value_ptr(model));
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glBufferSubData(GL_ARRAY_BUFFER, 0,  18 * sizeof(float), &vertices[0]);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }

    ~DisplayDrawer(){}
};


