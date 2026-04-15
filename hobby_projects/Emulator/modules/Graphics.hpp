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
    GLuint VAO, VBO, textureID;
    std::unique_ptr<Shader> shader;


    DisplayDrawer(char* argv[], uint8_t* Display){
        std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
        std::filesystem::path vertDir = exe_path / ".." / "shaders" / "chip.vert";
        std::filesystem::path fragDir = exe_path / ".." / "shaders" / "chip.frag";
        shader = std::make_unique<Shader>(vertDir.string(), fragDir.string());
        shader->use();

        setup_buffers();
        glGenTextures(1, &textureID);
        glBindTexture(GL_TEXTURE_2D, textureID);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 64, 32, 0, GL_RED, GL_UNSIGNED_BYTE, Display);
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

    void draw(uint8_t* display) {
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 64, 32, GL_RED, GL_UNSIGNED_BYTE, display);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glBufferSubData(GL_ARRAY_BUFFER, 0,  18 * sizeof(float), &vertices[0]);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }

    ~DisplayDrawer(){}
};









