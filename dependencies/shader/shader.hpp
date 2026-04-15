#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <glad/glad.h>

class Shader{
    private:
    uint shaderProgram;

    public:
    Shader() : shaderProgram(0) {}

    Shader(const std::string& vertex_filepath, const std::string fragment_filepath, bool isPhysicsShader = false){
        make_shader_program(vertex_filepath, fragment_filepath, isPhysicsShader);
    }

    uint make_shader_module(const std::string& filepath, unsigned int module_type){
        std::ifstream file(filepath);
        std::stringstream bufferedLines;
        std::string line;

        while(std::getline(file, line)){
           bufferedLines << line << "\n";
        }

        std::string sourceCode = bufferedLines.str();
        const char* shaderSource = sourceCode.c_str();
        bufferedLines.str("");
        file.close();
    
        unsigned int shaderModule = glCreateShader(module_type);
        glShaderSource(shaderModule, 1, &shaderSource, NULL);
        glCompileShader(shaderModule);

        int success;
        glGetShaderiv(shaderModule, GL_COMPILE_STATUS, &success);
        if (!success){
            char errorLog[1024];
            glGetShaderInfoLog(shaderModule, 1024, NULL, errorLog);
            std::cout << "Shader Module compilation error:\n" << errorLog << std::endl;
        }

        return shaderModule;
    }   

    virtual void make_shader_program(const std::string& vertex_filepath, const std::string fragment_filepath, bool isPhysicsShader){
        std::vector<unsigned int> modules;
        modules.push_back(make_shader_module(vertex_filepath, GL_VERTEX_SHADER));
        modules.push_back(make_shader_module(fragment_filepath, GL_FRAGMENT_SHADER));
        
        unsigned int shaderProgram = glCreateProgram();
        for (unsigned int shaderModule : modules){
            glAttachShader(shaderProgram, shaderModule);
        }

        if (isPhysicsShader) {
            glBindFragDataLocation(shaderProgram, 0, "nextPos");
            glBindFragDataLocation(shaderProgram, 1, "nextVel");
        } else {
            glBindFragDataLocation(shaderProgram, 0, "fragColor");
        }
        glLinkProgram(shaderProgram);

        int success;
        glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
        if (!success){
            char errorLog[1024];
            glGetProgramInfoLog(shaderProgram, 1024, NULL, errorLog);
            std::cout << "Shader Program compilation error:\n" << errorLog << std::endl;
        }

        for (unsigned int shaderModule : modules){
            glDeleteShader(shaderModule);
        }
        this->shaderProgram = shaderProgram;
    }

    void use(){
        glUseProgram(shaderProgram);
    }

    void set_shader_program(unsigned int program){
        this->shaderProgram = program;
    }

    uint getID() const { return shaderProgram; }

    ~Shader() {
        glDeleteProgram(shaderProgram);
    }
};

