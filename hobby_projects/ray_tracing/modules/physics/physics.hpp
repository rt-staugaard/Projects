#pragma once
#include <glm.hpp>
#include <glad/glad.h>
#include <gtc/type_ptr.hpp>
#include "shader.hpp"

struct State{
    glm::vec4 pos;
    glm::vec4 vel;

    State(glm::vec4 position, glm::vec4 velocity){
        this->pos = position;
        this->vel = velocity;
    }

    State(){}

    State operator+(const State& rhs) const {
        return State(pos + rhs.pos, vel + rhs.vel);
    }

    State operator*(const float& k) const {
        return State(k * pos, k * vel);
    }
};

// By default Minikowski
class spacetime{
    private:
    float G = 1.0;

    public:
    glm::mat4 gamma[4];

    virtual glm::mat4 get_metric(State state){
        glm::vec4 pos = state.pos;
        glm::mat4 g(0.0);
        g[0][0] = - pos.x * pos.x;
        g[1][1] =   pos.y * pos.y;
        g[2][2] =   pos.z * pos.z;
        g[3][3] =   pos.w * pos.w;

        return g;
    }

    void insert_christoffel_symbol(float value, int &mu, int &alpha, int &beta){
        gamma[mu][alpha][beta] = value;
        gamma[mu][beta][alpha] = value;
    }

    glm::mat4 get_christoffel_symbol(int &mu){
        glm::mat4 christoffel_symbol = gamma[mu];
        return christoffel_symbol;
    }

    virtual ~spacetime() {}
};

class Schwarzchild : public spacetime {
public:
    private:
    unsigned int VAO, VBO, EBO;
    std::shared_ptr<Shader> shader;
    glm::mat4 model;
    std::vector<float> sphereVertices;
    std::vector<unsigned int> indices;
    unsigned int modelLoc;

    public:
    float mass;
    float radius;

    Schwarzchild(float mass, State state, float radius = 2, unsigned int sectorCount = 30, unsigned int stackCount = 30){
        this->mass = mass;
        this->radius = radius;
        make_sphere(sectorCount, stackCount);
        setup_buffers();
    }

    glm::mat4 get_metric(State state) override {
        glm::vec4 pos = state.pos;
        float x = state.pos.y;
        float y = state.pos.z;
        float z = state.pos.w;

        float R = std::sqrt(x * x + y * y + z * z);
        float epsilon = 1e-8;

        glm::mat4 g(0.0);

        g[0][0] = -(1.0 - 2.0 * mass / R);

        float factor = (2.0 * mass) / (R * R * (R - 2.0 * mass + epsilon));

        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
                float delta_ij = (i == j) ? 1.0 : 0.0;
                g[i][j] = delta_ij + factor * pos[i] * pos[j];
            }
        }
        return g;
    }

    void set_Shader(std::shared_ptr<Shader>& shader){
        this->shader = shader;
        modelLoc = glGetUniformLocation(shader->getID(), "model"); 
    }

    void setup_buffers() {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, this->sphereVertices.size() * sizeof(float), this->sphereVertices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->indices.size() * sizeof(unsigned int), this->indices.data(),GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void make_sphere(unsigned int sectorCount = 30, unsigned int stackCount = 30){
        float sectorStep = 2 * M_PI / sectorCount;
        float stackStep = M_PI / stackCount;
        float sectorAngle, stackAngle;
        sphereVertices.clear();
        indices.clear();
        
        for (int i = 0; i <= stackCount; ++i){
           stackAngle = M_PI / 2 - i * stackStep;

           float xy = radius * cosf(stackAngle);

            for (int j = 0; j <= sectorCount; ++j){
                sectorAngle = j * sectorStep;

               // Inserting x, y, z
               this->sphereVertices.insert(this->sphereVertices.end(), { xy * cosf(sectorAngle), xy * sinf(sectorAngle), radius * sinf(stackAngle)});
            }
        }

        for (int i = 0; i < stackCount; ++i){
           int k1 = i * (sectorCount + 1);
           int k2 = k1 + sectorCount + 1;

            for (int j = 0; j < sectorCount; ++j, ++k1, ++k2){
                this->indices.push_back(k1);
                this->indices.push_back(k2);
                this->indices.push_back(k1 + 1);

                this->indices.push_back(k1 + 1);
                this->indices.push_back(k2);
                this->indices.push_back(k2 + 1);
            }
        }
    }

    void drawObject(glm::mat4 &model) {
        this->model = model;
        shader->use();
        glBindVertexArray(VAO);
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glDrawElements(GL_TRIANGLES, this->indices.size(), GL_UNSIGNED_INT, 0);
    }

    ~Schwarzchild(){
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
        glDeleteVertexArrays(1, &VAO);
    }
};