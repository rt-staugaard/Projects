#pragma once
#include <GLFW/glfw3.h>
#include "physics.hpp"
#include "shader.hpp"


class Gravitational_Object{
    public:
    State state;
    std::shared_ptr<spacetime> space;
    float step_size = 0.1;

    Gravitational_Object(std::shared_ptr<spacetime> &space, State state){
        this->state = state;
        this->space = space; 
    }

    Gravitational_Object(State state){
        this->state = state;
    }

    glm::mat4 get_partial_derivative(State state, int sigma){
        float h = 1e-4f;

        State state_plus = state;
        state_plus.pos[sigma] += h;

        State state_minus = state;
        state_minus.pos[sigma] -= h;
      
        return (space->get_metric(state_plus) - space->get_metric(state_minus))/(2*h);
    }

    void get_christoffel_symbol(State state){
        glm::mat4 dg[4];
        for(int i = 0; i < 4; ++i) dg[i] = get_partial_derivative(state, i);  
        
        glm::mat4 g = space->get_metric(state);
        glm::mat4 inv_g = glm::inverse(g);

        float christoffel_symbol = 0;

    
        for(int mu = 0; mu < 4; ++mu){
            for (int alpha = 0; alpha < 4; ++alpha){
                for (int beta = alpha; beta < 4; ++beta){
                    float christoffel_symbol = 0;
                    for (int lambda = 0; lambda < 4; lambda++){
                        christoffel_symbol += 0.5 * inv_g[mu][lambda] * (dg[alpha][lambda][beta] + dg[beta][alpha][lambda] - dg[lambda][alpha][beta]);
                    }
                    space->insert_christoffel_symbol(christoffel_symbol, mu, alpha, beta);
                }
            }
        }
    
    }

    State get_acceleration(State current){
        get_christoffel_symbol(current);
        State change;
        change.pos = current.vel;
        change.vel = glm::vec4(0.0f);

        for (int mu = 0; mu < 4; ++mu){
            for (int alpha = 0; alpha < 4; ++alpha){
                for (int beta = 0; beta < 4; ++beta){
                    change.vel[mu] -= space->get_christoffel_symbol(mu)[alpha][beta] * current.vel[alpha] * current.vel[beta];
                }
            }
        }

        glm::mat4 metric = space->get_metric(current);
        float ds_squared = metric[0][0] * change.vel[0] * change.vel[0];
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                ds_squared += metric[i + 1][j + 1] * change.vel[i + 1] * change.vel[j + 1];
            }
        }

        if (glm::abs(ds_squared) > 1.e-2){
            float magnitude = glm::sqrt(metric[1][1]*change.vel[1]*change.vel[1] + metric[2][2]*change.vel[2]*change.vel[2] + metric[3][3]*change.vel[3]*change.vel[3]);
            float target_magnitude = glm::sqrt(metric[0][0] * change.vel[0] * change.vel[0]);
            float scale = target_magnitude / magnitude;

            change.vel[1] = scale * change.vel[1];
            change.vel[2] = scale * change.vel[2];
            change.vel[3] = scale * change.vel[3];
        }

        return change;
    }

    void update(float h = 0.01){
        State k1 = get_acceleration(state);
        State k2 = get_acceleration(state + k1 * (h/2));
        State k3 = get_acceleration(state + k2 * (h/2));
        State k4 = get_acceleration(state + k3 * h);

        State new_state = state + (k1 + k2 * 2 + k3 * 2 + k4) * (h/6);
        state = new_state;
    }
};

class Ray : public Gravitational_Object{
    private:
    unsigned int VAO, VBO;
    std::shared_ptr<Shader> shader;
    glm::mat4 model;
    unsigned int timeLoc1, timeLoc2, modelLoc;

    public:
    std::vector<glm::vec4> history;
    bool is_inside_horizon = false;

    Ray(std::shared_ptr<spacetime> space , State state, unsigned int data_size = 200) : Gravitational_Object(space, state) {
        setup_buffers(VAO, VBO, data_size);
        this->state = state;
    }

    void set_Shader(std::shared_ptr<Shader>& shader){
        this->shader = shader;
        timeLoc1 = glGetUniformLocation(shader->getID(), "ray_time"); 
        timeLoc2 = glGetUniformLocation(shader->getID(), "real_time");
        modelLoc = glGetUniformLocation(shader->getID(), "model"); 
    }

    void setup_buffers(unsigned int& VAO, unsigned int& VBO, unsigned int& data_size) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, data_size * sizeof(glm::vec4), NULL, GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void propagate_ray(glm::mat4 &model){
        update();

        float radius = 2;
        if (state.pos.y * state.pos.y + state.pos.z * state.pos.z + state.pos.w * state.pos.w <= radius * radius) {
            is_inside_horizon = true;
        }

        if (!is_inside_horizon) {
            history.push_back(state.pos);
        }

        if (history.size() > 200) {
            history.erase(history.begin());
        }

        if (history.size() >= 2) {
            drawTrajectory(model);
        }
    }

    void drawTrajectory(glm::mat4 &model) {
        this->model = model;
        shader->use();
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glUniform1f(timeLoc1, this->state.pos[0]);
        glUniform1f(timeLoc2, (float)glfwGetTime());
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));


        glBufferSubData(GL_ARRAY_BUFFER, 0, history.size() * sizeof(glm::vec4), history.data());
        glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)history.size());
    }
    
    ~Ray() {
        glDeleteBuffers(1, &VBO);
        glDeleteVertexArrays(1, &VAO);
    }
};

class Star : public Gravitational_Object {
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

    Star(State state, float mass, float radius, unsigned int sectorCount = 30, unsigned int stackCount = 30) : Gravitational_Object(state) {
        this->state = state;
        this->radius = radius;
        this->mass = mass;
        make_sphere(sectorCount, stackCount);
        setup_buffers();
    }

    void set_Shader(std::shared_ptr<Shader>& shader){
        this->shader = shader;
        modelLoc = glGetUniformLocation(shader->getID(), "model"); 
    }

    void set_Spacetime(std::shared_ptr<spacetime> space){
        this->space = space;
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

    void make_sphere(unsigned int sectorCount, unsigned int stackCount){
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

    ~Star(){
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
        glDeleteVertexArrays(1, &VAO);
    }
};
