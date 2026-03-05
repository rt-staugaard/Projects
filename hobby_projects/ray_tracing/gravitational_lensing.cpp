#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>


class Shader{
    private:
    uint shaderProgram;

    public:

    Shader() : shaderProgram(0) {}

    Shader(const std::string& vertex_filepath, const std::string fragment_filepath){
        make_shader_program(vertex_filepath, fragment_filepath);
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

    void make_shader_program(const std::string& vertex_filepath, const std::string fragment_filepath){
        std::vector<unsigned int> modules;
        modules.push_back(make_shader_module(vertex_filepath, GL_VERTEX_SHADER));
        modules.push_back(make_shader_module(fragment_filepath, GL_FRAGMENT_SHADER));
        
        unsigned int shaderProgram = glCreateProgram();
        for (unsigned int shaderModule : modules){
            glAttachShader(shaderProgram, shaderModule);
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

    uint getID() const { return shaderProgram; }

    ~Shader() {
        glDeleteProgram(shaderProgram);
    }
};

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

    glm::vec4 down_scale(float scale = 100){
        glm::vec4 scalled_pos = this->pos;
        scalled_pos.y = this->pos.y/scale;
        scalled_pos.z = this->pos.z/scale;
        scalled_pos.w = this->pos.w/scale;
        return scalled_pos;
    }

};

// For now specifically Schwarzchild spacetime
struct spacetime{
    float M = 1.0;
    float G = 1.0;
    glm::mat4 gamma[4];

    glm::mat4 get_metric(State state){
        glm::vec4 pos = state.pos;
        float x = state.pos.y;
        float y = state.pos.z;
        float z = state.pos.w;

        float R = std::sqrt(x * x + y * y + z * z);

        glm::mat4 g(0.0);

        g[0][0] = -(1.0 - 2.0 * M / R);

        float factor = (2.0 * M) / (R * R * (R - 2.0 * M));

        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
                float delta_ij = (i == j) ? 1.0 : 0.0;
                g[i][j] = delta_ij + factor * pos[i] * pos[j];
            }
        }
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

};

class Gravitational_Object{
    public:
    State state;
    spacetime space;
    float step_size = 0.1;

    Gravitational_Object(State state){
        this->state = state;
    }

    glm::mat4 get_partial_derivative(State state, int sigma){
        float h = 1e-4f;

        State state_plus = state;
        state_plus.pos[sigma] += h;

        State state_minus = state;
        state_minus.pos[sigma] -= h;
      
        return (space.get_metric(state_plus) - space.get_metric(state_minus))/(2*h);
    }

    void get_christoffel_symbol(State state){
        glm::mat4 dg[4];
        for(int i = 0; i < 4; ++i) dg[i] = get_partial_derivative(state, i);  
        
        glm::mat4 g = space.get_metric(state);
        glm::mat4 inv_g = glm::inverse(g);

        float christoffel_symbol = 0;

    
        for(int mu = 0; mu < 4; ++mu){
            for (int alpha = 0; alpha < 4; ++alpha){
                for (int beta = alpha; beta < 4; ++beta){
                    float christoffel_symbol = 0;
                    for (int lambda = 0; lambda < 4; lambda++){
                        christoffel_symbol += 0.5 * inv_g[mu][lambda] * (dg[alpha][lambda][beta] + dg[beta][alpha][lambda] - dg[lambda][alpha][beta]);
                    }
                    space.insert_christoffel_symbol(christoffel_symbol, mu, alpha, beta);
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
                    change.vel[mu] -= space.get_christoffel_symbol(mu)[alpha][beta] * current.vel[alpha] * current.vel[beta];
                }
            }
        }

        glm::mat4 metric = space.get_metric(current);
        float ds_squared = metric[0][0] * change.vel[0] * change.vel[0];
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                ds_squared += metric[i + 1][j + 1] * change.vel[i + 1] * change.vel[j + 1];
            }
        }

        if (glm::abs(ds_squared) > 1.e5){
            float magnitude = glm::sqrt(metric[1][1]*change.vel[1]*change.vel[1] + metric[2][2]*change.vel[2]*change.vel[2] + metric[3][3]*change.vel[3]*change.vel[3]);
            float target_magnitude = glm::sqrt(metric[0][0] * change.vel[0] * change.vel[0]);
            float scale = target_magnitude / magnitude;

            change.vel[1] = scale * change.vel[1];
            change.vel[2] = scale * change.vel[2];
            change.vel[3] = scale * change.vel[3];
        }

        return change;
    }

    void update(State state, float h = 0.01){
        State k1 = get_acceleration(state);
        State k2 = get_acceleration(state + k1 * (h/2));
        State k3 = get_acceleration(state + k2 * (h/2));
        State k4 = get_acceleration(state + k3 * h);

        State new_state = state + (k1 + k2 * 2 + k3 * 2 + k4) * (h/6);
        this->state = new_state;
    }
};

class Stellar_Object : public Gravitational_Object {
public:
    private:
    uint VAO, VBO, EBO;
    std::shared_ptr<Shader> shader;
    glm::mat4 model;
    std::vector<float> sphereVertices;
    std::vector<uint> indices;
    uint modelLoc;

    public:
    float mass;
    float radius;

    Stellar_Object(float mass, State state, float radius = 2, uint sectorCount = 30, uint stackCount = 30) 
        : Gravitational_Object(state), 
          mass(mass), 
          radius(glm::max(radius, 2 * mass))
    {
        make_sphere(sectorCount, stackCount);
        setup_buffers();
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
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->indices.size() * sizeof(uint), this->indices.data(),GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void make_sphere(uint sectorCount = 30, uint stackCount = 30){
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

    ~Stellar_Object() {
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
        glDeleteVertexArrays(1, &VAO);
    }
};

class Ray : public Gravitational_Object{
    private:
    uint VAO, VBO;
    std::shared_ptr<Shader> shader;
    glm::mat4 model;
    uint timeLoc1, timeLoc2, modelLoc;

    public:
    Ray(State state, uint data_size = 200) : Gravitational_Object(state) {
        setup_buffers(VAO, VBO, data_size);
    }

    void set_Shader(std::shared_ptr<Shader>& shader){
        this->shader = shader;
        timeLoc1 = glGetUniformLocation(shader->getID(), "ray_time"); 
        timeLoc2 = glGetUniformLocation(shader->getID(), "real_time");
        modelLoc = glGetUniformLocation(shader->getID(), "model"); 
    }

    void setup_buffers(uint& VAO, uint& VBO, uint& data_size) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, data_size * sizeof(glm::vec4), NULL, GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void drawTrajectory(const std::vector<glm::vec4>& points, glm::mat4 &model) {
        this->model = model;
        shader->use();
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glUniform1f(timeLoc1, this->state.pos[0]);
        glUniform1f(timeLoc2, (float)glfwGetTime());
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));


        glBufferSubData(GL_ARRAY_BUFFER, 0, points.size() * sizeof(glm::vec4), points.data());
        glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)points.size());
    }
    
    ~Ray() {
        glDeleteBuffers(1, &VBO);
        glDeleteVertexArrays(1, &VAO);
    }
};

struct RayStruct{
    Ray ray;
    std::vector<glm::vec4> history;
    State current_state;
    bool is_inside_horizon = false;

    RayStruct(State &state, std::shared_ptr<Shader> &shader) 
        : ray(state), 
          current_state(state) 
    {
        ray.set_Shader(shader);
        history.push_back(state.pos);
    }
};

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


int main(void)
{   
    GLFWwindow* window;

    if (!glfwInit()){
        return -1;
    }
        
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


    window = glfwCreateWindow(600, 600, "Simulation Window", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);


    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    std::cout << glGetString(GL_VERSION) << std::endl;

    // Setting up camera
    glm::vec3 cameraPostion{0,0,0};
    glm::vec3 cameraFront{0,0,-4};
    glm::vec3 cameraUp{0,1,0};
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 600.0f / 600.0f, 0.1f, 100.0f);
    Camera camera(cameraPostion, cameraFront, cameraUp);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);  
    glfwSetWindowUserPointer(window, &camera);
    glfwSetCursorPosCallback(window, mouse_callback_bridge);
    glfwSetScrollCallback(window, scroll_callback_bridge); 

    // Making Spherical Central Object 
    glEnable(GL_PROGRAM_POINT_SIZE);

    uint sectorCount = 30;
    uint stackCount = 30;

    State object_state = State(glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f}, glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f});
    Stellar_Object obj1 = Stellar_Object(1, object_state);
    float radius = obj1.radius;
    auto circle_shader = std::make_shared<Shader>("shaders/circle.vert", "shaders/circle.frag");
    obj1.set_Shader(circle_shader);

    // Making Ray Bundle
    int number_of_rays = 10;
    std::vector<std::unique_ptr<RayStruct>> rayBundle;
    rayBundle.reserve(number_of_rays);
    auto ray_shader = std::make_shared<Shader>("shaders/line.vert", "shaders/line.frag");

    for (int i = 0; i < number_of_rays; ++i){
        for(int j = 0; j < number_of_rays; ++j){
            State initial_state(glm::vec4{0.0f, -10.0f, i-4, j-4}, glm::vec4{1.0f, 1, 0.0f, 0.0f});
            rayBundle.push_back(std::make_unique<RayStruct>(initial_state, ray_shader));
        }
    }

    glEnable(GL_DEPTH_TEST);
    float lastFrame = 0;
    glm::mat4 model;
    glm::mat4 object_matrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -30.0f));
    while (!glfwWindowShouldClose(window))
    {
        // Exportation of graphics to GPU
        float currentFrame = glfwGetTime();
        float deltaTime = currentFrame - lastFrame;
        processInput(window, camera, deltaTime);
        lastFrame = currentFrame;  

        glm::mat4 model = camera.transform_Model(object_matrix);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        obj1.drawObject(model);
        for (auto& ri : rayBundle) {
            ri->ray.update(ri->current_state);
            ri->current_state = ri->ray.state;

            if (!(ri->is_inside_horizon)){
                ri->history.push_back(ri->current_state.pos);
            }
            else if (ri->current_state.pos[1] * ri->current_state.pos[1] + ri->current_state.pos[2] * ri->current_state.pos[2] + ri->current_state.pos[3] * ri->current_state.pos[3] <= radius * radius){
                ri->is_inside_horizon = true;
            }
            
            if (ri->history.size() > 200){
                ri->history.erase(ri->history.begin());
            } 
            
            if (ri->history.size() >= 2) {
                ri->ray.drawTrajectory(ri->history, model);
                continue;
            }
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}