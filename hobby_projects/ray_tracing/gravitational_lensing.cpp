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

unsigned int make_shader(const std::string& vertex_filepath, const std::string fragment_filepath);
unsigned int make_module(const std::string& filepath, unsigned int module_type);
void send_geometry(const std::vector<glm::vec4>& vertices, unsigned int& VAO, unsigned int& VBO);
void send_geometry_with_size(const std::vector<glm::vec4>& vertices, unsigned int& VAO, unsigned int& VBO);

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
    float step_size = 0.1;

    Gravitational_Object(glm::vec4 &pos, glm::vec4 &vel){
        this->state.pos = pos;
        this->state.vel = vel;  
    }

    glm::mat4 get_partial_derivative(spacetime &space, State state, int sigma){
        float h = 1e-4f;

        State state_plus = state;
        state_plus.pos[sigma] += h;

        State state_minus = state;
        state_minus.pos[sigma] -= h;
      
        return (space.get_metric(state_plus) - space.get_metric(state_minus))/(2*h);
    }

    void get_christoffel_symbol(spacetime &space, State state){
        glm::mat4 dg[4];
        for(int i = 0; i < 4; ++i) dg[i] = get_partial_derivative(space, state, i);  
        
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

    State get_acceleration(spacetime &space, State current){
        get_christoffel_symbol(space, current);
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
        change.pos[3] = 0;
        change.vel[3] = 0;

        return change;
    }

    void update(spacetime &space, State state, float h = 0.01){
        State k1 = get_acceleration(space, state);
        State k2 = get_acceleration(space, state + k1 * (h/2));
        State k3 = get_acceleration(space, state + k2 * (h/2));
        State k4 = get_acceleration(space, state + k3 * h);

        State new_state = state + (k1 + k2 * 2 + k3 * 2 + k4) * (h/6);
        this->state = new_state;
    }
};

class Stellar_Object : public Gravitational_Object {
public:
    float mass;
    unsigned int shader;
    int colorLoc, projLoc;

    Stellar_Object(float mass, glm::vec4 pos, glm::vec4 vel) 
        : Gravitational_Object(pos, vel), mass(mass) 
    {
        shader = make_shader("shaders/circle.vert", "shaders/circle.frag");
        colorLoc = glGetUniformLocation(shader, "objectColor");
        projLoc = glGetUniformLocation(shader, "projection");
    }

    void drawObject(GLuint vao, GLuint vbo) {
        glUseProgram(shader);
        glUniform3f(colorLoc, 1.0f, 0.0f, 0.0f);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo); 
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec4), &state.pos);
        glDrawArrays(GL_POINTS, 0, 1);
    }
    
    void deleteProgram() { glDeleteProgram(shader); }
};

class Ray : public Gravitational_Object{
    public:
    unsigned int shader;
    int colorLoc, timeLoc;

    Ray(glm::vec4 pos, glm::vec4 vel, bool initialize_shader = true) : Gravitational_Object(pos, vel) {
        if (initialize_shader == true){
            shader = make_shader("shaders/line.vert", "shaders/line.frag");
            colorLoc = glGetUniformLocation(shader, "objectColor"); 
            timeLoc = glGetUniformLocation(shader, "currentTime");
        }
    }

    void drawTrajectory(const std::vector<glm::vec4>& points, GLuint vao, GLuint vbo) {
        glUseProgram(shader);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        //glm::mat4 projection = glm::ortho(-0.01f, 0.01f, -0.01f, 0.01f, -1.0f, 1.0f);
        glUniform3f(colorLoc, 0.0f, 1.0f, 1.0f);
        glUniform1f(timeLoc, (float)glfwGetTime());

        glBufferSubData(GL_ARRAY_BUFFER, 0, points.size() * sizeof(glm::vec4), points.data());
        glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)points.size());
    }
    
    void deleteProgram(){
        glDeleteProgram(shader);
    }
};

struct RayStruct{
    Ray *ray;
    unsigned int VAO, VBO;
    std::vector<glm::vec4> history;
    State current_state;
};


int main(void)
{   
    GLFWwindow* window;

    if (!glfwInit()){
        return -1;
    }
        
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


    window = glfwCreateWindow(640, 480, "Simulation Window", NULL, NULL);
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
    
    // Making central object 
    glEnable(GL_PROGRAM_POINT_SIZE);

    Stellar_Object obj1 = Stellar_Object(1, glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f}, glm::vec4 {0.0f, 0.0f, 0.0f, 0.0f});
        
    unsigned int VAO_circle;
    unsigned int VBO_circle;

    std::vector<glm::vec4> circle_data;
    circle_data.push_back({0.0f, 0.0f, 0.0f, 0.0f});

    send_geometry(circle_data, VAO_circle, VBO_circle);

    // Making ray bundle
    int ray_number = 20;
    std::vector<RayStruct> rayBundle;
    rayBundle.reserve(ray_number);

    for (int i = 0; i < ray_number; ++i){
        State initial_state(glm::vec4{0.0f, -20.0f, -10.0f + (i * 1), 0.0f}, 
                      glm::vec4{1.0f, 1, 0.0f, 0.0f});

        RayStruct ri;
        ri.ray = new Ray(initial_state.pos, initial_state.vel);
        ri.current_state = initial_state;
        
        ri.history.push_back(initial_state.down_scale(20));

        send_geometry_with_size(ri.history, ri.VAO, ri.VBO);
    
        rayBundle.push_back(ri);
    }

    spacetime space;
    while (!glfwWindowShouldClose(window))
    {
        // Exportation of graphics to GPU
        
        glClear(GL_COLOR_BUFFER_BIT);

        obj1.drawObject(VAO_circle, VBO_circle);

        for (auto& ri : rayBundle) {
            ri.ray->update(space, ri.current_state);
            ri.current_state = ri.ray->state;

            ri.history.push_back(ri.current_state.down_scale(20));
            if (ri.history.size() > 200){
                ri.history.erase(ri.history.begin());
            } 
            
            if (ri.history.size() >= 2) {
                ri.ray->drawTrajectory(ri.history, ri.VAO, ri.VBO);
            }
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    obj1.deleteProgram();
    for (auto ri : rayBundle){
        ri.ray->deleteProgram();
        delete ri.ray;
    }

    glfwTerminate();
    return 0;
}

unsigned int make_module(const std::string& filepath, unsigned int module_type){

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

unsigned int make_shader(const std::string& vertex_filepath, const std::string fragment_filepath){

    std::vector<unsigned int> modules;
    modules.push_back(make_module(vertex_filepath, GL_VERTEX_SHADER));
    modules.push_back(make_module(fragment_filepath, GL_FRAGMENT_SHADER));

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

    return shaderProgram;
}

void send_geometry(const std::vector<glm::vec4>& vertices, unsigned int& VAO, unsigned int& VBO) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec4), vertices.data(), GL_DYNAMIC_DRAW);

    // Position Attribute
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glEnableVertexAttribArray(0);
}

void send_geometry_with_size(const std::vector<glm::vec4>& vertices, unsigned int& VAO, unsigned int& VBO) {
    int size = 200;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, size * sizeof(glm::vec4), NULL, GL_DYNAMIC_DRAW);

    // Position Attribute
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glEnableVertexAttribArray(0);
}