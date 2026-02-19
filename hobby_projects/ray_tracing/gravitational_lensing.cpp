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
void send_data(std::vector<float> &data, unsigned int &VAO, unsigned int &VBO);


// For now specifically Schwarzchild spacetime
struct spacetime{
    float M = 1.0;
    float G = 1.0;

    glm::mat4 get_metric(glm::vec4 &pos){
        float r = pos.y;
        float sinTheta = glm::sin(pos.z);

        glm::mat4 metric(0.0f);

        float r_s = (2 * G * M)/r;

        metric[0][0] = -(1 - r_s);
        metric[1][1] = 1/(1 - r_s);
        metric[0][0] = r * r;
        metric[0][0] = r * r * sinTheta * sinTheta;

        return metric; 
    }

};

class Gravitational_Object{
    public:
    glm::vec3 pos;
    glm::vec3 vel;
    float step_size = 0.1;

    Gravitational_Object(glm::vec3 pos, glm::vec3 vel){
        this->pos = pos;
        this->vel = vel;
    }

    void update(glm::vec3 acc, float dt){
        vel += acc * dt;
        pos += vel * dt;
    }

    glm::mat4 get_partial_derivative(glm::vec4 pos, int sigma){
        spacetime spacetime;

        



        glm::mat4 metric = spacetime.get_metric(pos);


        glm::mat4 partial_derivative;


    }

};

class Stellar_Object : public Gravitational_Object{
    public:
    float mass;

    Stellar_Object(float mass, glm::vec3 pos, glm::vec3 vel): Gravitational_Object(pos, vel), mass(mass) 
    {this->mass = mass;}
};

class Ray : public Gravitational_Object{
    Ray(glm::vec3 pos, glm::vec3 vel) : Gravitational_Object(pos, vel) {}
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
    glEnable(GL_PROGRAM_POINT_SIZE);

    Stellar_Object obj1 = Stellar_Object(1, glm::vec3 {-0.8, 0.0f, 0.0f}, glm::vec3 {0.0f, 0.5, 0.0f});
    Stellar_Object obj2 = Stellar_Object(5, glm::vec3 {0.0f, 0.0f, 0.0f}, glm::vec3 {0.0f, 0.0f, 0.0f});

    std::vector<float> data = {0.0f, 0.0f, 0.0f};  
        
    unsigned int VAO;
    unsigned int VBO;
    
    send_data(data, VAO, VBO);
    unsigned int shader = make_shader("shader.vert", "shader.frag");
    int modelLoc = glGetUniformLocation(shader, "model");
    int colorLoc = glGetUniformLocation(shader, "objectColor");
    int vertexCount = 1;

    while (!glfwWindowShouldClose(window))
    {
        // Simulation logic computed on CPU
       glm::vec3 relative_vec = obj2.pos - obj1.pos;
       float r = glm::length(relative_vec);
       
       float G = 0.1f; 
       float epsilon = 0.1f;
       float force_mag = G * (obj1.mass * obj2.mass) / (r * r + epsilon);
       
       glm::vec3 force_dir = glm::normalize(relative_vec);
       
       glm::vec3 acc1 = (force_mag / obj1.mass) * force_dir;
       glm::vec3 acc2 = (force_mag / obj2.mass) * (-force_dir); 
       
       obj1.update(acc1,0.01);
       obj2.update(acc2,0.01);

       glm::vec3 positionA = obj1.pos;
       glm::vec3 positionB = obj2.pos;

       glm::vec3 centerOfMass = (obj1.pos * obj1.mass + obj2.pos * obj2.mass) / (obj1.mass + obj2.mass);

        // Exportation of graphics to GPU
        glClear(GL_COLOR_BUFFER_BIT);

        glEnable(GL_PROGRAM_POINT_SIZE);
        glUseProgram(shader);
        glBindVertexArray(VAO); 
        
        // Large object in the middle
        glm::mat4 modelA = glm::mat4(1.0f);
        modelA = glm::translate(modelA, positionA - centerOfMass); 
        modelA = glm::scale(modelA, glm::vec3(1.5f));
        glUniform3f(colorLoc, 1.0f, 0.8f, 0.0f); 
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(modelA));
        glDrawArrays(GL_POINTS, 0, vertexCount);
        
        // Smaller object on the side
        glm::mat4 modelB = glm::mat4(1.0f);
        modelB = glm::translate(modelB, positionB - centerOfMass); 
        modelB = glm::scale(modelB, glm::vec3(0.5f));
        glUniform3f(colorLoc, 0.0f, 0.5f, 1.0f); 
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(modelB));
        glDrawArrays(GL_POINTS, 0, vertexCount);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    glDeleteProgram(shader);
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

void send_data(std::vector<float> &data, unsigned int &VAO, unsigned int &VBO) {
    int size_of_data = 6;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(float),data.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, size_of_data * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    /// glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, size_of_data * sizeof(float), (void*)(3 * sizeof(float)));
    /// glEnableVertexAttribArray(1);
}