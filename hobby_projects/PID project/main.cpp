#include "modules/shader.hpp"
#include <array>
#include <cstdlib>	
#include <ctime>
#include <GLFW/glfw3.h>
#include <glm.hpp>

const float roomTemperature = 10.0f;

template <int Width, int Height>
struct Room{
    std::array<float, Width * Height> TempGrid;

    Room(){
        for (int i = 0; i < Width * Height; ++i){
            float noise = float(rand() % 21 - 10 ) / 5;  // [- 2 , 2 ]
            TempGrid[i] = roomTemperature + noise;
        }
    }

    float getTemperature(){
        float result = 0;
        for (int i = 0; i < Width * Height; ++i){
            result += TempGrid[i];
        }
        return result / (Width * Height); 
    }
};

template <int Width, int Height>
struct Heater{
    int index;
    float temp;
    float maxTemperature = 40.0f;

    Heater(int x, int y, float temp) : temp(temp){
        index = x + y * Width; 
    };

    void changeTemp(const float power){
        if (power == 0.0f){
            this->temp = roomTemperature;
        }
        else{
            this->temp = roomTemperature + power * maxTemperature;
        }
    }
};

template <int Width, int Height>
void update (Room<Width, Height> &room, Heater<Width, Height> heater, float alpha, float dt){
    std::array<float, Width * Height> nextGrid = room.TempGrid;

    for (int x = 0; x < Width; ++x){
        for (int y = 0; y <  Height; ++y){
            int idx = y * Width + x;

            float leakage = room.TempGrid[idx] * 0.995 + roomTemperature * 0.005;

            float T_up = (y > 0) ? room.TempGrid[idx - Width] : leakage;
            float T_down = (y < Height - 1) ?  room.TempGrid[idx + Width] : leakage;
            float T_right = (x < Width - 1) ? room.TempGrid[idx + 1] : leakage;
            float T_left = (x > 0) ? room.TempGrid[idx - 1] : leakage;

            float laplacian = (T_right + T_left - 2 * room.TempGrid[idx]) + 
                              (T_up + T_down - 2 * room.TempGrid[idx]);


            nextGrid[idx] += alpha * laplacian * dt;
        }
    }
    nextGrid[heater.index] = heater.temp;
    room.TempGrid = nextGrid;    
}

template <int Width, int Height>
float controlFunction(const float error, const float pastError, float futureChange, const std::array<float, 3> constants, float alpha, float dt){

    float PID_VALUE = constants[0] * error + constants[1] * pastError + constants[2] * futureChange;

    if (PID_VALUE > 100.0f) return 100.0f; 
    if (PID_VALUE < 0.0f)   return 0.0f;

    return PID_VALUE;
}

void setup_buffers(unsigned int &VAO, unsigned int &VBO) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 18 * sizeof(float), NULL, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
}

void draw(std::unique_ptr<Shader> &shader, const float vertices[18], unsigned int &VAO, unsigned int &VBO) {
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    glBufferSubData(GL_ARRAY_BUFFER, 0,  18 * sizeof(float), &vertices[0]);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}

constexpr int room_width = 20;
constexpr int room_height = 20;

constexpr float screenWidth = 600;
constexpr float screenHeight = 600;

int main(int argc, char* argv[]){
    srand(time(0));

    GLFWwindow* window;

    if (!glfwInit()){
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(screenWidth, screenHeight, "Simulation Window", NULL, NULL);
    if (!window){
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    std::filesystem::path exe_path = std::filesystem::absolute(argv[0]).parent_path();
    std::filesystem::path roomVertDir = exe_path / "shaders" / "room.vert";
    std::filesystem::path roomFragDir = exe_path / "shaders" / "room.frag";
    auto roomShader = std::make_unique<Shader>(roomVertDir.string(), roomFragDir.string());
    auto shaderID = roomShader->getID();
    roomShader->use();

    const float vertices[] = {
        -1.0f,  1.0f, 0.0f,  
        -1.0f, -1.0f, 0.0f,  
         1.0f, -1.0f, 0.0f,  

        -1.0f,  1.0f, 0.0f,  
        1.0f, -1.0f, 0.0f,  
        1.0f,  1.0f, 0.0f  
    };

    unsigned int VAO, VBO;
    setup_buffers(VAO, VBO);

    Room<room_width,room_height> room;
    Heater<room_width,room_height> heater(12, 7, 25);

    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, room_width, room_height, 0, GL_RED, GL_FLOAT, room.TempGrid.data());

    std::array<float, 5> errorHistory;
    float SP = 22.0;
    float alpha = 0.35;
    std::array<float, 3> constants = {15, 0.75, 5};
    
    int head = 0;
    float PV = room.getTemperature();
    float error = SP - PV;
    float lastError = 0;

    double lastTime = glfwGetTime();
    int nbFrames = 0;
    while (!glfwWindowShouldClose(window)){
        double currentTime = glfwGetTime();
        nbFrames++;

        float deltaTime = currentTime - lastTime;

        errorHistory[head] = error * deltaTime;
        head = (head + 1) % 5;
        float pastError = 0;
        for(int i = 0; i < errorHistory.size(); ++i){pastError += errorHistory[i];}

        update<room_width,room_height>(room, heater, alpha, deltaTime);

        PV = room.getTemperature();
        error = SP - PV;

        float futureChange = (error - lastError) / (deltaTime + 0.01);

        float PID_VALUE = controlFunction<room_width,room_height>(error, pastError, futureChange, constants, alpha, deltaTime);
        lastError = error;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, room_width, room_height, GL_RED, GL_FLOAT, room.TempGrid.data());

        if (deltaTime >= 1.0) { 
            double msPerFrame = 1000.0 / double(nbFrames);
            std::cout << msPerFrame << " ms/frame (" << nbFrames << " FPS)" << std::endl;
            std::cout << "Room Temperatur: " << PV << ", " << "PID_VALUE" << PID_VALUE <<"\n";
            nbFrames = 0;
            lastTime = currentTime; 
        }    
        
        glClear(GL_COLOR_BUFFER_BIT);
        draw(roomShader, vertices, VAO, VBO);

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (currentTime > 90.0) break;
    }

    glfwTerminate();
    return 0;
};