cmake_minimum_required(VERSION 3.23)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Chip-8 VERSION 1.0)

find_package(glfw3 REQUIRED)

# GLAD LIBRARY
add_library(GLAD ../../dependencies/glad/src/glad.c)
target_include_directories(GLAD PUBLIC ../../dependencies/glad/include)

# GRAPHICS MODULE
add_library(GraphicsModule INTERFACE)
target_sources(GraphicsModule INTERFACE modules/Graphics.hpp)
target_link_libraries(GraphicsModule INTERFACE GLAD)
target_include_directories(GraphicsModule INTERFACE ../../dependencies/shader)
target_include_directories(GraphicsModule INTERFACE ../../dependencies/glm)

# EMULATOR MODULE
add_library(EmulatorModule INTERFACE)
target_sources(EmulatorModule INTERFACE modules/Chip-8.hpp)

# GLFW MODULE
add_library(GLFWModule INTERFACE)
target_sources(GLFWModule INTERFACE modules/GLFW_setup.hpp)
target_link_libraries(GLFWModule INTERFACE GLAD)
target_link_libraries(GLFWModule INTERFACE glfw)

add_executable(Chip-8 main.cpp)
target_link_libraries(Chip-8 PRIVATE 
    GraphicsModule 
    GLFWModule 
    EmulatorModule
)
