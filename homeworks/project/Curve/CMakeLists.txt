cmake_minimum_required(VERSION 3.10)

project(Curve)

set(CMAKE_CXX_STANDARD 17)
add_executable(Curve "main.cpp")

add_subdirectory("${PROJECT_SOURCE_DIR}/gui")

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_SOURCE_DIR}/gui")

include_directories("D:/include")

#opengl library
list(APPEND EXTRA_LIBS "D:/lib/glfw/glfw3.lib" "GUI")

target_link_libraries(${PROJECT_NAME} PUBLIC ${EXTRA_LIBS})