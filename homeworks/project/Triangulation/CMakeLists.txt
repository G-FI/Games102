cmake_minimum_required(VERSION 3.10)

project(Triangulation)

set(CMAKE_CXX_STANDARD 17)
add_executable(Triangulation "main.cpp")

add_subdirectory("${PROJECT_SOURCE_DIR}/gui")

target_include_directories(Triangulation PUBLIC "${PROJECT_SOURCE_DIR}/gui")

include_directories("D:/CGAL-5.5.1/include" "D:/boost_1_80_0" "D:/include")

#opengl library
list(APPEND EXTRA_LIBS "D:/lib/glfw/glfw3.lib" "GUI")

target_link_libraries(${PROJECT_NAME} PUBLIC ${EXTRA_LIBS})