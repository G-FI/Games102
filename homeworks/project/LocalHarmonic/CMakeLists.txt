cmake_minimum_required(VERSION 3.10)

#设置项目名称
project(Curvature)

#生成可执行文件，输入为source files
add_executable(Curvature "main.cpp")

#添加include directories
include_directories("D:/include" "D:/boost_1_80_0" "D:/CGAL-5.5.1/include")

#start引入libigl viewer模块
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" ${LIBIGL_TOPLEVEL_PROJECT})

option(LIBIGL_OPENGL  "use opengl" ON)
option(LIBIGL_GLFW  "use GLFW" ON)
option(LIBIGL_IMGUI "use Imgui" ON)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(libigl)
igl_include(glfw)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl::imgui)
#end


#添加使用的库
list(APPEND EXTRA_LIBS "D:/lib/glfw/glfw3.lib")
target_link_libraries(Curvature PUBLIC ${EXTRA_LIBS})
