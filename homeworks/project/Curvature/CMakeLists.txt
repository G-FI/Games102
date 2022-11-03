cmake_minimum_required(VERSION 3.10)

#设置项目名称
project(Curvature)

#生成可执行文件，输入为source files
add_executable(Curvature main.cpp "glad.c")




#添加include directories
include_directories("E:/OpenGL/include")
#target_include_directories("...")为什么不行？



#添加使用的库
list(APPEND EXTRA_LIBS "E:/OpenGL/lib/glfw3.lib")
target_link_libraries(Curvature PUBLIC ${EXTRA_LIBS})