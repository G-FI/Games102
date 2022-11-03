#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include<glad/glad.h>
#include <GLFW/glfw3.h> 


#include"CanvasData.h"
#include<Eigen/Core>
#include<iostream>
#include<random>
using namespace Eigen;
using namespace std;

GLFWwindow* createWindow();

ImVec2 operator -(const ImVec2& lhs, const ImVec2& rhs) {
    return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y);
}

ImVec2 operator +(const ImVec2& lhs, const ImVec2& rhs) {
    return ImVec2(lhs.x + rhs.x, lhs.y + rhs.y);
}

ImVec2 operator *(float lhs, const ImVec2& rhs) {
    return ImVec2(lhs * rhs.x, lhs * rhs.y);
}

void NextFrame(GLFWwindow* window) {
    //轮询事件
    glfwPollEvents();
    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}
void EndFrame(GLFWwindow* window) {
    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);
}
void CleanContext(GLFWwindow* window) {
    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}

void Monitor(bool* p_open, CanvasData* data) {
    ImGuiIO& io = ImGui::GetIO();
    //ImGuiWindowFlags window_flags = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav;
    //const float PAD = 10.f;
    //ImVec2 work_pos = ImGui::GetMainViewport()->WorkPos;
    //work_pos.x += PAD;
    //work_pos.y += PAD;

    //ImVec2 work_pos_pivot(0.f, 0.f);
    //ImGui::SetNextWindowPos(work_pos, ImGuiCond_Always, work_pos_pivot);
    //window_flags |= ImGuiWindowFlags_NoMove;


    if (ImGui::Begin("moniter", p_open)) {
        ImGui::Text("Mouse pos(%f, %f)", io.MousePos.x, io.MousePos.y);
    }
    ImGui::End();
}

int main() {
    GLFWwindow* window = createWindow();

    {

        const char* glsl_version = "#version 130";
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;

        io.AddKeyEvent(ImGuiKey_Escape, true);


        // Setup Dear ImGui style
        ImGui::StyleColorsDark();

        // Setup Platform/Renderer backends
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init(glsl_version);
    }

    CanvasData data;
    while (!glfwWindowShouldClose(window)) {

        NextFrame(window);
        bool is_hoverd = false;
        bool is_clicked = false;
        bool is_actived = false;

        //逻辑代码
        if (ImGui::Begin("Args")) {
            ImGui::SliderInt("epochs", &data.epochs, 1, 50);
            ImGui::SliderInt("rbf number", &data.rbf_number, 1, 50);
            ImGui::SliderFloat("lr", &data.lr, 0.0001, 1, "%.4f");
            ImGui::Checkbox("draw line", &data.draw_line);
        }
        ImGui::End();
        
        //创建神经网络
        RBFNetwork model(data.rbf_number, data.epochs, data.lr);

        if (ImGui::Begin("Canvas")) {
            if (data.monitor_open) {
                Monitor(&data.monitor_open, &data);
            }

            // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
            ImVec2 origin = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
            ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
            if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
            if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;

            //绘制Canvas
            ImGuiIO& io = ImGui::GetIO();
            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(origin, origin + canvas_sz, IM_COL32(50, 50, 50, 255));
            draw_list->AddRect(origin, origin + canvas_sz, IM_COL32(255, 255, 255, 255));

            //获取是否点击（使用canvas_sz是，input之下的那些可用区域）
            ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft);
            ImVec2 mouse_in_canvas = io.MousePos - origin;

            if (!data.is_editing && ImGui::IsItemHovered() && ImGui::IsItemClicked()) {
                data.points.push_back(mouse_in_canvas);
                data.trained = false;
            }
            //交互修改控制点
            {
                if (data.is_editing && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    for (int i = 0; i < data.points.size(); ++i) {
                        ImVec2 min = data.points[i] + origin + ImVec2(-5.f, -5.f);
                        ImVec2 max = min + ImVec2(10.f, 10.f);
                        if (ImGui::IsMouseHoveringRect(min, max)) {
                            data.is_draging = true;
                            data.edit_idx = i;
                            break;
                        }
                    }
                }

                if (data.is_draging) {
                    data.points[data.edit_idx] = io.MousePos - origin;
                    data.trained = false;
                }
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                    data.is_draging = false;
                    data.trained = false;
                }
            }

            if (!data.trained && data.points.size() > 3) {
                model.train(data.points);
                data.trained = true;
            }

           

           
            //context popup
            {
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Right))
                    ImGui::OpenPopupContextItem("context");
                if (ImGui::BeginPopup("context")) {
                    ImGui::Selectable("edit", &data.is_editing);
                    if (ImGui::Button("remove one")) {
                        if (data.points.size() > 1) {
                            data.points.pop_back();
                        }
                    }
                    if (ImGui::Button("remove all")) {
                        data.Reset();
                    }

                    ImGui::EndPopup();
                }
            }

            //绘制
            for (const auto& point : data.points) {
                draw_list->AddCircleFilled(point + origin, 4.f, IM_COL32(255, 0, 0, 255));
            }
            if (data.points.size() > 1) {
                if (data.draw_line) {
                    for (int i = 0; i < data.points.size() - 1; ++i) {
                        draw_list->AddLine(origin + data.points[i], origin + data.points[i + 1], IM_COL32(255, 255, 255, 125), 1.f);
                    }
                }
            }

            //拟合曲线
            if (data.trained) {
                data.curve.clear();
                float startX = data.points[0].x;
                for (int i = 0; i < 1000; ++i){
                    float x = startX + i;
                    float y = model.forward(x, 0.f);
                    data.curve.push_back(ImVec2(x, y));
                }

                for (int i = 0; i < data.curve.size() - 1; ++i) {
                    draw_list->AddLine(origin + data.curve[i], origin + data.curve[i + 1], IM_COL32(255, 0, 0, 255));
                }
            }
        
        }
        ImGui::End();
        EndFrame(window);
    }

    CleanContext(window);
    return 0;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

GLFWwindow* createWindow()
{
    // glfw: initialize and configure
 // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    //msaa
    glfwWindowHint(GLFW_SAMPLES, 4);
    // glfw window creation
 // --------------------
    GLFWwindow* window = glfwCreateWindow(1280, 800, "Subdivision", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return NULL;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync


    //鼠标不可见，内嵌进应用中了
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    //ȫ�ֻص�����
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    //glfwSetCursorPosCallback(window, mouse_callback);
    //glfwSetScrollCallback(window, scroll_callback);

    // glad: load all OpenGL function pointers
    //---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return NULL;
    }
    return window;
}
