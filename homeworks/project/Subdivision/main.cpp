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

void Chaiukin2(CanvasData* data) {
    vector<ImVec2> tmp;
    data->chaiukin2_points = data->points;
    for (int i = 0; i < data->k; ++i) {
        tmp.resize(data->chaiukin2_points.size() * 2 - 2);

        for (int j = 0; j < data->chaiukin2_points.size()-1; ++j) {
            tmp[2*j] = 0.75 * data->chaiukin2_points[j] + 0.25 * data->chaiukin2_points[j + 1];
            tmp[2*j+1] = 0.25 * data->chaiukin2_points[j] + 0.75 * data->chaiukin2_points[j + 1];
        }
        if (data->closure) {
            tmp.push_back(0.25 * data->chaiukin2_points[0] + 0.75 * data->chaiukin2_points[data->chaiukin2_points.size() - 1]);
            tmp.push_back(0.25 * data->chaiukin2_points[data->chaiukin2_points.size() - 1] + 0.75 * data->chaiukin2_points[0]);
        }
        data->chaiukin2_points = tmp;
    }

}
void Chaiukin3(CanvasData* data) {
    vector<ImVec2> tmp;
    data->chaiukin3_points = data->points;
    for (int i = 0; i < data->k; ++i) {
        tmp.resize(data->chaiukin3_points.size() * 2 - 3);
        ImVec2 v0 = 0.5 * (data->chaiukin3_points[0] + data->chaiukin3_points[1]);

        for (int j = 1; j < data->chaiukin3_points.size() - 1; ++j) {
            ImVec2 v1 = 0.5 * (data->chaiukin3_points[j] + data->chaiukin3_points[j + 1]);
            ImVec2 v = 0.5 * (0.5 * (v0 + v1) + data->chaiukin3_points[j]);
            tmp[2 * (j-1)] = v0;
            tmp[2 *j -1] = v;
            v0 = v1;
        }
        *(tmp.end()-1) = v0;
        if (data->closure) {
            //将最后一条线细分为两条线，新增三个顶点
            int n = data->chaiukin3_points.size()-1;
            ImVec2 v1 = 0.5 * (data->chaiukin3_points[n] + data->chaiukin3_points[0]);
            ImVec2 v = 0.5 * (0.5 * (v0 + v1) + data->chaiukin3_points[n]);
            tmp.push_back(v);
            v0 = v1;
            tmp.push_back(v0);
            v1 = 0.5 * (data->chaiukin3_points[0] + data->chaiukin3_points[1]);
            v = 0.5 * (0.5 * (v0 + v1) + data->chaiukin3_points[0]);
            tmp.push_back(v);
        }
        data->chaiukin3_points = tmp;
    }
}
void FourPointSubDivision(CanvasData* data) {
    vector<ImVec2> tmp;
    data->fpsubdiv_points = data->points;
    for (int i = 0; i < data->k; ++i) {
        tmp.resize(data->fpsubdiv_points.size() * 2 - 3);

        tmp[0] = data->fpsubdiv_points[0];
        for (int j = 1; j < data->fpsubdiv_points.size() - 2; ++j) {
            tmp[2 * j - 1] = data->fpsubdiv_points[j];
            ImVec2 v1 = 0.5 * (data->fpsubdiv_points[j] + data->fpsubdiv_points[j + 1]);
            ImVec2 v2 = 0.5 * (data->fpsubdiv_points[j - 1] + data->fpsubdiv_points[j + 2]);
            tmp[2 * j] = v1 + data->alpha * (v1 - v2);
        }
        
        if (!data->closure) {
            *(tmp.end() - 2) = *(data->fpsubdiv_points.end() - 2);
            *(tmp.end() - 1) = *(data->fpsubdiv_points.end() - 1);
            data->fpsubdiv_points = tmp;
        }
        else {
            int n = data->fpsubdiv_points.size() - 1;
            *(tmp.end() - 2) = data->fpsubdiv_points[n - 1];

            //P(n-2) 与P(n-1)之间的细分点
            ImVec2 v1 = 0.5 * (data->fpsubdiv_points[n-1] + data->fpsubdiv_points[n]);
            ImVec2 v2 = 0.5 * (data->fpsubdiv_points[n - 2] + data->fpsubdiv_points[0]);
            *(tmp.end() - 1) = (v1 + data->alpha * (v1 - v2));

            tmp.push_back(data->fpsubdiv_points[n]);
          
            //P(n-1) 与P0之间的细分点
            v1 = 0.5 * (data->fpsubdiv_points[n] + data->fpsubdiv_points[0]);
            v2 = 0.5 * (data->fpsubdiv_points[n-1] + data->fpsubdiv_points[1]);
            tmp.push_back(v1 + data->alpha * (v1 - v2));

            //加入第0个控制点，可是已经加过了
            tmp.push_back(data->fpsubdiv_points[0]);
            
            //P0和P1之间的细分点
            v1 = 0.5 * (data->fpsubdiv_points[0] + data->fpsubdiv_points[1]);
            v2 = 0.5 * (data->fpsubdiv_points[n] + data->fpsubdiv_points[2]);
            tmp.push_back(v1 + data->alpha * (v1 - v2));

            data->fpsubdiv_points.clear();
            data->fpsubdiv_points.insert(data->fpsubdiv_points.cbegin(), (tmp.cbegin() + 1), tmp.cend());
        }
    }
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
            ImGui::InputInt("k", &data.k);
            ImGui::SliderFloat("alpha", &data.alpha, 0, 1);
            ImGui::Checkbox("draw line", &data.draw_line);
            ImGui::Checkbox("closure", &data.closure);
            ImGui::Checkbox("chaiukin2", &data.chaiukin2_draw);
            ImGui::Checkbox("chaiukin3", &data.chaiukin3_draw);
            ImGui::Checkbox("fpsubdiv", &data.fpsubdiv_draw);
        }
        ImGui::End();

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
                }
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                    data.is_draging = false;
                }
            }

            if (data.points.size() >= 2 && data.chaiukin2_draw) {
                Chaiukin2(&data);
            }
            if (data.points.size() >= 3 && data.chaiukin3_draw) {
                Chaiukin3(&data);
            }
            if (data.points.size() >= 4 && data.fpsubdiv_draw) {
                FourPointSubDivision(&data);
            }
            //context popup
            {
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Right))
                    ImGui::OpenPopupContextItem("context");
                if (ImGui::BeginPopup("context")) {
                    ImGui::Selectable("edit", &data.is_editing);
                    if (ImGui::Button("remove one")) {
                        if(data.points.size()> 1){
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
                    if (data.closure) {
                        draw_list->AddLine(origin + data.points[data.points.size() - 1], origin + data.points[0], IM_COL32(255, 255, 255, 125), 1.f);
                    }
                }
            }

            if (data.chaiukin2_draw && data.chaiukin2_points.size() > 1) {
                for (int i = 0; i < data.chaiukin2_points.size() - 1; ++i) {
                    draw_list->AddLine(origin + data.chaiukin2_points[i], origin + data.chaiukin2_points[i + 1], IM_COL32(255, 0, 0, 255),2.f);
                }
                if (data.closure) {
                    draw_list->AddLine(origin + data.chaiukin2_points[data.chaiukin2_points.size()-1], origin + data.chaiukin2_points[0], IM_COL32(255, 0, 0, 255), 2.f);
                }
            }
            if (data.chaiukin3_draw && data.chaiukin3_points.size() > 1) {
                for (int i = 0; i < data.chaiukin3_points.size() - 1; ++i) {
                    draw_list->AddLine(origin + data.chaiukin3_points[i], origin + data.chaiukin3_points[i + 1], IM_COL32(0, 255, 0, 255), 2.f);
                }
                if (data.closure) {
                    draw_list->AddLine(origin + data.chaiukin3_points[data.chaiukin2_points.size()-1], origin + data.chaiukin3_points[0], IM_COL32(0, 255, 0, 255), 2.f);
                }
            }
            if (data.fpsubdiv_draw && data.fpsubdiv_points.size() > 1) {
                for (int i = 0; i < data.fpsubdiv_points.size() - 1; ++i) {
                    draw_list->AddLine(origin + data.fpsubdiv_points[i], origin + data.fpsubdiv_points[i + 1], IM_COL32(0, 0, 255, 255), 2.f);
                }
                if (data.closure) {
                    draw_list->AddLine(origin + data.fpsubdiv_points[data.fpsubdiv_points.size()-1], origin + data.fpsubdiv_points[0], IM_COL32(0, 0, 255, 255), 2.f);
                }
            }
        }
        ImGui::End();
        EndFrame(window);
    }

    CleanContext(window);
    return 0;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height){
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
