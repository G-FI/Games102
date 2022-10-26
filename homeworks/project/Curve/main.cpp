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

void Monitor(bool *p_open, CanvasData* data) {
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

void Bezier(CanvasData* data) {
    //de  Casteljau algorithm
    //n个控制点生成，n-1阶的Bezier
    int n = data->points.size();
    data->bezier_curve.clear();
    vector<ImVec2> tmp1;
    vector<ImVec2> tmp2;
    for (int num = 0; n < data->num; ++num) {
        tmp1 = data->points;
        tmp2.resize(n);

        float t = (float)num / data->num;
       
        //n个控制点，需要迭代n-1次，生成一个曲线再t处的点
        for (int r = 1; r< n; ++r) {
            //总共n个控制点，经r次迭代之后只剩n-r个控制点
            for (int i = 0; i < n - r; ++i) {
                tmp2[i] = (1 - t) * tmp1[i] + t * tmp1[i+1];
            }
            swap(tmp1, tmp2);
        }
        data->bezier_curve.push_back(tmp1[0]);
    }
    data->bezier_draw = true;
}

float Nbase(int i, int k, float t, float dt) {
    if (k == 1) {
        if (t >= dt * i && t < dt * (i + 1)) {
            return 1;
        }
        else {
            return 0;
        }
    }
    float first = Nbase(i, k - 1, t, dt);
    float second = Nbase(i + 1, k - 1, t, dt);
    first = first * (t - i * dt) / ((k - 1) * dt);
    second = second * ((i + k) * dt - t) / ((k - 1) * dt);
    return first + second;
}
void Bspline(CanvasData* data) {
    //绘制k阶B样条曲线，至少需要k个顶点
    data->b_curve.clear();

    if (data->points.size() < data->k) {
        return;
    }
    int k = data->k;
    float dt = (float)1 / (data->points.size() -1 + data->k);
    for (int i = k - 1; i < data->points.size(); ++i) {
        //对在[ti..ti+1]区间的曲线进行采样
        for (int m = 0; m < data->num; ++m) {
             ImVec2 p(0.f, 0.f);
            float tot = 0.f;
            //这里的采样是每一个[ti,ti+1]段采num个样，t的增量是相对全局的，因此还需要除以曲线的段数[n-1+k]
            float t = i * dt + (float)m/(data->num * (data->points.size() - 1 + data->k));
            for (int j = i + 1 - k; j <= i; ++j) {
                 p = p + Nbase(j, k, t, dt) * data->points[j];
                float tmp = Nbase(j, k, t, dt);
                tot += tmp;
            }
            data->b_curve.push_back(p);
        }

    }

    data->b_draw = true;
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
                draw_list->AddRect(origin, origin+canvas_sz, IM_COL32(255, 255, 255, 255));

                //获取是否点击（使用canvas_sz是，input之下的那些可用区域）
                ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft);
                ImVec2 mouse_in_canvas = io.MousePos - origin;

                if (!data.editing && ImGui::IsItemHovered() && ImGui::IsItemClicked()) {
                    data.points.push_back(mouse_in_canvas);
                }
                //交互修改控制点
                if (data.editing && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    for(int i =0 ;i < data.points.size(); ++i){
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

                //生成Bezeir曲线
                /*if (data.points.size() > 1) {
                    Bezier(&data);
                }*/
                if (data.points.size() > 1) {
                    Bspline(&data);
                }
                //context popup
                {
                    if (ImGui::IsMouseReleased(ImGuiMouseButton_Right))
                        ImGui::OpenPopupContextItem("context");
                    if (ImGui::BeginPopup("context")) {
                        ImGui::Selectable("monitor", &data.monitor_open);
                        ImGui::Selectable("edit", &data.editing);

                        if (ImGui::Button("remove one")) {
                            data.points.pop_back();
                        }
                        if (ImGui::Button("remove all")) {
                            data.Reset();
                        }
                        ImGui::InputInt("k", &data.k);

                        ImGui::EndPopup();
                    }

                }

                //绘制
                for (const auto& point : data.points) {
                    draw_list->AddCircleFilled(point + origin, 4.f, IM_COL32(255, 0, 0, 255));
                }

                if (data.bezier_draw && data.bezier_curve.size() > 1) {
                    for (int i = 0; i < data.bezier_curve.size()-1; ++i) {
                        draw_list->AddLine(ImVec2(origin + data.bezier_curve[i]), ImVec2(origin + data.bezier_curve[i + 1]), IM_COL32(255, 255, 0, 255));
                    }
                }
                if (data.b_draw && data.b_curve.size() > 1) {
                    for (int i = 0; i < data.b_curve.size() - 1; ++i) {
                            draw_list->AddLine(ImVec2(origin + data.b_curve[i]), ImVec2(origin + data.b_curve[i + 1]), IM_COL32(255, 255, 0, 255));
                    }
                    /*for (const auto& point : data.b_curve) {
                        draw_list->AddCircleFilled(point + origin, 1.f, IM_COL32(255, 255, 0, 255));
                    }*/
                }

                if (data.editing && (data.bezier_draw || data.b_draw)) {
                    for (int i = 0; i < data.points.size() - 1; ++i) {
                        draw_list->AddLine((origin + data.points[i]), (origin + data.points[i+1]), IM_COL32(0, 127, 0, 127), 0.5);
                    }
                }

               /* if (data.points.size() == 4) {
                    draw_list->AddBezierCurve((origin + data.points[0]), (origin + data.points[1]), (origin + data.points[2]), (origin + data.points[3]), IM_COL32(0, 127, 0, 127), 0.5);
                }*/
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
    GLFWwindow* window = glfwCreateWindow(1280, 800, "LearnOpenGL", NULL, NULL);
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
