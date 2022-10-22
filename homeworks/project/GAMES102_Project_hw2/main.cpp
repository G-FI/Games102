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

struct RBFNetwork_ {
    
    //X, p3, p4, p5, p6, w1, b1, w2, b2
    //p3_grad, p4_grad, p5_grad, p6_grad,...
    RowVectorXf p3, p4,/* p5, p6,*/ w1, b1, w2/*, b2*/;
    RowVectorXf p3_grad, p4_grad, p5_grad, /* p6_grad, */ w1_grad, b1_grad, w2_grad/*, b2_grad*/;
    float p5, p6, b2;
    float  p6_grad, b2_grad;

    float x, y = 0;

    float Pi = 3.1415926;

    float lr = 0.01;
    float epochs = 50;
    //float rbf_size = 100;

    //初始化size，还有值
    void init(int rbf_size){

        default_random_engine e;        
        normal_distribution<float> normal_dist(0, 2.f);
        auto random_init = [&](float dummy) {return normal_dist(e); };

        p3 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);
        p3_grad = RowVectorXf::Zero(rbf_size);

        p4 = RowVectorXf::Zero(rbf_size);
        p4_grad = RowVectorXf::Zero(rbf_size);

        w1 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);;
        w1_grad = RowVectorXf::Zero(rbf_size);

        b1 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);;
        b1_grad = RowVectorXf::Zero(rbf_size);

        w2 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);;
        w2_grad = RowVectorXf::Zero(rbf_size);

        p3 = RowVectorXf::Zero(rbf_size);
        p3_grad = RowVectorXf::Zero(rbf_size);

        b2 = random_init(0);
        b2_grad = 0;

        p5 = 0;
        p5_grad = RowVectorXf::Zero(rbf_size);

        p6 = 0;
        p6_grad = 0;


        /*b2 = RowVectorXf::Zero(1).unaryExpr(random_init);
        b2_grad = RowVectorXf::Zero(1);

      

        p5 = RowVectorXf::Zero(1);
        p5_grad = RowVectorXf::Zero(1);

        p6 = RowVectorXf::Zero(1);
        p6_grad = RowVectorXf::Zero(1);*/
      
    }
    void debug() {
        cout << "p3 = " << p3 << endl;
        cout << "p4 = " << p4 << endl;
        cout << "p5 = " << p5 << endl;
        cout << "p6 = " << p6 << endl;
        cout << "w1 = " << w1 << endl;
        cout << "b1 = " << b1 << endl;
        cout << "w2 = " << w2 << endl;
        cout << "b2 = " << b2 << endl;

        cout << "b1_grad = " << b1_grad << endl;
        cout << "w1_grad = " << w1_grad << endl;
        cout << "w2_grad = " << w2_grad << endl;
        cout << "b2_grad = " << b2_grad << endl;
        cout <<"===============================" << endl << endl;
    }

    float predict(float x) {
        auto gauss_kernel = [=](float x)->float {return exp(-0.5 * x * x) / sqrt(2 * Pi); };

        auto get_p4 = [&](const RowVectorXf& p3)->RowVectorXf {
            RowVectorXf tmp = RowVectorXf::Zero(p3.cols());
            for (int i = 0; i < p3.cols(); ++i) {
                tmp(i) = gauss_kernel(p3(i));
            }
            return tmp;
        };
        p3 = x * w1 + b1;
        p4 = get_p4(p3);
        p5 = p4.dot(w2) + b2;
        return p5; 
    }
    float forward(float x, float y) {
        this->x = x;
        this->y = y;

        auto gauss_kernel = [=](float x)->float {return exp(-0.5 * x * x) / sqrt(2 * Pi); };

        auto get_p4 = [&](const RowVectorXf& p3)->RowVectorXf{
            RowVectorXf tmp = RowVectorXf::Zero(p3.cols());
            for (int i = 0; i < p3.cols(); ++i) {
                tmp(i) = gauss_kernel(p3(i));
            }
            return tmp;
        };
        p3 = x * w1 + b1;
       // cout << "p3 = "<<p3 << endl;
        p4 = get_p4(p3);
        //cout << "p4 = " << p4 << endl;
        p5 = p4.dot(w2) + b2;
        //cout << "p5 = " << p5 << endl;
        p6 = (y - p5) * (y - p5);
        //cout << "p6 = " << p6 << endl;
        return p6;
    }
  
    void backward() {
        p6_grad = -2 * (y - p5);
        //cout << "p6_grad = " << p6_grad << endl;


        p5_grad = p6_grad * w2; //相对于p4的偏导
        //cout << "p5_grad = " << p5_grad << endl;

        w2_grad = p6_grad * p4;
        //cout << "w2_grad = " << w2_grad << endl;

        b2_grad = p6_grad * 1;
        //cout << "b2_grad = " << b2_grad << endl;


        //辅助计算
        RowVectorXf tmp = -0.5 * p3.array() * p3.array();
        tmp = tmp.unaryExpr([](float x) { return exp(x); });
        p4_grad = p5_grad.array() * p3.array() * tmp.array() / (-1 * sqrt(2 * Pi));
        //cout << "p4_grad = " << p4_grad << endl;

        w1_grad = p4_grad * x;
        //cout << "w1_grad = " << w1_grad << endl;
        b1_grad = p4_grad * 1;
        //cout << "b1_grad = " << b1_grad << endl;
    }
    void optimial() {
        w1 -= w1_grad * lr;
        b1 -= b1_grad * lr;
        w2 -= w2_grad * lr;
        b2 -= b2_grad * lr;
    }

    void train(const vector<ImVec2>& data) {
        for (int i = 0; i < epochs; ++i) {
            for (int j = 0; j < data.size(); ++j) {
                forward(data[j].x, data[j].y);
                backward();
                optimial();
            }
        }
    }
    
};

struct RBFNetwork{
    //super parameters
    int epochs_;
    float lr_;

    //forward时缓存的参数，在backward时需要用到
    float x, y, y_pred, loss;

    //Parameters
    ArrayXf W1, W2, B1;
    float b2;

    //对于Parameters的偏导数
    ArrayXf DW1, DW2, DB1;
    float db2;

    //存储forward时变量
    ArrayXf Z, A; //Z = W1 * x + B1,  A = gauss(Z)
    RBFNetwork(int rbf_size, int epochs, float lr) {
        epochs_ = epochs;
        lr_ = lr;
        
        random_device e;
        normal_distribution<float> normal_dist(100, 10.f);
        auto random_init = [&](float dummy) {return normal_dist(e); };

        W1 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);
        W2 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);

        B1 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);

        b2 = random_init(0);
       
        DW1 = RowVectorXf::Zero(rbf_size);
        DW2 = RowVectorXf::Zero(rbf_size);
        DB1 = RowVectorXf::Zero(rbf_size);
        db2 = 0.f;
        
    }

    float forward(float x_input, float y_input) {
        //auto gauss_kernel = [](float x)->float {return exp(-0.5 * x * x) / sqrt(2 * 3.1415926); };
        x = x_input;
        y = y_input;
        
        Z = x * W1 + B1;
        //cout <<"Z = "<< Z << endl;
        A =  ((-0.5f) * Z.pow(2) / sqrt(2 * 3.1415926)).exp();
        //cout << "A = " << A << endl;

        y_pred = (A * W2).sum() + b2;
        //cout << "y_pred = " << y_pred << endl;
        return y_pred;
    }
    void backward() {
        loss = (y - y_pred) * (y - y_pred);

        cout << "loss = " << loss << endl;

        float DYpred = -2 * (y - y_pred);
        
        DW2 = DYpred * A;
        //cout << "DW2 = " << DW2 << endl;

        db2 = DYpred  * 1;
        //cout << "db2 = " << db2 << endl;
        ArrayXf DA = DYpred * W2;
        ArrayXf DZ = DA * A * (-1 * Z);

        DW1 = DZ * x;
        //cout << "DW1 = " << DW1 << endl;
        DB1 = DZ * 1;
        //cout << "DB1 = " << DB1 << endl;
   }

    void optimial() {
        W1 -= DW1 * lr_;
        B1 -= DB1 * lr_;
        W2 -= DW2 * lr_;
        b2 -= db2 * lr_;
    }

    void train(const vector<ImVec2>& data) {


        for (int epoch = 0; epoch < epochs_; ++epoch) {
            for (const ImVec2& point : data) {
                forward(point.x, point.y);
                backward();
                optimial();
            }
        }
    }
};

void main_() {
    using namespace Eigen;
    RBFNetwork model(50, 20, 0.1f);

    for (int i = 0; i < 1000; i++) {
        float loss = model.forward(i, i);
        model.backward();
        model.optimial();
        //model.debug();
        cout << "loss = " << loss << endl;
    }
    /*vector<ImVec2> data{ {1, 1}, {10,10}, {100, 100} };
    model.train(data);  */

    cout << "pridict" << endl;
    for (int x = 1000; x < 1010; x++) {
        float y_pred =  model.forward(x, 0.f);
        cout << "x = " << x << ", y_pred = " << y_pred << endl;
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

    float learning_rate = 0.1f;
    int epochs = 1;
    int rbf_number = 1;
    bool train = false;

    RBFNetwork model(rbf_number, epochs, learning_rate);
 

    //for test
    ImVec2 sampel_pos(0.f, 0.f);
    while (!glfwWindowShouldClose(window)) {
        {//轮询事件
            glfwPollEvents();
            // Start the Dear ImGui frame
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
                glfwSetWindowShouldClose(window, true);
        }

        bool is_hoverd = false;
        bool is_clicked = false;
        bool is_actived = false;

        //逻辑代码
        {
            if (ImGui::Begin("Canvas")) {
                //TODO: 获取window宽度，动态设置windth和space
                ImGui::SetNextItemWidth(100.f);  ImGui::InputFloat("lr", &learning_rate, 0.f, 0.f, "%.5f");

                ImGui::SameLine(0.f, 40.f); ImGui::SetNextItemWidth(100.f); ImGui::InputInt("epochs", &epochs);
       
                ImGui::SameLine(0.f, 40.f); ImGui::SetNextItemWidth(100.f);ImGui::InputInt("Rbf number", &rbf_number);

                ImGui::SameLine(0.f, 40.f); ImGui::SetNextItemWidth(100.f); 
                if (ImGui::Button("train")) {
                    train = true;
                    model.train(data.sample_data);
                    train = false;
                    data.draw_line = true;
                }

                ImGui::SameLine(0.f, 40.f); ImGui::SetNextItemWidth(100.f);
                if (ImGui::Button("clear")) {
                    data.draw_data.clear();
                    data.draw_line = false;
                    data.sample_data.clear();
                }
                   

                // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
                ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
                ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
                if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
                if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
                ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

                // Draw border and background color
                ImGuiIO& io = ImGui::GetIO();
                ImDrawList* draw_list = ImGui::GetWindowDrawList(); //我要绘制的图元
                draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
                draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));


                ImVec2 origin = canvas_p0;
                
                //获取是否点击（使用canvas_sz是，input之下的那些可用区域）
                ImGui::InvisibleButton("canvas", canvas_sz - ImVec2(0.f, 100.f), ImGuiButtonFlags_MouseButtonLeft);
                is_hoverd = ImGui::IsItemHovered();
                is_clicked = ImGui::IsItemClicked();
                if (is_hoverd && is_clicked) {
                    sampel_pos = ImVec2(io.MousePos.x - origin.x, io.MousePos.y - origin.y);
                    data.sample_data.push_back(sampel_pos);
                }

                //数据监控面板
                {
                    if (ImGui::Begin("sample pos")) {
                        ImVec2 pos = io.MousePos;

                        ImGui::Text("cursor pos: (%f, %f), canvas_sz(%f, %f)", canvas_p0.x, canvas_p0.y, canvas_sz.x, canvas_sz.y);
                        ImGui::Text("Mouse pos(%f, %f)", pos.x, pos.y);
                        ImGui::Text("sample pos(%f, %f)", sampel_pos.x, sampel_pos.y);
                    }
                    ImGui::End();
                }

                //绘制
                //绘制采样点
                for (const auto& point : data.sample_data) {
                    draw_list->AddCircleFilled(point + origin, 4.f, IM_COL32(255, 0, 0, 255));
                }
                //绘制预测曲线,                    //起点x坐标为 data.sample_data[0].x - 100.f
                if (data.draw_line) {
                    float start_x = data.sample_data[0].x - 100.f;
                    for (int i = 0; i < data.draw_data_num; ++i) {
                        float x = start_x + i * 1.f;
                        float y = model.forward(x, 0.f);
                        //x,y是window 局部坐标，绘制是相对于screen坐标的
                        data.draw_data.push_back(ImVec2(x, y) + origin);
                    }
                    draw_list->AddPolyline(data.draw_data.data(), data.draw_data_num, IM_COL32(0, 255, 0, 255), ImDrawFlags_None, 1.f);
                    ImGui::Text("first data: (%f, %f)", data.draw_data[0].x, data.draw_data[0].y);
                }
               
            }

            ImGui::End();

        }
        // Rendering
        {
            ImGui::Render();
            int display_w, display_h;
            glfwGetFramebufferSize(window, &display_w, &display_h);
            glViewport(0, 0, display_w, display_h);
            glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        glfwSwapBuffers(window);
    }


        // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
//int main() {
//    /*
//        虚拟环境下，导入torch包等，出现问题,正常环境下应该可以
//    */
//    //Py_SetPythonHome(L"E:\\Anaconda3\\envs\\rl");
//    Py_SetPythonHome(L"E:\\Python");
//
//
//    Py_Initialize();
//
//    //添加当前路径 E:\code\GAMES102\homeworks\project\GAMES102_Project_hw2
//    PyRun_SimpleString("import sys\nsys.path.append('E:/code/GAMES102/homeworks/project/GAMES102_Project_hw2/')");
//
//    PyRun_SimpleString("import torch");
//    PyObject* pMoudle = PyImport_ImportModule("hw2");
//
//
//    PyObject* pFunc = PyObject_GetAttrString(pMoudle, "mytest");
//    PyObject* res = PyObject_CallObject(pFunc, NULL);
//    int result;
//    PyArg_Parse(res, "i", &result);
//
//    Py_Finalize();
//}


int test() {
    const char* glsl_version = "#version 130";
    GLFWwindow* window = createWindow();

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableSetMousePos;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    bool show_some_window = false;
    bool show_another_window = false;

    float color[] = { 0.1f, 0.1f, 0.1f, 1.f };
    int counter = 0;
    while (!glfwWindowShouldClose(window)) {
        //轮询事件
        io.WantCaptureMouse = false;
        io.WantCaptureKeyboard = false;  //imgui 监听事件，不发送给主程序
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        //ImGui::ShowDemoWindow((bool*)1);
        {//手动配置UI
            float slider_value = 0.f;
            if (!ImGui::Begin("h1")) {
                ImGui::Begin("info");
                ImGui::Text("window h1 is hiden");
                ImGui::End();
            }
            ImGui::Text("first text");
            ImGui::Checkbox("first checkbox", &show_some_window);
            {//显示在同一行
                ImGui::SameLine(0.f, 10.f);
                ImGui::Spacing();
               
                ImGui::Checkbox("show another window", &show_another_window);
            }
            ImGui::SliderFloat("float slider", &slider_value, 0.f, 10.f, "%.1f");
            ImGui::ColorEdit4("eidt color", &color[0]);

            if (ImGui::Button("counter")) {
                ++counter;
            }

            //
            ImGui::SameLine();
            ImGui::Text("counter: %d", counter);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();

            if (show_another_window) {
                ImGui::Begin("another window", &show_another_window);
                ImGui::Text("Hello from another window!");
                if (ImGui::Button("click me")) {
                    show_another_window = false;
                }
                ImGui::End();
            }
        }

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(color[0] * color[3], color[1] * color[3], color[2] * color[3], color[3]);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());


        glfwSwapBuffers(window);
    }


    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

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
