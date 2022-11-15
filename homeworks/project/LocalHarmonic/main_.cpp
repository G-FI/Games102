#include<iostream>
#include<algorithm>
#include<vector>

#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "Eigen/Core"
#include "Eigen/Sparse"




using namespace Eigen;
using namespace std;

using T= Eigen::Triplet<int>;


GLFWwindow* CreateWindow();
void NextFrame(GLFWwindow* window);
void EndFrame(GLFWwindow* window);
void CleanContext(GLFWwindow* window);

int DenseMatrixTest()
{
    MatrixXi m = MatrixXi::Random(10, 10);
    Index maxRow, maxCol;
    int maxCoeff = m.maxCoeff(&maxRow, &maxCol);
    cout << m << endl;
    cout << "maxCoeff: " << maxCoeff << "row: " << maxRow << "col: " << maxCol << endl;

    cout << "partial reduction operator" << endl;
    cout << "colwise: " << endl;
    cout << m.colwise().maxCoeff() << endl;
    
    cout << "reshape and resize" << endl;
    //auto m2 = MatrixXi::Random(4, 4);
    auto m2 = Matrix4i::Random();
    cout << "m2: " << m2 << endl;
    auto m2_reshape = m2.reshaped(16, 1);
    cout << "m2_reshped: " << m2_reshape << endl;

    cout << "interator of 2D matrix=============" << endl;
    m = MatrixXi::Random(4, 4);
    cout << "m:" << endl << m << endl;
    for (auto row : m.rowwise())
    {
        sort(row.begin(), row.end());
    }
    cout << "after sort for rowwise:" << endl << m << endl;
    
    return 0;
}
int SparseMatrixTest()
{
    SparseMatrix<int> mat(100, 200);
    cout << "outersize of mat" << mat.outerSize() << endl;
    cout << "innersize of mat" << mat.innerSize() << endl;
    vector<T> tripletList(10);
    tripletList[0] = T(0, 0, 1);
    tripletList[1] = T(5, 5, 2);
    tripletList[2] = T(10, 10, 3);
    tripletList[3] = T(10, 11, 4);
    tripletList[4] = T(20, 20, 5);
    tripletList[5] = T(49, 49, 6);
    tripletList[6] = T(80, 89, 7);
    tripletList[7] = T(99, 99, 8);

    mat.reserve(10);
    //mat.setFromTriplets(tripletList.cbegin(), tripletList.cend());
    for (auto& t : tripletList) {
        mat.insert(t.row(), t.col()) = t.value();
    }

    for (int k = 0; k < mat.outerSize(); ++k) {
        for (SparseMatrix<int>::InnerIterator it(mat, k); it; ++it) {
            cout << it.value() << endl;
            cout << it.index() << endl;
            cout << it.row() << endl;
            cout << it.col() << endl;
            cout << "=====================" << endl << endl;
        }
    }

    

    return 0;
}

int main_() 
{
    GLFWwindow* window = CreateWindow();

    while (!glfwWindowShouldClose(window))
    {
        NextFrame(window);

        EndFrame(window);
    }

    CleanContext(window);
    return 0;
}
void NextFrame(GLFWwindow* window) {
    //轮询事件
    glfwPollEvents();
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}
void EndFrame(GLFWwindow* window) {
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glfwSwapBuffers(window);
}
void CleanContext(GLFWwindow* window) {
    // Cleanup

    glfwDestroyWindow(window);
    glfwTerminate();
}



void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

GLFWwindow* CreateWindow()
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
