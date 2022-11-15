#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

//surface mesh
#include "CGAL/Simple_cartesian.h"
#include "CGAL/Surface_mesh.h"

//读入网格库
#include "CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h"

//libigl显示
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <Eigen/Sparse>
#include<Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>



typedef igl::opengl::glfw::Viewer Viewer;

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef Kernel::Vector_3                    Vector_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Mesh::Halfedge_index                HalfedgeDescirptor;
typedef Mesh::Vertex_index                  VertexDescriptor;
typedef Mesh::Face_index                    FaceDescriptor;


class DenoiseData {
public:
    std::vector<Mesh> meshs;
    std::vector<Mesh> copy_meshs;
    bool toggle = true;
    int idx = 0;

    void ReadMesh(const std::string mesh_path)
    {
        meshs.emplace_back();
        CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(mesh_path, meshs[meshs.size() - 1]);
        copy_meshs.push_back(meshs[meshs.size() - 1]);
    }

    void Restore(int idx)
    {
        meshs[idx] = copy_meshs[idx];
    }

};

class CullPlane 
{
public:
    std::vector<Point_2> border;
};

DenoiseData data;
Eigen::MatrixXd V;
Eigen::MatrixXi F;


void HalfEdge2Metrix(Mesh& m, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    V.resize(m.number_of_vertices(), 3);
    F.resize(m.number_of_faces(), 3);
    //顺序加入顶点
    int i = 0;
    for (VertexDescriptor vd : m.vertices())
    {
        auto point = m.point(vd);
        V.row(i) = Eigen::RowVector3d(point.x(), point.y(), point.z());
        ++i;
    }

    std::cout << "shape of V: " << V.rows() << " , " << V.cols() << std::endl;
    std::cout << "shape of F: " << F.rows() << " , " << F.cols() << std::endl;

    //顺序加入face
    i = 0;
    for (FaceDescriptor fd : m.faces())
    {
        //std::cout << i << std::endl;
        int j = 0;
        for (VertexDescriptor vd : m.vertices_around_face(m.halfedge(fd)))
        {

            F(i, j) = vd.idx();
            ++j;
        }
        ++i;
    }

}

//从一个边界半边开始，按顺序遍历下一个半边，添加对应顶点
std::vector<VertexDescriptor> SortBorderVertex(Mesh& m)
{
    std::vector<VertexDescriptor> vbs;
    VertexDescriptor v_start;
    vbs.reserve(m.number_of_vertices() / 5);
    
    HalfedgeDescirptor he;
    //找到边界点
    for (auto tmp : m.halfedges()) {
        if (m.is_border(tmp)) {
            he = tmp;
            break;
        }
    }
    auto done(he);

    do {
        vbs.push_back(m.target(he));
        he = m.next(he);
        if (!m.is_border(he)) {
            do {
                he = m.next(m.opposite(he));
            } while (!m.is_border(he));
        }
 
    } while (he != done);

    return vbs;
}



Eigen::MatrixXd Parameterize(Mesh& m)
{
    auto Cot = [](Vector_3 v1, Vector_3 v2)->double {
        auto v2_norm = v2 / sqrt(v2.squared_length());

        auto p = v1 * v2_norm;//v1在v2上投影长度
        auto pv = v2_norm * p;

        auto a = sqrt((v1 - pv).squared_length());
        /*
              v1
             /| \
            / |a \
           /__|___\v2
            p
        */
        return p / a;
    };

    Eigen::SparseMatrix<double> L(m.number_of_vertices(), m.number_of_vertices());
    Eigen::SparseMatrix<double> delta(m.number_of_vertices(), 3);
    Eigen::MatrixXd X(m.number_of_vertices(), 3);

    std::vector<Eigen::Triplet<double>> tripL;
    std::vector<Eigen::Triplet<double>> tripDelta;
    tripL.reserve(m.number_of_halfedges());
    tripDelta.reserve(m.number_of_vertices());

    auto s_border = SortBorderVertex(m);
    //构造delta矩阵
    int bv_number = (int)s_border.size();
    double dt = (double)2*3.1415926 / bv_number;
    double theta = 0.0;
    for (auto vb : s_border)
    {
        auto u= cos(theta);
        auto v = sin(theta);

        tripDelta.emplace_back(vb.idx(), 0, u);
        tripDelta.emplace_back(vb.idx(), 1, v);
        theta += dt;
    }
    
    
    //构建稀疏矩阵
    for (VertexDescriptor v : m.vertices())
    {
        Point_3 p_i = m.point(v);

       
        if (!m.is_border(v))
        {
            auto he = m.halfedge(v), done(he);
            double total_weight = 0.0;
            do {
                //v是he的target， vj是he的source
                auto vj = m.source(he);
                auto v_alpha = m.target(m.next(he));
                auto v_beta = m.source(m.prev(m.opposite(he)));
                auto p_j = m.point(vj);
                auto p_alpha = m.point(v_alpha);
                auto p_beta = m.point(v_beta);

                auto wij = Cot(p_alpha - p_i, p_alpha - p_j) + Cot(p_beta - p_i, p_beta - p_j);
                total_weight += wij;

                tripL.emplace_back(v.idx(), vj.idx(), -wij);
                he = m.opposite(m.next(he));
            } while (he != done);
            tripL.emplace_back(v.idx(), v.idx(), total_weight);
        }
        else
        {
            tripL.emplace_back(v.idx(), v.idx(), 1);
        }
    }

    L.setFromTriplets(tripL.begin(), tripL.end());
    delta.setFromTriplets(tripDelta.begin(), tripDelta.end());
    //求解
    Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(L);

    solver.factorize(L);

    X = solver.solve(delta);

    Eigen::MatrixXd V_uv(m.number_of_vertices(), 2);
    for (auto v : m.vertices())
    {
        V_uv(v.idx(), 0) = X(v.idx(), 0);
        V_uv(v.idx(), 1) = X(v.idx(), 1);
    }
    return V_uv;
}

void HarmonicCotangentWeight(Mesh& m)
{
    auto Cot = [](Vector_3 v1, Vector_3 v2)->double {
        auto v2_norm = v2 / sqrt(v2.squared_length());

        auto p = v1 * v2_norm;//v1在v2上投影长度
        auto pv = v2_norm * p;

        auto a = sqrt((v1 - pv).squared_length());
        /*
              v1
             /| \
            / |a \
           /__|___\v2
            p
        */
        return p / a;
    };

    Eigen::SparseMatrix<double> L(m.number_of_vertices(), m.number_of_vertices());
    Eigen::SparseMatrix<double> delta(m.number_of_vertices(), 3);
    Eigen::MatrixXd X(m.number_of_vertices(), 3);

    std::vector<Eigen::Triplet<double>> tripL;
    std::vector<Eigen::Triplet<double>> tripDelta;
    tripL.reserve(m.number_of_halfedges());
    tripDelta.reserve(m.number_of_vertices());

    //构建稀疏矩阵
    for (VertexDescriptor v : m.vertices())
    {
        Point_3 p_i = m.point(v);

        if (!m.is_border(v))
        {
            auto he = m.halfedge(v), done(he);
            double total_weight = 0.0;
            do {
                //v是he的target， vj是he的source
                auto vj = m.source(he);
                auto v_alpha = m.target(m.next(he));
                auto v_beta = m.source(m.prev(m.opposite(he)));
                auto p_j = m.point(vj);
                auto p_alpha = m.point(v_alpha);
                auto p_beta = m.point(v_beta);

                auto wij = Cot(p_alpha - p_i, p_alpha - p_j) + Cot(p_beta - p_i, p_beta - p_j);
                total_weight += wij;

                tripL.emplace_back(v.idx(), vj.idx(), -wij);
                he = m.opposite(m.next(he));
            } while (he != done);
            tripL.emplace_back(v.idx(), v.idx(), total_weight);
        }
        else
        {
            tripL.emplace_back(v.idx(), v.idx(), 1);
            auto x = p_i.x();
            auto y = p_i.y();
            auto z = p_i.z();
            tripDelta.emplace_back(v.idx(), 0, x);
            tripDelta.emplace_back(v.idx(), 1, y);
            tripDelta.emplace_back(v.idx(), 2, z);
        }
    }

    L.setFromTriplets(tripL.begin(), tripL.end());
    delta.setFromTriplets(tripDelta.begin(), tripDelta.end());
    //求解
    Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(L);
    solver.factorize(L);

    X = solver.solve(delta);

    for (auto v : m.vertices())
    {
        auto x = X(v.idx(), 0);
        auto y = X(v.idx(), 1);
        auto z = X(v.idx(), 2);

        m.point(v) = Point_3(X(v.idx(), 0), X(v.idx(), 1), X(v.idx(), 2));
    }
}

void HarmonicUniformWeight(Mesh& m)
{
    auto Cot = [](Vector_3 v1, Vector_3 v2)->double {
        auto v2_norm = v2 / sqrt(v2.squared_length());

        auto p = v1 * v2_norm;//v1在v2上投影长度
        auto pv = v2_norm * p;

        auto a = sqrt((v1 - pv).squared_length());
        /*
              v1
             /| \
            / |a \
           /__|___\v2
            p
        */
        return p / a;
    };

    Eigen::SparseMatrix<double> L(m.number_of_vertices(), m.number_of_vertices());
    Eigen::SparseMatrix<double> delta(m.number_of_vertices(), 3);
    Eigen::MatrixXd X(m.number_of_vertices(), 3);

    std::vector<Eigen::Triplet<int>> tripL;
    std::vector<Eigen::Triplet<double>> tripDelta;
    tripL.reserve(m.number_of_halfedges());
    tripDelta.reserve(m.number_of_vertices());

    //构建稀疏矩阵
    for (VertexDescriptor v : m.vertices())
    {
        Point_3 p_i = m.point(v);
        X(v.idx(), 0) = p_i.x();
        X(v.idx(), 1) = p_i.y();
        X(v.idx(), 2) = p_i.z();


        if (!m.is_border(v))
        {
            auto he = m.halfedge(v), done(he);
            /*double total_weight = 0.0;
            do {
                //v是he的target， vj是he的source
                auto vj = m.source(he);
                auto v_alpha = m.target(m.next(he));
                auto v_beta = m.source(m.prev(m.opposite(he)));
                auto p_j = m.point(vj);
                auto p_alpha = m.point(v_alpha);
                auto p_beta = m.point(v_beta);

                total_weight += Cot(p_alpha - p_i, p_alpha - p_j) + Cot(p_beta - p_i, p_beta - p_j);

                he = m.opposite(m.next(he));
            } while (he != done);
            */
            int n = 0;

            do {
                ++n;
                auto vj = m.source(he);

                tripL.emplace_back(v.idx(), vj.idx(), -1);

                he = m.opposite(m.next(he));
            } while (he != done);
            tripL.emplace_back(v.idx(), v.idx(), n);
        }
        else
        {
            tripL.emplace_back(v.idx(), v.idx(), 1);
            auto x = p_i.x();
            auto y = p_i.y();
            auto z = p_i.z();
            tripDelta.emplace_back(v.idx(), 0, x);
            tripDelta.emplace_back(v.idx(), 1, y);
            tripDelta.emplace_back(v.idx(), 2, z);
        }
    }

    L.setFromTriplets(tripL.begin(), tripL.end());
    delta.setFromTriplets(tripDelta.begin(), tripDelta.end());
    //L.makeCompressed();
    //求解
    Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(L);
    solver.factorize(L);

    X = solver.solve(delta);

    for (auto v : m.vertices())
    {
        auto x = X(v.idx(), 0);
        auto y = X(v.idx(), 1);
        auto z = X(v.idx(), 2);

        m.point(v) = Point_3(X(v.idx(), 0), X(v.idx(), 1), X(v.idx(), 2));
    }
}

//igl viewer koyboard callback
//1, 2 global harmonic
//3, 4 Parameterize
bool KeyDown(Viewer& viewer, unsigned char key, int  modifier)
{
    if (key == '1') 
    {
        //std::cout << "Denoise callback" << std::endl;
        HarmonicCotangentWeight(data.meshs[data.idx]);
        HalfEdge2Metrix(data.meshs[data.idx], V, F);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.core().align_camera_center(V, F);
    }
    else if (key == '2')
    {
        data.Restore(data.idx);
        HalfEdge2Metrix(data.meshs[data.idx], V, F);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.core().align_camera_center(V, F);

    }
    else if (key == '3')
    {
        //显示2D uv图
        auto V_uv = Parameterize(data.meshs[data.idx]);

        V_uv *= 5;
        viewer.data().set_vertices(V_uv);
        viewer.data().set_uv(V_uv);
        viewer.core().align_camera_center(V_uv, F);
    }
    else if (key == '4')
    {
        //显3Dmodel
        auto V_uv =Parameterize(data.meshs[data.idx]);        
        HalfEdge2Metrix(data.meshs[data.idx], V, F);

        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.data().show_texture = true;

        V_uv *= 5;
        viewer.data().set_uv(V_uv);
        viewer.core().align_camera_center(V, F);
    }
    
    return false;
}



int main()
{
    std::string bunny_head = "D:/Games102/homeworks/project/assets/BunnyHead.off";
    std::string cat_head = "D:/Games102/homeworks/project/assets/CatHead.off";
    std::string david = "D:/Games102/homeworks/project/assets/David.off";
    std::string Nefer = "D:/Games102/homeworks/project/assets/NefertitiFace.off";
    int w = 800, h = 600;
    data.ReadMesh(bunny_head);
    data.ReadMesh(cat_head);
    data.ReadMesh(david);
    data.ReadMesh(Nefer);


    Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]() {
        menu.draw_viewer_menu();
    };

    std::vector<std::string> meshName = { "bunny_head", "cat_head","David","NefertitiFace" };
   
    menu.callback_draw_custom_window = [&]() {
        if (ImGui::Begin("args"))
        {
            if (ImGui::Combo("meshIdx", &data.idx, meshName))
            {
                HalfEdge2Metrix(data.meshs[data.idx], V, F); 
                viewer.data().clear();
                viewer.data().set_mesh(V, F);

                viewer.core().align_camera_center(V, F);
            }
        }
        ImGui::End();

    };

    //默认使用第一个网格
    HalfEdge2Metrix(data.meshs[data.idx], V, F);
    viewer.data().set_mesh(V, F);

    viewer.callback_key_down = KeyDown;
    viewer.launch();
}