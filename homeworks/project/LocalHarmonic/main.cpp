#include <iostream>
#include <vector>
#include <string>

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


typedef igl::opengl::glfw::Viewer Viewer;

typedef CGAL::Simple_cartesian<double>      Kernel;
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
    int idx = 0;
    double lambda = 1e-6;
    int k = 1;

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

DenoiseData data;
Eigen::MatrixXd V;
Eigen::MatrixXi F;

int test_range_iterator()
{
    Mesh m;
    // u            x
    // +------------+
    // |            |
    // |            |
    // |      f     |
    // |            |
    // |            |
    // +------------+
    // v            w
    // Add the points as vertices
    VertexDescriptor u = m.add_vertex(Kernel::Point_3(0, 1, 0));
    VertexDescriptor v = m.add_vertex(Kernel::Point_3(0, 0, 0));
    VertexDescriptor w = m.add_vertex(Kernel::Point_3(1, 0, 0));
    VertexDescriptor x = m.add_vertex(Kernel::Point_3(1, 1, 0));
    /* FaceDescriptor f = */ m.add_face(u, v, w, x);
    {
        std::cout << "all vertices " << std::endl;
        // The vertex iterator type is a nested type of the Vertex_range
        Mesh::Vertex_range::iterator  vb, ve;
        Mesh::Vertex_range r = m.vertices();
        // The iterators can be accessed through the C++ range API
        vb = r.begin();
        ve = r.end();
        // or the boost Range API
        vb = boost::begin(r);
        ve = boost::end(r);
        // or with boost::tie, as the CGAL range derives from std::pair
        for (boost::tie(vb, ve) = m.vertices(); vb != ve; ++vb) {
            std::cout << *vb << std::endl;
        }
        // Instead of the classical for loop one can use
        // the boost macro for a range
        for (VertexDescriptor vd : m.vertices()) {
            
            std::cout << vd.idx() << std::endl;
        }
        // or the C++11 for loop. Note that there is a ':' and not a ',' as in BOOST_FOREACH
        for (VertexDescriptor vd : m.vertices()) {
            std::cout << vd.idx() << std::endl;
        }
    }
    return 0;
}
int test()
{
    Mesh m;
    // Add the points as vertices
    VertexDescriptor u = m.add_vertex(Kernel::Point_3(0, 1, 0));
    VertexDescriptor v = m.add_vertex(Kernel::Point_3(0, 0, 0));
    VertexDescriptor w = m.add_vertex(Kernel::Point_3(1, 1, 0));
    VertexDescriptor x = m.add_vertex(Kernel::Point_3(1, 0, 0));
    m.add_face(u, v, w);
    FaceDescriptor f = m.add_face(u, v, x);
    if (f == Mesh::null_face())
    {
        std::cerr << "The face could not be added because of an orientation error." << std::endl;
        f = m.add_face(u, x, v);
        assert(f != Mesh::null_face());
    }

    auto vit = m.vertices();
    
    return 0;
}




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

double GetAreaMix(Mesh& m, VertexDescriptor vp)
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
    auto GetArea = [](Vector_3 v1, Vector_3 v2) {
        auto v2_norm = v2 / sqrt(v2.squared_length());

        auto p = v1 * v2_norm;//v1在v2上投影长度
        auto pv = v2_norm * p;
        auto a = (v1 - pv).squared_length();

        return  0.5 * sqrt(v2.squared_length() * a);
    };

    auto GetVoronoiArea = [&Cot](Point_3 P, Point_3 Q, Point_3 R) {
        auto CotQ = Cot(P - Q, R - Q);
        auto CotR = Cot(P - R, Q - R);
        return (Q - P).squared_length() * CotR + (R - P).squared_length() * CotQ;
    };
    auto P = m.point(vp);
    //v关联的是halfedge的target， target(halfedge(v)) == v
    HalfedgeDescirptor he = m.halfedge(vp), done(he);
    double area = 0.0;

    do {
        //v是he的target， vj是he的source
        auto vr = m.source(he);
        auto vq = m.target(m.next(he));
        auto Q = m.point(vq);
        auto R = m.point(vr);

        //P点夹边向量
        auto v1 = R - P;
        auto v2 = Q - P;
        //Q点夹边向量
        auto v3 = -v2;
        auto v4 = R - Q;
        //R点夹边向量
        auto v5 = -v1;
        auto v6 = -v4;

        if (v1 * v2 > 0 && v2 * v3 > 0 && v5 * v6 > 0)
        {//使用vronoio面积
            area += GetVoronoiArea(P, Q, R);
        }
        else
        {
            if (v1 * v2 < 0)
            {
                area += GetArea(v1,v2) / 2.0;
            }
            else
            {
                area += GetArea(v1,v2) / 4.0;
            }
        }

        he = m.opposite(m.next(he));
    } while (he != done);

    return area;
}

//return K (也叫Laplac-Beltrami算子)
Vector_3 GetMeanCurvatureOpertor(Mesh& m, VertexDescriptor v)
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
    auto p_i = m.point(v);
    //v关联的是halfedge的target， target(halfedge(v)) == v
    HalfedgeDescirptor he = m.halfedge(v), done(he);

    Vector_3 K(0.0, 0.0, 0.0);
    double area = GetAreaMix(m, v);
    
    double total_weight = 0.0;
    Point_3 avg_p(0, 0, 0);
    do {
        //v是he的target， vj是he的source
        auto vj = m.source(he);
        auto v_alpha = m.target(m.next(he));
        auto v_beta = m.source(m.prev(m.opposite(he)));
        auto p_j = m.point(vj);
        auto p_alpha = m.point(v_alpha);
        auto p_beta = m.point(v_beta);

        auto cot_alpha = Cot(p_i - p_alpha,p_j - p_alpha);
        auto cot_beta = Cot(p_i - p_beta, p_j - p_beta);

        auto w = (cot_alpha + cot_beta);
        
        //此处直接使用Voronoi面积，论文上说误差边界太大，所以使用上面的GetAreaMix
        //area += (cot_alpha + cot_beta) * (p_i - p_j).squared_length();
        avg_p = Point_3(w*p_j.x() + avg_p.x(), w * p_j.y() + avg_p.y(), w * p_j.z() + avg_p.z());
        //K += (cot_alpha + cot_beta) * (p_i - p_j);
        total_weight += w;
        he = m.opposite(m.next(he));

    } while (he != done);
    avg_p = Point_3(avg_p.x()/total_weight,avg_p.y()/total_weight, avg_p.z()/total_weight);

    //使用 cot(alpha) + cot(beta) / total_weight 作为权重，权重满足和为1，那平均值点应该在1-ring neibour所构成的平面中
    //使用这个作为法向更新顶点的话，按理不应该生成反向的凸起啊？
    //return  K / total_weight;
     
    Vector_3 normal = p_i - avg_p;
    return normal;
    //return K * 1.0 / (2 * area);
    
}

void Denoise(Mesh& m)
{
    //保存点
    std::vector<Point_3> newPoints(m.number_of_vertices());

    for (int n = 0; n < data.k; ++n)
    {
        //计算新的顶点
        for (VertexDescriptor v : m.vertices())
        {
            int i = v.idx();
            if (!m.is_border(v) /*&& i <=5*/)
            {
                Vector_3 K = GetMeanCurvatureOpertor(m, v);
                Vector_3 Hn = K;// / 2.0;
                /* auto oldPoint = m.point(v);
                 auto newPoint = oldPoint - lambda * Hn;*/
                newPoints[i] = m.point(v) - data.lambda * Hn;

            }
            else
            {
                newPoints[i] = m.point(v);
            }
        }

        //更新顶点
        for (VertexDescriptor v : m.vertices())
        {
            m.point(v) = newPoints[v.idx()];
        }
    }
}

//void Denoise(Mesh& m)
//{
//    //保存点
//    std::vector<Point_3> newPoints(m.number_of_vertices());
//
//    for (int k = 0; k < data.k; ++k)
//    {
//        //计算新的顶点
//        for (VertexDescriptor v : m.vertices())
//        {
//            int i = v.idx();
//            if (!m.is_border(v) /*&& i <=5*/)
//            {
//                
//                auto he = m.halfedge(v), done(he);
//                int n = 0;
//                Point_3 avg_p(0, 0, 0);
//                do {
//                    ++n;
//                    auto vj = m.source(he);
//                    auto p = m.point(vj);
//                    avg_p = Point_3(p.x() + avg_p.x(), p.y() + avg_p.y(), p.z() + avg_p.z());
//                    he = m.opposite(m.next(he));
//                } while (he != done);
//
//                avg_p = Point_3(avg_p.x() / n, avg_p.y() / n, avg_p.z() / n);
//
//                auto p0 = m.point(v);
//                Vector_3 normal = p0 - avg_p;
//                
//                newPoints[i] = p0 - data.lambda * normal;
//            }
//            else
//            {
//                newPoints[i] = m.point(v);
//            }
//        }
//
//        //更新顶点
//        for (VertexDescriptor v : m.vertices())
//        {
//            m.point(v) = newPoints[v.idx()];
//        }
//    }
//}

//igl viewer koyboard callback
bool KeyDown(Viewer& viewer, unsigned char key, int  modifier)
{
    if (key == '1') //进行迭代去噪
    {
        //std::cout << "Denoise callback" << std::endl;
        Denoise(data.meshs[data.idx]);
        HalfEdge2Metrix(data.meshs[data.idx], V, F);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.core().align_camera_center(V, F);
    }
    else if (key == '2')
    {
        std::cout << "restore mesh" << std::endl;

        data.Restore(data.idx);
        HalfEdge2Metrix(data.meshs[data.idx], V, F);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
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

    std::vector<std::string> meshName = { "bunny_head", "cat_head","David","NefertitiFace"};

    double a1, b1, c1, a2, b2, c2, a3, b3, c3 = 0;
    auto Cot = [](Vector_3 v1, Vector_3 v2)->double {
        auto v2_norm = v2 / sqrt(v2.squared_length());
        std::cout << "v2_norm = " << v2_norm << std::endl;

        auto p = v1 * v2_norm;//v1在v2上投影长度
        auto pv = v2_norm * p;
        std::cout << "projection vector = " << pv << std::endl;
        std::cout << "length of projection = " << p << std::endl;
        
        auto av = (v1 - pv);
        std::cout << "height vector = " << av << std::endl;
        
        auto a = sqrt(av.squared_length());
        std::cout << "length of height vector = " << a << std::endl;
        /*
              v1
             /| \
            / |a \
           /__|___\v2
            p
        */
        std::cout << "cot = " << p / a << std::endl;
        return p / a;
    };

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

            ImGui::InputDouble("lambda", &data.lambda, 0.00000001, 0.0, "%.8f");
            ImGui::SliderInt("iterator num", &data.k, 1, 100);
        }
        ImGui::End();
        
        if (ImGui::Begin("test")) {
            ImGui::InputDouble("a1", &a1);
            ImGui::InputDouble("b1", &b1);
            ImGui::InputDouble("c1", &c1);
            ImGui::InputDouble("a2", &a2);
            ImGui::InputDouble("b2", &b2);
            ImGui::InputDouble("c2", &c2);
            ImGui::InputDouble("a3", &a3);
            ImGui::InputDouble("b3", &b3);
            ImGui::InputDouble("c3", &c3);

            if (ImGui::Button("calculate")) {
                Point_3 p1(a1, b1, c1);
                Point_3 p2(a2, b2, c2);
                Point_3 p3(a3, b3, c3);
                Cot(p2 - p1, p3 - p1);
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