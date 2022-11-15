#include "igl/readOBJ.h"
#include "igl/writeOFF.h"
//#include "igl/opengl/glfw/Viewer.h"

#include "Eigen/Core"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Polyhedron_3.h"
#include "CGAL/IO/Polyhedron_iostream.h"
#include "CGAL/HalfedgeDS_vector.h"
//#include "CGAL/draw_polyhedron.h"
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>     Polyhedron;
//#include "igl/opengl/glfw/Viewer.h"

#include <string>
std::string assetsDir = "D:/Games102/homeworks/project/assets";
std::string BunnyPath = "D:/Games102/homeworks/project/assets/BunnyHead.off";

void Obj2Off()
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	std::string bunny = assetsDir + "/Bunny_head.obj";
	std::string ball = assetsDir + "/Balls.obj";
	std::string cat = assetsDir + "/Cat_head.obj";
	std::string David328 = assetsDir + "/David328.obj";
	std::string Nefertiti = assetsDir + "/Nefertiti_face.obj";
	
	igl::readOBJ(bunny, V, F);
	igl::writeOFF(assetsDir + "/BunnyHead.off", V, F);

	igl::readOBJ(Nefertiti, V, F);
	igl::writeOFF(assetsDir + "/NefertitiFace.off", V, F);

	igl::readOBJ(ball, V, F);
	igl::writeOFF(assetsDir + "/Balls.off", V, F);

	igl::readOBJ(cat, V, F);
	igl::writeOFF(assetsDir + "/CatHead.off", V, F);

	igl::readOBJ(David328, V, F);
	igl::writeOFF(assetsDir + "/David.off", V, F);

}


//将半边结构中的数据转换到矩阵表示，用与Igl显示, TODO无法找到顶点对应的索引，也就无法给F添加顶点索引
void HalfEdge2Matrix(const Polyhedron& Mesh, Eigen::MatrixXd &V, Eigen::MatrixXi &F){
	auto he = Mesh.halfedges_begin();
	
}


int main_use_ployhedron()
{
	Polyhedron Bunny;
	std::ifstream in(CGAL::data_file_path(BunnyPath));

	CGAL::IO::set_pretty_mode(std::cout);
	in >> Bunny;

	std::cout << "sizof border edges: " << Bunny.size_of_border_edges() << std::endl;
	std::cout << "sizof border halfedges: " << Bunny.size_of_border_halfedges() << std::endl;
	std::cout << "sizof facets: " << Bunny.size_of_facets() << std::endl;
	std::cout << "sizof halfedges: " << Bunny.size_of_halfedges() << std::endl;
	std::cout << "sizof vertices: " << Bunny.size_of_vertices() << std::endl;
	for (auto bhe = Bunny.border_halfedges_begin(); bhe != nullptr; ++bhe) {
		std::cout << "he in border: " << bhe->vertex()->point() << std::endl;
	}

	for (auto he = Bunny.halfedges_begin(); he != Bunny.halfedges_end(); ++he)
	{
		/*
		{auto f = he->facet();
		
		if (f == nullptr) {
			std::cout << "boder edge" << std::endl;
			std::cout << "he->is_border():" << he->is_border() << std::endl;
			std::cout << "incident vertex postion: " << he->vertex()->point() << std::endl;
		}
		}*/
		if (he->is_border()) {
			//该半边位于border
			//记录位置
		}
	}
	return 0;
}