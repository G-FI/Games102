#include <iostream>
#include <algorithm>

#include "CGAL/HalfedgeDS_default.h"
#include "CGAL/HalfedgeDS_min_items.h"
#include "CGAL/HalfedgeDS_decorator.h"
#include "CGAL/HalfedgeDS_items_2.h"
#include "CGAL/IO/Color.h"

#include "CGAL/Simple_cartesian.h"
#include"CGAL/Polyhedron_3.h"
#include "CGAL/draw_polyhedron.h"

using namespace std;

typedef CGAL::Simple_cartesian<double> Kernel;

int test_kernel() 
{
	Kernel::Point_2 p1(10, 10), p2(1, 1), q;
	Kernel::Line_2 l(p2,p1);
	cout << l.has_on_positive_side(Kernel::Point_2(9, 5)) << endl;
	cout << p1 << endl << p2 << endl;
	q = p1 + (p2 - p1) / 2.0;
	cout << (p2-p1).squared_length()<< endl;
	cout << "middle point: " << q << endl;

	cout << "affine transformation: " << endl;
	Kernel::Aff_transformation_2 trans;
	return 0;
}



/// ================================
template<class Refs>
struct MyFace : public CGAL::HalfedgeDS_face_base<Refs> {
	CGAL::IO::Color color;
	MyFace() {}
	MyFace(CGAL::IO::Color c) : color(c) {}
};

// An items type using my face.
struct My_items : public CGAL::HalfedgeDS_items_2 {
	template <class Refs, class Traits>
	struct Face_wrapper {
		typedef MyFace<Refs> Face;
	};
};
struct My_traits { // arbitrary point type, not used here.
	typedef int  Point_2;
}; 

struct PlaneEquation {
	template<class Facet>
	typename Facet::Plane_3 operator()(Facet& f) {
		typename Facet::Halfedge_handle h = f.halfedge();
		typedef typename Facet::Plane_3  Plane;
		return Plane(h->vertex()->point(),
			h->next()->vertex()->point(),
			h->next()->next()->vertex()->point());
	}
};


template <class Poly>
typename Poly::Halfedge_handle make_cube_3(Poly& P) {
	// appends a cube of size [0,1]^3 to the polyhedron P.
	CGAL_precondition(P.is_valid());
	typedef typename Poly::Point_3         Point;
	typedef typename Poly::Halfedge_handle Halfedge_handle;
	Halfedge_handle h = P.make_tetrahedron(Point(1, 0, 0),
		Point(0, 0, 1),
		Point(0, 0, 0),
		Point(0, 1, 0));
	Halfedge_handle g = h->next()->opposite()->next();             // Fig. (a)
	P.split_edge(h->next());
	P.split_edge(g->next());
	P.split_edge(g);                                              // Fig. (b)
	h->next()->vertex()->point() = Point(1, 0, 1);
	g->next()->vertex()->point() = Point(0, 1, 1);
	g->opposite()->vertex()->point() = Point(1, 1, 0);            // Fig. (c)
	Halfedge_handle f = P.split_facet(g->next(),
		g->next()->next()->next()); // Fig. (d)
	Halfedge_handle e = P.split_edge(f);
	e->vertex()->point() = Point(1, 1, 1);                        // Fig. (e)
	P.split_facet(e, f->next()->next());                          // Fig. (f)
	CGAL_postcondition(P.is_valid());
	return h;
}
int test7()
{
	typedef CGAL::Simple_cartesian<double>     Kernel;
	typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
	typedef Polyhedron::Halfedge_handle        Halfedge_handle;
	Polyhedron P;
	Halfedge_handle h = make_cube_3(P);
	CGAL::draw(P);
	return 0;
}

int test6()
{
	typedef CGAL::Simple_cartesian<double>               Kernel;
	typedef Kernel::Point_3                              Point_3;
	typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron;		//Kernel 组为Traits class
	typedef Polyhedron::Facet_iterator                   Facet_iterator;
	typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
	
	Point_3 p(0.0, 0.0, 0.0);
	Point_3 q(1.0, 0.0, 0.0);
	Point_3 r(0.0, 1.0, 0.0);
	Point_3 s(0.0, 0.0, 1.0);
	Polyhedron P;
	P.make_tetrahedron(p, q, r, s);
	// Write polyhedron in Object File Format (OFF).
	CGAL::IO::set_ascii_mode(std::cout);
	std::cout << "OFF" << std::endl << P.size_of_vertices() << ' '
		<< P.size_of_facets() << " 0" << std::endl;
	std::copy(P.points_begin(), P.points_end(),
		std::ostream_iterator<Point_3>(std::cout, "\n"));
	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		Halfedge_facet_circulator j = i->facet_begin();
		// Facets in polyhedral surfaces are at least triangles.
		std::cout << CGAL::circulator_size(j) << ' ';
		do {
			
			//顶点被顺序存储，计算当前遍历的顶点到起点的距离，
			std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());
		} while (++j != i->facet_begin());
		std::cout << std::endl;
	}
	return 0;
	
}


int test5()
{
	

	typedef CGAL::Simple_cartesian<double>     Kernel;
	typedef Kernel::Point_3                    Point_3;
	typedef Kernel::Plane_3                Plane_3;
	typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

	Point_3 p(1.0, 0.0, 0.0);
	Point_3 q(0.0, 1.0, 0.0);
	Point_3 r(0.0, 0.0, 1.0);
	Point_3 s(0.0, 0.0, 0.0);

	Polyhedron P;
	P.make_tetrahedron(p, q, r, s);

	std::transform(P.facets_begin(), P.facets_end(), P.planes_begin(), PlaneEquation());

	CGAL::IO::set_pretty_mode(std::cout);
	std::copy(P.planes_begin(), P.planes_end(),
		std::ostream_iterator<Plane_3>(std::cout, "\n"));
	return 0;
	


}

int test4()
{
	typedef CGAL::Simple_cartesian<double>     Kernel;
	typedef Kernel::Point_3                    Point_3;
	typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
	typedef Polyhedron::Vertex_iterator        Vertex_iterator;

	Point_3 p(1.0, 0.0, 0.0);
	Point_3 q(0.0, 1.0, 0.0);
	Point_3 r(0.0, 0.0, 1.0);
	Point_3 s(0.0, 0.0, 0.0);

	Polyhedron P;
	P.make_tetrahedron(p, q, r, s);
	CGAL::IO::set_ascii_mode(cout);

	/*for (Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
	{
		cout << v->point() << endl;
	}*/
	std::copy(P.points_begin(), P.points_end(),
		std::ostream_iterator<Point_3>(std::cout, "\n"));
	return 0;
}

int test3()
{
	//double 是fieldtype Simple_cartesian<FieldNumberType>
	typedef CGAL::Simple_cartesian<double>     Kernel;
	typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
	typedef Polyhedron::Halfedge_handle        Halfedge_handle;

	Kernel::FT a = 10.0;
	cout << a << endl;

	Polyhedron p;
	Halfedge_handle h= p.make_tetrahedron();
	cout << p.is_tetrahedron(h) << endl;
	return 0;
}

int test2()
{
	typedef CGAL::HalfedgeDS_default<My_traits, My_items> HDS;
	typedef HDS::Face	Face;
	typedef HDS::Face_handle FaceHandle;

	HDS hds;
	FaceHandle f = hds.faces_push_back(Face(CGAL::IO::red()));
	f->color = CGAL::IO::blue();

	assert(f->color == CGAL::IO::blue());
	return 0;
}

int test1()
{
	typedef CGAL::HalfedgeDS_default<My_traits> HDS;

	typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

	typedef CGAL::HalfedgeDS_default<int, CGAL::HalfedgeDS_min_items> HDS2;

	HDS hds;
	Decorator decorator(hds);
	decorator.create_loop();
	assert(decorator.is_valid());
	return 0;
}
