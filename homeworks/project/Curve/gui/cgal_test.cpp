//#include <iostream>
//#include<vector>
//
//#include "CGAL/Simple_cartesian.h"
//#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
//#include"CGAL/Exact_predicates_inexact_constructions_kernel.h"
//#include "CGAL/convex_hull_2.h"
//
//using namespace std;
//
//int PointAndSegment()
//{
//    typedef CGAL::Simple_cartesian<double> Kernel;
//    typedef Kernel::Point_2 Point_2;
//    typedef Kernel::Segment_2 Segment_2;
//    Point_2 p(1, 1), q(10, 10);
//    std::cout << "p = " << p << std::endl;
//    std::cout << "q = " << q.x() << " " << q.y() << std::endl;
//    std::cout << "sqdist(p,q) = "
//        << CGAL::squared_distance(p, q) << std::endl;
//    Segment_2 s(p, q);
//    Point_2 m(5, 9);
//    std::cout << "m = " << m << std::endl;
//    std::cout << "sqdist(Segment_2(p,q), m) = "
//        << CGAL::squared_distance(s, m) << std::endl;
//    std::cout << "p, q, and m ";
//    switch (CGAL::orientation(p, q, m)) {
//    case CGAL::COLLINEAR:
//        std::cout << "are collinear\n";
//        break;
//    case CGAL::LEFT_TURN:
//        std::cout << "make a left turn\n";
//        break;
//    case CGAL::RIGHT_TURN:
//        std::cout << "make a right turn\n";
//        break;
//    }
//    std::cout << " midpoint(p,q) = " << CGAL::midpoint(p, q) << std::endl;
//    return 0;
//}
//void Pricies()
//{
//    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//    typedef Kernel::Point_2 Point_2;
//    typedef Kernel::Segment_2 Segment_2;
//    {
//        Point_2 p(0, 0.3), q(1, 0.6), r(2, 0.9);
//        std::cout << (CGAL::collinear(p, q, r) ? "collinear" : "non_colleanear") << std::endl;
//    }
//    {
//        Point_2 p(0, 0.3), q(1, 0.6), r(2, 0.9);
//
//    }
//}
//void ConvexHull1()
//{
//    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//    typedef Kernel::Point_2 Point_2;
//    Point_2 points[5] = { Point_2(0,0), Point_2(10,0), Point_2(10,10), Point_2(6,5), Point_2(4,1) };
//    Point_2 result[5];
//    Point_2* ptr = CGAL::convex_hull_2(points, points + 5, result);
//    cout << "number of points in convex hull: " << ptr - result << endl;
//    for (int i = 0; i < ptr - result; ++i)
//    {
//        cout << result[i] << endl;
//    }
//}
//void ConvexHull2()
//{
//    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//    typedef Kernel::Point_2 Point_2;
//    vector<Point_2> points;
//    vector<Point_2> result;
//    points.emplace_back(0, 0);
//    points.emplace_back(10, 0);
//    points.emplace_back(10, 10);
//    points.emplace_back(6, 5);
//    points.emplace_back(4, 1);
//
//    CGAL::convex_hull_2(points.begin(), points.end(), back_inserter(result));
//
//    cout << "number of points in convex hull: " << result.size() << endl;
//    for (auto& point : result)
//    {
//        cout << point << endl;
//    }
//}
//
//
//
//// standard includes
//#include <iostream>
//#include <fstream>
//#include <cassert>
//// example that uses the filtered traits and
//// the segment Delaunay graph hierarchy
//// choose the kernel
//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;
//typedef K::Point_2                     Point_2;
//
//
//// typedefs for the traits and the algorithm
//#include <vector>
//#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
//#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
//typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<K> Gt;
//typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt>  SDG2;
//int main()
//{
//
//    Point_2 /*p1(0, 0), */p1(1, 1), p2(2, 0);
//    std::vector<SDG2::Site_2> sites;
//    sites.emplace_back(SDG2::Site_2::construct_site_2(p1));
//    sites.emplace_back(SDG2::Site_2::construct_site_2(p2));
//
//    SDG2          sdg;
//    SDG2::Site_2  site =
//        SDG2::Site_2::construct_site_2(K::Point_2(CGAL::ORIGIN));
//
//    sites.push_back(site);
//
//    for (auto s : sites) {
//        sdg.insert(s);
//    }
//    // validate the segment Delaunay graph
//    assert(sdg.is_valid(true, 1));
//
//    cout << endl << endl;
//    // print the number of input and output sites
//    cout << "# of input sites : " << sdg.number_of_input_sites() << endl;
//    cout << "# of output sites: " << sdg.number_of_output_sites() << endl;
//    unsigned int n_ipt(0), n_iseg(0), n_opt(0), n_oseg(0), n_ptx(0);
//    // count the number of input points and input segments
//    SDG2::Input_sites_iterator iit;
//    for (iit = sdg.input_sites_begin(); iit != sdg.input_sites_end(); ++iit)
//    {
//        if (iit->is_point()) { n_ipt++; }
//        else { n_iseg++; }
//    }
//    // count the number of output points and output segments, as well
//    // as the number of points that are points of intersection of pairs
//    // of strongly intersecting sites
//    SDG2::Output_sites_iterator oit;
//    for (oit = sdg.output_sites_begin(); oit != sdg.output_sites_end(); ++oit)
//    {
//        if (oit->is_segment()) { n_oseg++; }
//        else {
//            n_opt++;
//            if (!oit->is_input()) { n_ptx++; }
//        }
//    }
//    cout << endl << "# of input segments:  " << n_iseg << endl;
//    cout << "# of input points:    " << n_ipt << endl << endl;
//    cout << "# of output segments: " << n_oseg << endl;
//    cout << "# of output points:   " << n_opt << endl << endl;
//    cout << "# of intersection points: " << n_ptx << endl;
//    return 0;
//    return 0;
//}