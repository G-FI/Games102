#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>

//检查重心坐标是否在convex之外
#include <vector>
#include<random>
#include<optional>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef Delaunay::Edge_iterator Edge_iterator;

typedef Delaunay::Point Point;
typedef K::Vector_2 Vector2;
typedef K::Segment_2 Segment;
typedef K::Triangle_2 Triangle;
typedef K::Ray_2 Ray;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef Delaunay::Vertex_handle Vertex_handle;

class Data {
public:
	int point_num = 20;
	std::vector<Point> filoyd_points;

	std::vector<Point> points;
	std::vector<ImVec2> rays_sampled;
	std::vector<ImVec2> seg_sampled;
	std::vector<ImVec2> convex_border;
	std::vector<ImVec2> triangles;

	bool show_voronoi = false;
	bool show_triangles = false;
	bool show_centrials = false;

	std::vector<Point> checker;

	Iso_rectangle_2 bbox;

	void Clear() {
		points.clear();
		filoyd_points.clear();
		rays_sampled.clear();
		seg_sampled.clear();
		convex_border.clear();
		triangles.clear();
	}
	void GenVoronoi(const Delaunay& dt) {
		std::vector<Segment> segments;
		std::vector<Ray> rays;


		for (Edge_iterator eit = dt.edges_begin(); eit != dt.edges_end(); ++eit) {
			CGAL::Object o = dt.dual(eit);
			if (CGAL::object_cast<K::Segment_2>(&o)) {
				//auto seg = CGAL::object_cast<Segment>(o);
				//auto so = CGAL::intersection(bbox, seg);
				//Segment s = CGAL::object_cast<Segment>(so);
				segments.push_back(CGAL::object_cast<Segment>(o));

			}
			else if (CGAL::object_cast<K::Ray_2>(&o)) {
				rays.push_back(CGAL::object_cast<Ray>(o));
			}
		}

		seg_sampled.clear();
		seg_sampled.resize(segments.size() * 2);
		for (int i = 0; i < (int)segments.size(); ++i) {
			seg_sampled[i * 2].x = segments[i].vertex(0).x();
			seg_sampled[i * 2].y = segments[i].vertex(0).y();
			seg_sampled[i * 2 + 1].x = segments[i].vertex(1).x();
			seg_sampled[i * 2 + 1].y = segments[i].vertex(1).y();
		}

		rays_sampled.clear();
		rays_sampled.resize(rays.size() * 2);
		for (int i = 0; i < (int)rays.size(); ++i) {
			rays_sampled[i * 2].x = rays[i].point(0).x();
			rays_sampled[i * 2].y = rays[i].point(0).y();
			rays_sampled[i * 2 + 1].x = rays[i].point(100).x();
			rays_sampled[i * 2 + 1].y = rays[i].point(100).y();
		}
	}

	void GenConvexBorder(const Delaunay& dt) {
		convex_border.clear();
		///循环遍历infinite vertex的邻接顶点就是convex border
		auto iv = dt.infinite_vertex();
		auto cv = iv->incident_vertices(), done(cv);
		do {
			auto x = cv->point().x();
			auto y = cv->point().y();
			convex_border.emplace_back(x, y);
			++cv;
		} while (cv != done);
		convex_border.emplace_back(cv->point().x(), cv->point().y());
	}

	void GenDelaunayTri(const Delaunay& dt) {
		//TODO 如何只是通过边，添加边上的顶点，能够减少顶点的绘制
		/*for (auto fei = dt.finite_edges_begin(); fei != dt.finite_edges_end(); ++fei) {
			fei->first->vertex
		}*/
		triangles.clear();

		for (auto ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end(); ++ffi) {
			triangles.emplace_back(ffi->vertex(0)->point().x(), ffi->vertex(0)->point().y());
			triangles.emplace_back(ffi->vertex(1)->point().x(), ffi->vertex(1)->point().y());
			triangles.emplace_back(ffi->vertex(2)->point().x(), ffi->vertex(2)->point().y());
		}
	}
	void FIloyd(Delaunay& dt) {
		//找到被voronoi区域包围的顶点
		std::vector<Vertex_handle> old_vertices;
		std::vector<Point> new_points;
		//找到封闭的面的边界点
		//使用边界点计算重心坐标
		//一种中间顶点到中心坐标

		//finite face这些概念都表示的Delanuay 三角形区域的

		for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
			//当前正在遍历的顶点坐标
			auto p = vit->point();


			auto fi = vit->incident_faces(), done(fi);
			bool is_border = false;
			do {
				if (dt.is_infinite(fi)) {
					//顶点相邻infinite face，不更新
					new_points.push_back(p);
					is_border = true;
					break;
				}
				++fi;
			} while (fi != done);


			if (is_border) {
				continue;
			}
			Point centroid = GetCentroid(dt, vit);

			new_points.push_back(centroid);
			old_vertices.emplace_back(vit->handle());
		}

		filoyd_points = new_points;

	}
	void Move() {
		points = filoyd_points;
	}
private:
	/// <summary>
	/// 将vit邻接Delauney三角形的外接圆心作为voronoi边界的顶点，计算这些顶点的重心坐标
	/// 有些vit的邻接三角形为钝角三角形，钝角的外界圆心在三角形之外，因此voronoi边界顶点在bbox之外，因此中心坐标也容易在bbox之外
	/// 
	///Bug: 当vit有两个相邻的voronoi顶点都位于bbox之外时，这两个顶点的线段可能与bbox有0，2个交点
	/// 当0个顶点时OK
	/// 当2个交点时，相交算法默认返回线段与bbox的上边，左边的顶点(实际应该返回靠近当前正在操作的这个顶点的交点)
	/// 因此会导致重心计算错误
	/// </summary>
	/// <param name="dt"></param>
	/// <param name="vit">vit周围都是finite face</param>
	/// <returns></returns>
	Point GetCentroid(const Delaunay& dt, Vertex_handle vit) {
		auto p = vit->point();
		//找到顶点对应的Voronoi顶点
		std::vector<Point> v_neighbours;

		auto fi = vit->incident_faces(), done(fi);
		do {

			auto voronoi_p = dt.dual(fi);
			//判断voronoi_p 是否超出边界
			Segment s0(p, voronoi_p);
			/*{
				CGAL::Object o = CGAL::intersection(s0, bbox);
				const Point* inter_p = CGAL::object_cast<Point>(&o);
				if (inter_p) {
					//相交，添加两个顶点
					auto p1 = dt.dual(fi->neighbor(fi->ccw(0)));
					auto p2 = dt.dual(fi->neighbor(fi->cw(0)));
					Segment s1(p1, voronoi_p);
					Segment	s2(p2, voronoi_p);
					CGAL::Object o1 = CGAL::intersection(bbox, s1);
					CGAL::Object o2 = CGAL::intersection(bbox, s1);
					const Point* new_p1 = CGAL::object_cast<Point>(&o1);
					const Point* new_p2 = CGAL::object_cast<Point>(&o2);
					v_neighbours.push_back(*new_p1);
					v_neighbours.push_back(*new_p2);

				}
				else {
					//不相交
					v_neighbours.push_back(voronoi_p);
				}
			}
			*/
			auto opt_p = Intersect_RL(bbox, s0);
			if (opt_p.has_value()) {
				//相交
				//不应该时face的neibors，这个顶点的邻接面，而这两个邻接面又与这个钝角三角形相邻
				std::vector<Point> ps;
				for (int i = 0; i < 3; ++i) {
					//相邻的面不能是infinite
					if (dt.is_infinite(fi->neighbor(fi->ccw(i)))) {
						continue;
					}
					ps.push_back(dt.dual(fi->neighbor(fi->ccw(i))));
				}
				Segment s1(ps[0], voronoi_p);
				Segment	s2(ps[1], voronoi_p);
				//问题，超出边界的voronoi点，的相邻的voronoi点也超出边界的话就会与bbox产生两个交点，甚至没有交点
				auto opt_p1 = Intersect_RL(bbox, s1);
				auto opt_p2 = Intersect_RL(bbox, s2);
				if (opt_p1.has_value()) {
					v_neighbours.push_back(opt_p1.value());
				}
				if (opt_p2.has_value()) {
					v_neighbours.push_back(opt_p2.value());

				}
			}
			else {
				//不相交
				v_neighbours.push_back(voronoi_p);
			}
			
			//v_neighbours.push_back(voronoi_p);
			++fi;
		} while (fi != done);

		double avg_x = 0.0, avg_y = 0.0;
		for (auto v : v_neighbours) {
			avg_x += v.x();
			avg_y += v.y();
		}
		avg_x /= v_neighbours.size();
		avg_y /= v_neighbours.size();
		return Point(avg_x, avg_y);
	}
	double Cross2d(Vector2 p1, Vector2 p2) {
		return p1.x() * p2.y() - p1.y() * p2.x();
	}
	bool SegmentIntersect(const Segment& s1, const Segment& s2) {
		//使用跨线段计算
		Point a = s1.source();
		Point b = s1.target();
		Point c = s2.source();
		Point d = s2.target();
		if (Cross2d(c - a, b - a) * Cross2d(b - a, d - a) < 0)
			return false;
		if (Cross2d(a - c, d - c) * Cross2d(d - c, b - c) < 0)
			return false;
		return true;
	}
	std::optional<Point> Intersect_RL(Iso_rectangle_2& r, const Segment& s) {
		Segment up(Point(r.xmin(), r.ymax()), Point(r.xmax(), r.ymax()));
		Segment bottom(Point(r.xmin(), r.ymin()), Point(r.xmax(), r.ymin()));
		Segment left(Point(r.xmin(), r.ymin()), Point(r.xmin(), r.ymax()));
		Segment right(Point(r.xmax(), r.ymin()), Point(r.xmax(), r.ymax()));
		double x = 0.0, y = 0.0, dx = 0.0, dy = 0.0;
		double k = (s.target().y() - s.source().y()) / (s.target().x() - s.source().x());
		if (SegmentIntersect(up, s)) {
			y = r.ymax();
			dy = y - s.source().y();
			dx = dy / k;
			x = s.source().x() + dx;
			return Point(x, y);
		}
		else if (SegmentIntersect(bottom, s)) {
			y = r.ymin();
			dy = y - s.source().y();
			dx = dy / k;
			x = s.source().x() + dx;
			return Point(x, y);

		}
		else if (SegmentIntersect(left, s)) {
			x = r.xmin();
			dx = x - s.source().x();
			dy = dx * k;
			y = s.source().y() + dy;
			return Point(x, y);
		}
		else if (SegmentIntersect(right, s)) {
			x = r.xmax();
			dx = x - s.source().x();
			dy = dx * k;
			y = s.source().y() + dy;
			return Point(x, y);
		}
		return {};
	}
};

GLFWwindow* createWindow();
void NextFrame(GLFWwindow* window);
void EndFrame(GLFWwindow* window);
void CleanContext(GLFWwindow* window);
ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) {
	return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y);
}

ImVec2 operator+(const ImVec2& lhs, const ImVec2& rhs) {
	return ImVec2(lhs.x + rhs.x, lhs.y + rhs.y);
}

ImVec2 operator+(const Point& lhs, const ImVec2& rhs) {
	return ImVec2(lhs.x() + rhs.x, lhs.y() + rhs.y);
}
//The verticesand faces of the triangulations are accessed through handles, iteratorsand circulators.

//void test(Data &data, Point p){
//	data.checker.clear();
//
//	Delaunay dt(data.points.begin(), data.points.end());
//	std::vector<Point> tmp_convex_border;
//	///循环遍历infinite vertex的邻接顶点就是convex border
//	auto iv = dt.infinite_vertex();
//	auto cv = iv->incident_vertices(), done(cv);
//	do {
//		auto x = cv->point().x();
//		auto y = cv->point().y();
//		tmp_convex_border.emplace_back(x, y);
//		++cv;
//	} while (cv != done);
//	
//	tmp_convex_border.push_back(p);
//	CGAL::ch_graham_andrew(tmp_convex_border.begin(), tmp_convex_border.end(), std::back_inserter(data.checker));
//	for (auto p : data.checker) {
//		std::cout << p << std::endl;
//	}
//	std::cout<<"========================"<<std::endl;
//}

void Voronoi(Data& data) {
	Delaunay dt(data.points.begin(), data.points.end());

	data.GenDelaunayTri(dt);
	data.GenConvexBorder(dt);
	data.GenVoronoi(dt);
}
void FIloyd(Data& data) {
	Delaunay dt(data.points.begin(), data.points.end());
	data.FIloyd(dt);
}

void GenPoints(Data& data) {
	std::random_device r;
	std::default_random_engine e(r());
	std::uniform_real_distribution<float> uniform_dist(0.f, 700.f);

	data.Clear();
	data.points.resize(data.point_num);
	for (int i = 0; i < data.point_num; ++i) {
		data.points[i] = Point(uniform_dist(e), uniform_dist(e));
	}


}

void Monitor(Data& data, ImVec2 mouse_in_canvas) {
	ImGuiIO& io = ImGui::GetIO();
	// ImGuiWindowFlags window_flags = ImGuiWindowFlags_AlwaysAutoResize |
	// ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
	// ImGuiWindowFlags_NoNav; const float PAD = 10.f; ImVec2 work_pos =
	// ImGui::GetMainViewport()->WorkPos; work_pos.x += PAD; work_pos.y += PAD;

	// ImVec2 work_pos_pivot(0.f, 0.f);
	// ImGui::SetNextWindowPos(work_pos, ImGuiCond_Always, work_pos_pivot);
	// window_flags |= ImGuiWindowFlags_NoMove;

	if (ImGui::Begin("moniter")) {
		ImGui::InputInt("number of point", &data.point_num, 1);
		if (ImGui::Button("generate point")) {
			GenPoints(data);
		}
		if (ImGui::Button("calculate")) {
			Voronoi(data);
			FIloyd(data);
		}
		ImGui::Checkbox("Voronoi", &data.show_voronoi);
		ImGui::Checkbox("Triangles", &data.show_triangles);
		ImGui::Checkbox("Centrials", &data.show_centrials);

		if (ImGui::Button("move")) {
			data.Move();
		}


		ImGui::SameLine();
		if (ImGui::Button("clear")) {
			data.Clear();
		}


		ImGui::Text("Mouse pos(%.2f, %.2f)", io.MousePos.x, io.MousePos.y);
		ImGui::Text("pos in canvas(%.2f, %2f)", mouse_in_canvas.x, mouse_in_canvas.y);
	}
	ImGui::End();
}

/// <summary>
/// TODO
/// CGAL三角化
/// 封装Imgui常用功能
/// 使用cmake生成链接库，加速编译
/// </summary>
/// <returns></returns>
int main() {
	auto window = createWindow();

	Iso_rectangle_2 bbox(-1, -1, 2, 2);
	Data data;
	data.bbox = Iso_rectangle_2(0, 0, 700, 700);

	while (!glfwWindowShouldClose(window)) {
		NextFrame(window);
		if (ImGui::Begin("Canvas")) {
			// Using InvisibleButton() as a convenience 1) it will advance the layout
			// cursor and 2) allows us to use IsItemHovered()/IsItemActive()
			ImVec2 origin = ImGui::GetCursorScreenPos(); // ImDrawList API uses screen
			// coordinates!
			ImVec2 canvas_sz(700.f, 700.f);
			//ImVec2 canvas_sz =
			//	ImGui::GetContentRegionAvail(); // Resize canvas to what's available
			//if (canvas_sz.x < 50.0f)
			//	canvas_sz.x = 50.0f;
			//if (canvas_sz.y < 50.0f)
			//	canvas_sz.y = 50.0f;

			// 绘制Canvas
			ImGuiIO& io = ImGui::GetIO();
			ImDrawList* draw_list = ImGui::GetWindowDrawList();
			draw_list->AddRectFilled(origin, origin + canvas_sz, IM_COL32(50, 50, 50, 255));
			draw_list->AddRect(origin, origin + canvas_sz, IM_COL32(255, 255, 255, 255));


			// 获取是否点击（使用canvas_sz是，input之下的那些可用区域）
			ImGui::InvisibleButton("canvas", canvas_sz,
				ImGuiButtonFlags_MouseButtonLeft);
			ImVec2 mouse_in_canvas = io.MousePos - origin;
			Monitor(data, mouse_in_canvas);



			if (ImGui::IsItemActive() && ImGui::IsItemClicked()) {
				data.points.emplace_back(mouse_in_canvas.x, mouse_in_canvas.y);
			}

			for (const auto& point : data.points) {
				draw_list->AddCircleFilled(point + origin, 4.f,
					IM_COL32(255, 0, 0, 255));
			}

			if (data.show_centrials) {
				for (const auto& point : data.filoyd_points) {
					draw_list->AddCircleFilled(point + origin, 2.f,
						IM_COL32(255, 255, 255, 255));
				}
			}

			//Voronoi 区域线段
			if (data.show_voronoi) {
				for (int i = 0; i < (int)data.seg_sampled.size(); i += 2) {
					draw_list->AddLine(origin + data.seg_sampled[i], origin + data.seg_sampled[i + 1], IM_COL32(128, 128, 0, 255));

				}
				//Voronoi 区域生成的射线
				for (int i = 0; i < (int)data.rays_sampled.size(); i += 2) {
					draw_list->AddLine(origin + data.rays_sampled[i], origin + data.rays_sampled[i + 1], IM_COL32(255, 0, 0, 255));
				}
			}
			if (data.show_triangles) {

				//Delaunay 三角形
				for (int i = 0; i < (int)data.triangles.size(); i += 3) {
					draw_list->AddLine(origin + data.triangles[i], origin + data.triangles[i + 1], IM_COL32(255, 255, 0, 255));
					draw_list->AddLine(origin + data.triangles[i + 1], origin + data.triangles[i + 2], IM_COL32(255, 255, 0, 255));
					draw_list->AddLine(origin + data.triangles[i + 2], origin + data.triangles[i], IM_COL32(255, 255, 0, 255));
				}
				//凸包边界
				for (int i = 1; i < (int)data.convex_border.size(); i += 1) {

					draw_list->AddLine(origin + data.convex_border[i - 1], origin + data.convex_border[i], IM_COL32(0, 0, 255, 255));
				}
			}
		}
		ImGui::End();

		EndFrame(window);
	}

	CleanContext(window);
}

void NextFrame(GLFWwindow* window) {
	// 轮询事件
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

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

GLFWwindow* createWindow() {
	// glfw: initialize and configure
	// ------------------------------
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// msaa
	glfwWindowHint(GLFW_SAMPLES, 4);
	// glfw window creation
	// --------------------
	GLFWwindow* window = glfwCreateWindow(1280, 800, "LearnOpenGL", NULL, NULL);
	if (window == NULL) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return NULL;
	}
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync

	// 鼠标不可见，内嵌进应用中了
	// glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	// ȫ�ֻص�����
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	// glfwSetCursorPosCallback(window, mouse_callback);
	// glfwSetScrollCallback(window, scroll_callback);

	// glad: load all OpenGL function pointers
	//---------------------------------------
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cout << "Failed to initialize GLAD" << std::endl;
		return NULL;
	}
	{

		const char* glsl_version = "#version 130";
		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		ImGuiIO& io = ImGui::GetIO();
		(void)io;

		io.AddKeyEvent(ImGuiKey_Escape, true);

		// Setup Dear ImGui style
		ImGui::StyleColorsDark();

		// Setup Platform/Renderer backends
		ImGui_ImplGlfw_InitForOpenGL(window, true);
		ImGui_ImplOpenGL3_Init(glsl_version);
	}
	return window;
}
