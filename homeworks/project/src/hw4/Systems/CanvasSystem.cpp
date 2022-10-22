#include "CanvasSystem.h"

#include "../Components/CanvasData.h"

#include <_deps/imgui/imgui.h>

using namespace Ubpa;
using namespace std;
#include<Eigen/Core>
#include<Eigen/Dense>
using namespace Eigen;

#include<vector>
#include<sstream>
#include <spdlog/spdlog.h>


std::vector<float> uniform_parameterization(const std::vector<Ubpa::pointf2>& points) {
	int N = points.size();
	std::vector<float> t(N);
	if (N > 1) {
		float start_t = 0.f;
		for (int i = 0; i < N; ++i) {
			t[i] = (float)i / (N - 1);
		}
	}
	return t;
}

std::vector<float> chordal_parameterization(const std::vector<Ubpa::pointf2>& points) {
	int N = points.size();
	std::vector<float> t(N);
	if (N > 1) {
		auto length = [](const auto& p1, const auto& p2)->float {
			return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]));
		};

		//限定采样点 t在[0, 1]之间
		float start_t = 0.f;

		t[0] = start_t;
		float total_chrod = 0.f;
		std::vector<float> chord_length(points.size() - 1);

		for (int i = 0; i < points.size() - 1; ++i) {
			chord_length[i] = length(points[i + 1], points[i]);
			total_chrod += chord_length[i];
		}

		for (int i = 0; i < chord_length.size(); ++i) {
			t[i + 1] = chord_length[i] / total_chrod + t[i];
		}
	}
	return t;
}

MatrixXf spline(const vector<float>& x, const vector<float>& y) {
	int N = x.size();
	VectorXf dx = VectorXf::Zero(N - 1);
	VectorXf dy = VectorXf::Zero(N - 1);
	MatrixXf coeff = MatrixXf::Zero(N - 1, 3);	//[b, c, d]

	for (int i = 0; i < x.size() -1 ; ++i) {
		dx(i) = x[i + 1] - x[i];
		dy(i) = y[i + 1] - y[i];
	}
	
	//计算ci
	/*for (int i = 1; i < N; ++i) {
		
	}*/
	//TODO 追赶法
	MatrixXf A = MatrixXf::Zero(N, N);
	VectorXf r = VectorXf::Zero(N);
	A(0, 0) = 1;
	A(N - 1, N - 1) = 1;
	for (int i = 1; i < N - 1; ++i) {
		A(i, i - 1) = dx(i-1);
		A(i, i) = 2 * (dx(i - 1) + dx(i));
		A(i, i + 1) = dx(i);
		r(i) = 3 * (dy(i) / dx(i) - dy(i - 1) / dx(i - 1));
	}

	VectorXf c = A.colPivHouseholderQr().solve(r);

	coeff.col(1) = c.head(N - 1);
	for (int i = 0; i < N - 1; ++i) {
		coeff(i, 2) = (c(i+1) - c(i)) / (3 * dx(i));
		coeff(i, 0) = dy(i) / dx(i) - dx(i) * (2 * c(i) + c(i + 1)) / 3;
	}
	return coeff;
}
void sampling(CanvasData* data, int k, const MatrixXf& coeff) {
	//每两个控制点之间采样
	data->sample_points.clear();
	for (int i = 0; i < data->points.size() - 1; ++i) {
		float len = data->points[i + 1][0] - data->points[i][0];
		for (int j = 0; j < k; ++j) {
			float dx = (float)j / k * len;
			float x = data->points[i][0] + dx;
			float y = coeff(i, 2) * dx;
			y = (y + coeff(i, 1)) * dx;
			y = (y + coeff(i, 0)) * dx + data->points[i][1];
			data->sample_points.push_back(Ubpa::pointf2(x, y));
		}
	}
}
void sampling(CanvasData* data, const vector<float>& ts, int k, const MatrixXf& coeffX, const MatrixXf& coeffY, bool uniform) {
	if (uniform) {
		data->sample_points1.clear();
	}
	else {
		data->sample_points2.clear();
	}
	for (int i = 0; i < data->points.size() - 1; ++i) {
		float len = ts[i + 1] - ts[i];

		for (int j = 0; j < k; ++j) {
			float dt = (float)j / k * len;
			float t = ts[i] + dt;

			float x = coeffX(i, 2) * dt;
			x = (x + coeffX(i, 1)) * dt;
			x = (x + coeffX(i, 0)) * dt + data->points[i][0];

			float y = coeffY(i, 2) * dt;
			y = (y + coeffY(i, 1)) * dt;
			y = (y + coeffY(i, 0)) * dt + data->points[i][1];
			if (uniform) {
				data->sample_points1.push_back(Ubpa::pointf2(x, y));
			}
			else {
				data->sample_points2.push_back(Ubpa::pointf2(x, y));
			}
		}
	}
}
void GenControlPoint(CanvasData* data) {
	//uniform参数化
	//n个性质点，对应n-1段曲线，每段曲线两端各绘制一个UI线段
	float len = 0.1f;

	//存储导数
	data->dx_p.resize(data->points.size() - 1);
	data->dy_p.resize(data->points.size() - 1);
	data->dx_n.resize(data->points.size() - 1);
	data->dy_n.resize(data->points.size() - 1);

	for (int i = 0; i < data->points.size() - 1; ++i) {
		//导数方向，目前只修改均匀参数化曲线
		float dx1 = data->coeffXu(i, 0);  //bi时dx/dt；
		float dy1 = data->coeffYu(i, 0);
		data->dx_p[i] = dx1;
		data->dy_p[i] = dy1;



		float dt = data->tu[i + 1] - data->tu[i];
		float dx2 = data->coeffXu(i, 0) + 2 * data->coeffXu(i, 1) * dt + 3 * data->coeffXu(i, 2) * dt * dt;
		float dy2 = data->coeffYu(i, 0) + 2 * data->coeffYu(i, 1) * dt + 3 * data->coeffYu(i, 2) * dt * dt;
		data->dx_n[i] = dx2;
		data->dy_n[i] = dy2;
	}
}

void CheckPoint(CanvasData* data, const ImVec2& origin) {
	ImGuiIO& io = ImGui::GetIO();
	//ImVec2 click_pos = io.MouseClickedPos[0];
	float len = 0.1f;
	for (int i = 0; i < data->points.size() - 1; ++i) {
		ImVec2 min(origin.x + data->points[i][0] - 6.f, origin.y + data->points[i][1] - 6.f);
		ImVec2 max(origin.x + data->points[i][0] + 6.f, origin.y + data->points[i][1] + 6.f);
		if (ImGui::IsMouseHoveringRect(min, max)) {
			data->is_darging = true;
			data->draging_idx = i;
			break;
		}

		float dx1 = data->dx_p[i];
		float dy1 = data->dy_p[i];
		min.x = origin.x + data->points[i][0] + dx1 * len - 5.f;
		min.y = origin.y + data->points[i][1] + dy1 * len - 5.f;
		max.x = min.x + 10.f;
		max.y = min.y + 10.f;
		if (ImGui::IsMouseHoveringRect(min, max)){
			data->edit_idx = i;
			data->changing_point = true;
			data->positive = true;
			break;
		}

		float dx2 = data->dx_n[i]; 
		float dy2 = data->dy_n[i];
		min.x = origin.x + data->points[i+1][0] - dx2 * len - 5.f;
		min.y = origin.y + data->points[i+1][1] - dy2 * len - 5.f;
		max.x = min.x + 10.f;
		max.y = min.y + 10.f;
		if (ImGui::IsMouseHoveringRect(min, max)) {
			data->edit_idx = i + 1;
			data->changing_point = true;
			data->positive = false;
			break;
		}
	}
	if (!data->is_darging && !data->changing_point) {//检查最后一个点，因为它没有在上面的循环中
		int i = data->points.size() - 1;
		ImVec2 min(origin.x + data->points[i][0] - 6.f, origin.y + data->points[i][1] - 6.f);
		ImVec2 max(origin.x + data->points[i][0] + 6.f, origin.y + data->points[i][1] + 6.f);
		if (ImGui::IsMouseHoveringRect(min, max)) {
			data->is_darging = true;
			data->draging_idx = i;
		}
	}
}
void SmoothVertex(CanvasData* data, float dx_new, float dy_new) {
	//调整第一个值型点或最后一个，只需要求解一段曲线，两组方程
	if (data->edit_idx == 0 || data->edit_idx == data->tu.size() - 1) {
		int i = 0;
		Matrix3f A;
		Vector3f coeff;
		Vector3f B;
		if (data->edit_idx == 0) {//更新第0段曲线
			i = 0;
			data->dx_p[i] = dx_new;
			data->dy_p[i] = dy_new;
		}
		else {//更新第n-1段曲线
			i = data->edit_idx - 1;
			data->dx_n[i] = dx_new;
			data->dy_n[i] = dy_new;
		}
		float dt = data->tu[i + 1] - data->tu[i];
		float sxi = data->points[i][0];
		float syi = data->points[i][1];

		float sxi1 = data->points[i + 1][0];
		float syi1 = data->points[i + 1][1];

		float dx = data->dx_p[i];
		float dy = data->dy_p[i];
		float dx1 = data->dx_n[i];
		float dy1 = data->dy_n[i];

		A << dt, dt* dt, dt* dt* dt,
			1, 0, 0,
			1, 2 * dt, 3 * dt * dt;
		B << sxi1 - sxi, dx, dx1;
		auto tmp = A.colPivHouseholderQr();
		data->coeffXu.row(i) = tmp.solve(B).transpose();
		B << syi1 - syi, dy, dy1;
		data->coeffYu.row(i) = tmp.solve(B).transpose();
	}
	else {//更新第i-1段曲线和第i段曲线，更新4组参数 0 < i < n-1;
		data->dx_n[data->edit_idx - 1] = data->dx_p[data->edit_idx] = dx_new;
		data->dy_n[data->edit_idx - 1] = data->dy_p[data->edit_idx] = dy_new;
		//先假定只更新中间的值型点
		Matrix3f A;
		Vector3f coeff;
		Vector3f B;
		float sxi_1 = data->points[data->edit_idx - 1][0];
		float syi_1 = data->points[data->edit_idx - 1][1];
		float sxi = data->points[data->edit_idx][0];
		float syi = data->points[data->edit_idx][1];
		float sxi1 = data->points[data->edit_idx + 1][0];
		float syi1 = data->points[data->edit_idx + 1][1];

		float dx1 = data->dx_p[data->edit_idx - 1];
		float dy1 = data->dy_p[data->edit_idx - 1];
		float dx2 = data->dx_n[data->edit_idx];
		float dy2 = data->dy_n[data->edit_idx];

		float dt1 = data->tu[data->edit_idx] - data->tu[data->edit_idx - 1];
		float dt2 = data->tu[data->edit_idx + 1] - data->tu[data->edit_idx];

		//更新si-1段参数
		A << dt1, dt1* dt1, dt1* dt1* dt1,
			1, 0, 0,
			1, 2 * dt1, 3 * dt1 * dt1;
		B << sxi - sxi_1, dx1, dx_new;
		auto tmp = A.colPivHouseholderQr();
		coeff = tmp.solve(B);
		data->coeffXu.row(data->edit_idx - 1) = coeff.transpose();

		B << syi - syi_1, dy1, dy_new;
		coeff = tmp.solve(B);
		data->coeffYu.row(data->edit_idx - 1) = coeff.transpose();

		//更新si段参数
		A << dt2, dt2* dt2, dt2* dt2* dt2,
			1, 0, 0,
			1, 2 * dt2, 3 * dt2 * dt2;
		B << sxi1 - sxi, dx_new, dx2;

		tmp = A.colPivHouseholderQr();
		coeff = tmp.solve(B);
		data->coeffXu.row(data->edit_idx) = coeff.transpose();

		B << syi1 - syi, dy_new, dy2;
		coeff = tmp.solve(B);
		data->coeffYu.row(data->edit_idx) = coeff.transpose();
	}

}
void LineVertex(CanvasData* data, float dx_new, float dy_new) {
	if (data->edit_idx == 0 || data->edit_idx == data->tu.size() - 1) {//只更新一段曲线
		//更新端点导数
		int i = -1;
		Matrix3f A;
		Vector3f coeff;
		Vector3f B;
		if (data->edit_idx == 0) {//更新第0段曲线
			i = 0;
			data->dx_p[i] = dx_new;
			data->dy_p[i] = dy_new;
		}
		else {//更新第n-1段曲线
			i = data->edit_idx - 1;
			data->dx_n[i] = dx_new;
			data->dy_n[i] = dy_new;
		}
		float dt = data->tu[i + 1] - data->tu[i];
		float sxi = data->points[i][0];
		float syi = data->points[i][1];

		float sxi1 = data->points[i + 1][0];
		float syi1 = data->points[i + 1][1];

		float dx = data->dx_p[i];
		float dy = data->dy_p[i];
		float dx1 = data->dx_n[i];
		float dy1 = data->dy_n[i];

		A << dt, dt* dt, dt* dt* dt,
			1, 0, 0,
			1, 2 * dt, 3 * dt * dt;
		B << sxi1 - sxi, dx, dx1;
		auto tmp = A.colPivHouseholderQr();
		data->coeffXu.row(i) = tmp.solve(B).transpose();
		B << syi1 - syi, dy, dy1;
		data->coeffYu.row(i) = tmp.solve(B).transpose();
		//更新更新系数
	}
	else {//更新Si-1和Si两端曲线
		//更新端点导数
		int i = data->edit_idx;
		float length_new = sqrt(dx_new*dx_new + dy_new*dy_new);
		float cos = dx_new / length_new;
		float sin = dy_new / length_new;
		if (data->positive) {
			data->dx_p[i] = dx_new;
			data->dy_p[i] = dy_new;
			float length_old = sqrt(data->dx_n[i - 1] * data->dx_n[i - 1] + data->dy_n[i - 1] * data->dy_n[i - 1]);
			data->dx_n[i - 1] = length_old * cos;
			data->dy_n[i - 1] = length_old * sin;
		}
		else {
			data->dx_n[i-1] = dx_new;
			data->dy_n[i-1] = dy_new;
			float length_old = sqrt(data->dx_p[i] * data->dx_p[i ] + data->dy_p[i] * data->dy_p[i]);
			data->dx_p[i] = length_old * cos;
			data->dy_p[i] = length_old * sin;
		}
		
		float dt1 = data->tu[i] - data->tu[i - 1];
		float dt2 = data->tu[i + 1] - data->tu[i];


		//sx表示t关于x的参数，i_1表示t(i-1)， i， i1表示t(i+1)点
		float sxi_1 = data ->points[i - 1][0];
		float sxi = data->points[i][0];
		float sxi1 = data->points[i + 1][0];

		float syi_1 = data->points[i - 1][1];
		float syi = data->points[i][1];
		float syi1 = data->points[i + 1][1];
		

		
		//dx1,2,3,4顺序从左到右分别为dx_p[i-1], dx_n[i-1], dx_p[i], dx_n[i]
		float dx1 = data->dx_p[i - 1];
		float dx2 = data->dx_n[i - 1];
		float dx3 = data->dx_p[i];
		float dx4 = data->dx_n[i];

		float dy1 = data->dy_p[i - 1];
		float dy2 = data->dy_n[i - 1];
		float dy3 = data->dy_p[i];
		float dy4 = data->dy_n[i];

		Matrix3f A;
		Vector3f B;
		Vector3f coeff;

		A << dt1, dt1* dt1, dt1* dt1* dt1,
			1, 0, 0,
			1, 2 * dt1, 3 * dt1 * dt1;
		B << sxi - sxi_1, dx1, dx2;
		auto tmp = A.colPivHouseholderQr();
		data->coeffXu.row(i-1) = tmp.solve(B).transpose();

		B << syi - syi_1, dy1, dy2;
		data->coeffYu.row(i - 1) = tmp.solve(B).transpose();

		A << dt2, dt2* dt2, dt2* dt2* dt2,
			1, 0, 0,
			1, 2 * dt2, 3 * dt2 * dt2;
		tmp = A.colPivHouseholderQr();
		B << sxi1 - sxi, dx3, dx4;
		data->coeffXu.row(i) = tmp.solve(B).transpose();

		B << syi1- syi, dy3, dy4;
		data->coeffYu.row(i) = tmp.solve(B).transpose();

		//更新更新系数
	}
}
void CornerVertex(CanvasData* data, float dx_new, float dy_new) {
	//若选择的值型点为ti，左边为Si-1段，右边为Si段
	int i = -1;

	if (data->positive) {//更新Si段左端点导数
		i = data->edit_idx;
		data->dx_p[i] = dx_new;
		data->dy_p[i] = dy_new;
	}
	else {//更新Si-1段右端点导数
		i = data->edit_idx - 1;
		data->dx_n[i] = dx_new;
		data->dy_n[i] = dy_new;
	}
	//只修改一条曲线，解两个方程组
	float dt = data->tu[i + 1] - data->tu[i];
	float sxi = data->points[i][0];
	float syi = data->points[i][1];

	float sxi1 = data->points[i + 1][0];
	float syi1 = data->points[i + 1][1];

	float dx = data->dx_p[i];
	float dy = data->dy_p[i];
	float dx1 = data->dx_n[i];
	float dy1 = data->dy_n[i];

	Matrix3f A;
	Vector3f coeff;
	Vector3f B;
	A << dt, dt* dt, dt* dt* dt,
		1, 0, 0,
		1, 2 * dt, 3 * dt * dt;
	B << sxi1 - sxi, dx, dx1;
	auto tmp = A.colPivHouseholderQr();
	data->coeffXu.row(i) = tmp.solve(B).transpose();
	B << syi1 - syi, dy, dy1;
	data->coeffYu.row(i) = tmp.solve(B).transpose();
}

void CanvasSystem::OnUpdate(Ubpa::UECS::Schedule& schedule) {
	schedule.RegisterCommand([](Ubpa::UECS::World* w) {
		auto data = w->entityMngr.GetSingleton<CanvasData>();
		if (!data)
			return;

		if (ImGui::Begin("Canvas")) {
			ImGui::Checkbox("Enable grid", &data->opt_enable_grid);
			ImGui::Checkbox("Enable context menu", &data->opt_enable_context_menu);
			ImGui::Text("Mouse Left: drag to add lines,\nMouse Right: drag to scroll, click for context menu.");

			// Typically you would use a BeginChild()/EndChild() pair to benefit from a clipping region + own scrolling.
			// Here we demonstrate that this can be replaced by simple offsetting + custom drawing + PushClipRect/PopClipRect() calls.
			// To use a child window instead we could use, e.g:
			//      ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));      // Disable padding
			//      ImGui::PushStyleColor(ImGuiCol_ChildBg, IM_COL32(50, 50, 50, 255));  // Set a background color
			//      ImGui::BeginChild("canvas", ImVec2(0.0f, 0.0f), true, ImGuiWindowFlags_NoMove);
			//      ImGui::PopStyleColor();
			//      ImGui::PopStyleVar();
			//      [...]
			//      ImGui::EndChild();

			// Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
			ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
			ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
			if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
			if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
			ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

			// Draw border and background color
			ImGuiIO& io = ImGui::GetIO();
			ImDrawList* draw_list = ImGui::GetWindowDrawList();
			draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
			draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

			// This will catch our interactions
			ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
			const bool is_hovered = ImGui::IsItemHovered(); // Hovered
			const bool is_active = ImGui::IsItemActive();   // Held
			const ImVec2 origin(canvas_p0.x + data->scrolling[0], canvas_p0.y + data->scrolling[1]); // Lock scrolled origin
			const pointf2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

			//未开启edit点击为添加值点，开启则检测是否，修改值型点或值型点平滑属性
			if (!data->is_editing && is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left) ) {
				data->points.push_back(mouse_pos_in_canvas);
				data->updatad = true;
			}
			else if (data->is_editing && is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
				if (data->points.size() > 1) {
					CheckPoint(data, origin);
				}
			}

			//拖动值型点 || 更改值型点平滑属性
			if (data->is_darging) {
				data->points[data->draging_idx] = mouse_pos_in_canvas;
				data->updatad = true;
				//data->is_darging = false;
			}else if(data->changing_point) {
				float len = 0.1f;
				//拖动
				//dx_new = （io.mousepos - 值型点.x)/len;
				float dx_new = (io.MousePos.x - (origin.x + data->points[data->edit_idx][0])) / len;
				float dy_new = (io.MousePos.y - (origin.y + data->points[data->edit_idx][1])) / len;
				if (!data->positive) {//反向的控制点
					dx_new = -dx_new;
					dy_new = -dy_new;
				}
				if (data->edit_type == 0) {
					SmoothVertex(data, dx_new, dy_new);
				}
				else if (data->edit_type == 1) {
					LineVertex(data, dx_new, dy_new);
				}
				else {
					CornerVertex(data, dx_new, dy_new);
				}
				sampling(data, data->tu, 1000, data->coeffXu, data->coeffYu, true);
			}
			
			if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
				data->changing_point = false;
				data->is_darging = false;
			}

			//添加值型点之后，进行重新绘制三次样条
			if (data->points.size() > 1) {
				if (data->updatad) {
					data->updatad = false;

					//使用x，y坐标作为控制点
					vector<float> x(data->points.size());
					vector<float> y(data->points.size());
					for (int i = 0; i < data->points.size(); ++i) {
						x[i] = data->points[i][0];
						y[i] = data->points[i][1];
					}
					MatrixXf coeff = spline(x, y);
					sampling(data, 1000, coeff);
					data->draw_curve = true;
					//使用参数坐标作为控制点
						//均匀化参数
					vector<float> t = uniform_parameterization(data->points);
					data->tu = t;
					data->coeffXu = spline(t, x);
					data->coeffYu = spline(t, y);
					GenControlPoint(data);

					sampling(data, t, 1000, data->coeffXu,data->coeffYu, true);
					data->draw_curve1 = true;
						//弦长参数
					t = chordal_parameterization(data->points);
					data->coeffXc = spline(t, x);
					data->coeffYc = spline(t, y);
					sampling(data, t, 1000, data->coeffXc, data->coeffYc, false);
					data->draw_curve2 = true;
				}
			}

			

			//  mouse threshold when there's no context menu)
			// You may decide to make that threshold dynamic based on whether the mouse is hovering something etc.
			const float mouse_threshold_for_pan = data->opt_enable_context_menu ? -1.0f : 0.0f;
			if (is_active && ImGui::IsMouseDragging(ImGuiMouseButton_Right, mouse_threshold_for_pan))
			{
				data->scrolling[0] += io.MouseDelta.x;
				data->scrolling[1] += io.MouseDelta.y;
			}

			// Context menu (under default mouse threshold)
			ImVec2 drag_delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
			if (data->opt_enable_context_menu && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && drag_delta.x == 0.0f && drag_delta.y == 0.0f)
				ImGui::OpenPopupContextItem("context");

			//清除数据位置
			if (ImGui::BeginPopup("context"))
			{
				data->adding_line = false;
				if (ImGui::MenuItem("Remove one", NULL, false, data->points.size() > 0)) {
					data->points.pop_back();
					data->updatad = true;
				}
				if (ImGui::MenuItem("Remove all", NULL, false, data->points.size() > 0)) {
					data->points.clear();
					data->sample_points.clear();
					data->sample_points1.clear();
					data->sample_points2.clear();
					data->tu.clear();
					data->dx_p.clear();
					data->dy_p.clear();
					data->dx_n.clear();
					data->dy_n.clear();
					data->draw_curve = false;
					data->draw_curve1 = false;
					data->draw_curve2 = false;
				}
				ImGui::Checkbox("edit vertex", &data->is_editing);
				//更改值型点平滑属性
				const char* type[] = { "smooth", "line", "corner" };
				for (int i = 0; i < 3; ++i) {
					if (ImGui::Selectable(type[i], data->edit_type == i)) {
						data->edit_type = i;
					}
				}
				ImGui::EndPopup();
			}
			if (ImGui::Begin("debug")) {
				float len = 0.1f;
				if (data->points.size() > 1) {
					//n个性质点，对应n-1段曲线，每段曲线两端各绘制一个UI线段
					for (int i = 0; i < data->points.size() - 1; ++i) {
						//导数方向，目前只修改均匀参数化曲线
						float dx1 = data->dx_p[i];
						float dy1 = data->dy_p[i];
						ImVec2 p1(origin.x + data->points[i][0] + dx1 * len, origin.y + data->points[i][1] + dy1 * len);

						float dx2 = data->dx_n[i];
						float dy2 = data->dy_n[i];
						ImVec2 p2(origin.x + data->points[i + 1][0] - dx2 * len, origin.y + data->points[i + 1][1] - dy2 * len);
						ImGui::Text("idx: %d, p1(%f, %f), p2(%f, %f)", i, p1.x, p1.y, p2.x, p2.y);
					}
				}
			}
			ImGui::End();
			//控制面板
			if (ImGui::Begin("control panle")) {
				ImGui::Text("(%.2f, %.2f)", io.MousePos.x, io.MousePos.y);
				ImGui::Checkbox("normal", &data->normal);
				ImGui::Checkbox("uniform", &data->unifrom);
				ImGui::Checkbox("choral", &data->choral);
				for (int i = 0; i < data->points.size(); ++i) {
					ImGui::Text("(%.2f, %.2f)", origin.x + data->points[i][0], origin.y+data->points[i][1]);
				}
			}
			ImGui::End();

			// Draw grid + all lines in the canvas
			draw_list->PushClipRect(canvas_p0, canvas_p1, true);
			if (data->opt_enable_grid)
			{
				const float GRID_STEP = 64.0f;
				for (float x = fmodf(data->scrolling[0], GRID_STEP); x < canvas_sz.x; x += GRID_STEP)
					draw_list->AddLine(ImVec2(canvas_p0.x + x, canvas_p0.y), ImVec2(canvas_p0.x + x, canvas_p1.y), IM_COL32(200, 200, 200, 40));
				for (float y = fmodf(data->scrolling[1], GRID_STEP); y < canvas_sz.y; y += GRID_STEP)
					draw_list->AddLine(ImVec2(canvas_p0.x, canvas_p0.y + y), ImVec2(canvas_p1.x, canvas_p0.y + y), IM_COL32(200, 200, 200, 40));
			}
			/*for (int n = 0; n < data->points.size(); n += 2)
				draw_list->AddLine(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n + 1][0], origin.y + data->points[n + 1][1]), IM_COL32(255, 255, 0, 255), 2.0f);
			draw_list->PopClipRect();*/

			//绘制位置
			for (int n = 0; n < data->points.size(); ++n) {
				draw_list->AddCircleFilled(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), 3.f, IM_COL32(255, 0, 0, 255));
			}
			//修改值型点时显示UI
			if (data->is_editing) {
				float len = 0.1f;
				if (data->points.size() > 1) {
					//n个性质点，对应n-1段曲线，每段曲线两端各绘制一个UI线段
					for (int i = 0; i < data->points.size() - 1; ++i) {
						//导数方向，目前只修改均匀参数化曲线
						float dx1 = data->dx_p[i];
						float dy1 = data->dy_p[i];
						ImVec2 p1(origin.x + data->points[i][0] + dx1 * len, origin.y + data->points[i][1] + dy1 * len);
						draw_list->AddRectFilled(ImVec2(p1.x - 2.f, p1.y - 2.f), ImVec2(p1.x + 2.f, p1.y + 2.f), IM_COL32(255, 255, 255, 255));
						draw_list->AddLine(ImVec2(origin.x + data->points[i][0], origin.y + data->points[i][1]), p1, IM_COL32(255, 255, 255, 255));

						float dx2 = data->dx_n[i];
						float dy2 = data->dy_n[i];
						ImVec2 p2(origin.x + data->points[i + 1][0] - dx2 * len, origin.y + data->points[i + 1][1] - dy2 * len);
						draw_list->AddRectFilled(ImVec2(p2.x - 2.f, p2.y - 2.f), ImVec2(p2.x + 2.f, p2.y + 2.f), IM_COL32(255, 255, 255, 255));
						draw_list->AddLine(p2, ImVec2(origin.x + data->points[i + 1][0], origin.y + data->points[i + 1][1]), IM_COL32(255, 255, 255, 255));
					}
				}
			}

			if (data->normal && data->draw_curve) {
				for (int n = 0; n < data->sample_points.size() - 1; ++n) {
					draw_list->AddLine(ImVec2(origin.x + data->sample_points[n][0], origin.y + data->sample_points[n][1]), ImVec2(origin.x + data->sample_points[n + 1][0], origin.y + data->sample_points[n + 1][1]), IM_COL32(255, 0, 0, 255), 2.0f);
				}
			}
			if (data->unifrom && data->draw_curve1) {
				for (int n = 0; n < data->sample_points.size() - 1; ++n) {
					draw_list->AddLine(ImVec2(origin.x + data->sample_points1[n][0], origin.y + data->sample_points1[n][1]), ImVec2(origin.x + data->sample_points1[n + 1][0], origin.y + data->sample_points1[n + 1][1]), IM_COL32(0, 255, 0, 255), 2.0f);
				}
			}
			if (data->choral && data->draw_curve2) {
				for (int n = 0; n < data->sample_points.size() - 1; ++n) {
					draw_list->AddLine(ImVec2(origin.x + data->sample_points2[n][0], origin.y + data->sample_points2[n][1]), ImVec2(origin.x + data->sample_points2[n + 1][0], origin.y + data->sample_points2[n + 1][1]), IM_COL32(0, 0, 255, 255), 2.0f);
				}
			}
			if (data->is_darging) {
				draw_list->AddCircleFilled(io.MousePos, 6.f, IM_COL32(255, 255, 0, 255));
			}
		}

		ImGui::End();
		});
}