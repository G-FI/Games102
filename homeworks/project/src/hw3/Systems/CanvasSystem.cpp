#include "CanvasSystem.h"

#include "../Components/CanvasData.h"

#include <_deps/imgui/imgui.h>

#include<Eigen/Dense>
#include<Eigen/Core>
#include<numeric>

void uniform_parameterization(CanvasData* data);
void chordal_parameterization(CanvasData* data);

using namespace Ubpa;

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

			//添加数据位置
			if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
				data->points.push_back(mouse_pos_in_canvas);
			}

			//参数化 + 采样
			if (ImGui::Begin("control panle")) {
				//均匀参数化
				if (ImGui::Button("uniform parameterization")){
					uniform_parameterization(data);
				}
				if (ImGui::Button("Chordal parameterization")) {
					chordal_parameterization(data);
				}
			}
			ImGui::End();

		
		
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
				}
				if (ImGui::MenuItem("Remove all", NULL, false, data->points.size() > 0)) { 
					data->points.clear(); 
					data->sample_points.clear();
					data->uniform_param = false;
					data->sample_points2.clear();
					data->chordal_param = false;
				}
				ImGui::EndPopup();
			}

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

			if (data->uniform_param) {
				for (int n = 0; n < data->sample_points.size() - 1; ++n) {
					draw_list->AddLine(ImVec2(origin.x + data->sample_points[n][0], origin.y + data->sample_points[n][1]), ImVec2(origin.x + data->sample_points[n + 1][0], origin.y + data->sample_points[n + 1][1]), IM_COL32(255, 255, 0, 255), 2.0f);
				}
			}
			if (data->chordal_param) {
				for (int n = 0; n < data->sample_points2.size() - 1; ++n) {
					draw_list->AddLine(ImVec2(origin.x + data->sample_points2[n][0], origin.y + data->sample_points2[n][1]), ImVec2(origin.x + data->sample_points2[n + 1][0], origin.y + data->sample_points2[n + 1][1]), IM_COL32(0, 255, 0, 255), 2.0f);
				}
			}
		}

		ImGui::End();
	});
}

void chordal_parameterization(CanvasData* data) {
	int N = data->points.size();
	if (N > 1) {
		auto length = [](const auto& p1, const auto& p2)->float {
			return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]));
		};

		//限定采样点 t在[0, 1]之间
		float start_t = 0.f;
		
		std::vector<float> t(data->points.size());
		t[0] = start_t;
		float total_chrod = 0.f;
		std::vector<float> chord_length(data->points.size() - 1);
		for (int i = 0; i < data->points.size() - 1; ++i) {
			chord_length[i] = length(data->points[i + 1], data->points[i]);
			total_chrod += chord_length[i];
		}
		for (int i = 0; i < chord_length.size(); ++i) {
			t[i + 1] = chord_length[i] / total_chrod + t[i];
		}


		Eigen::VectorXd coefficientX(N);
		Eigen::VectorXd coefficientY(N);
		Eigen::MatrixXd A(N, N);

		Eigen::VectorXd X(N);
		Eigen::VectorXd Y(N);

		for (int i = 0; i < N; ++i) {
			float x = data->points[i][0];
			float y = data->points[i][1];
			X(i) = x;
			Y(i) = y;
			for (int j = 0; j < N; ++j) {
				A(i, j) = std::pow(t[i], j);
			}
		}
		coefficientX = A.colPivHouseholderQr().solve(X);
		coefficientY = A.colPivHouseholderQr().solve(Y);

		auto calculateX = [&](float t)->float {
			float res = 0.f;
			for (int i = 0; i < coefficientX.size(); ++i) {
				res += coefficientX[i] * std::pow(t, i);
			}
			return res;
		};

		auto calculateY = [&](float t)->float {
			float res = 0.f;
			for (int i = 0; i < coefficientX.size(); ++i) {
				res += coefficientY[i] * std::pow(t, i);
			}
			return res;
		};

		//采样
		start_t -= 0.1f;
		for (int i = 0; i < data->sample_num; ++i) {
			float t = start_t + i * data->stride;
			float x = calculateX(t);
			float y = calculateY(t);
			data->sample_points2.push_back(Ubpa::pointf2{ x,  y });
		}
		data->chordal_param = true;
	}
}

void uniform_parameterization(CanvasData* data) {
	int N = data->points.size();
	if (N > 0) {
		float start_t = 0.f;
		float t;

		Eigen::VectorXd coefficientX(N);
		Eigen::VectorXd coefficientY(N);
		Eigen::MatrixXd A(N, N);

		Eigen::VectorXd X(N);
		Eigen::VectorXd Y(N);

		for (int i = 0; i < N; ++i) {
			float x = data->points[i][0];
			float y = data->points[i][1];
			X(i) = x;
			Y(i) = y;
			t = (float)i / N ;

			for (int j = 0; j < N; ++j) {
				A(i, j) = std::pow(t, j);
			}
		}
		coefficientX = A.colPivHouseholderQr().solve(X);
		coefficientY = A.colPivHouseholderQr().solve(Y);

		auto calculateX = [&](float t)->float {
			float res = 0.f;
			for (int i = 0; i < coefficientX.size(); ++i) {
				res += coefficientX[i] * std::pow(t, i);
			}
			return res;
		};

		auto calculateY = [&](float t)->float {
			float res = 0.f;
			for (int i = 0; i < coefficientX.size(); ++i) {
				res += coefficientY[i] * std::pow(t, i);
			}
			return res;
		};

		//采样
		start_t -= 0.1f;
		for (int i = 0; i < data->sample_num; ++i) {
			float t = start_t + i * data->stride;
			float x = calculateX(t);
			float y = calculateY(t);
			data->sample_points.push_back(Ubpa::pointf2{ x,  y });
		}

		data->uniform_param = true;
	}
}

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


void powerbasefit(CanvasData* data) {
	int N = data->points.size();
	if (N > 0) {
		Eigen::VectorXd coefficient(N);
		Eigen::MatrixXd A(N, N);
		Eigen::VectorXd Y(N);

		//幂基矩阵
		for (int i = 0; i < N; ++i) {
			float x = data->points[i][0];
			float y = data->points[i][1];
			Y(i) = y;
			for (int j = 0; j < N; ++j) {
				A(i, j) = std::pow(x, j);
			}
		}
		//
		coefficient = A.colPivHouseholderQr().solve(Y);

		auto calculateY = [&](float x)->float {
			float res = 0.f;
			for (int i = 0; i < coefficient.size(); ++i) {
				res += coefficient[i] * std::pow(x, i);
			}
			return res;
		};

		//采样
		float start_x = data->points[0][0] - 50.f;
		for (int i = 0; i < data->sample_num; ++i) {
			float x = start_x + i * data->stride;
			float y = calculateY(x);
			data->sample_points.push_back(Ubpa::pointf2{  x,  y });
		}
	}
}

//void gaussbasefit(CanvasData* data) {
//
//
//	//点拟合函数
//	int N = data->points.size();
//	if (N > 0) {
//		Eigen::VectorXd coefficient(N);
//		Eigen::MatrixXd A(N, N);
//		Eigen::VectorXd Y(N);
//
//		float delt = data->delt;
//		//计算矩阵中的一个系数
//		auto normalDistrb = [&](int x, int xi)->float {
//			return std::exp(-((x - xi) * (x - xi)) / (2 * delt * delt));
//		};
//		//初始化μ（为采样点X值）
//		for (int i = 1; i < N; ++i) {
//			data->interploate_xi.push_back(data->points[i][0]);
//		}
//		//gauss基矩阵
//		for (int i = 0; i < N; ++i) {
//			float x = data->points[i][0];
//			float y = data->points[i][1];
//			Y(i) = y;
//			A(i, 0) = 1;
//			for (int j = 1; j < N; ++j) {
//				A(i, j) = normalDistrb(x, data->interploate_xi[j - 1]);
//			}
//		}
//		//求解
//		coefficient = A.colPivHouseholderQr().solve(Y);
//		/*{
//			std::stringstream ss;
//			ss << coefficient;
//			spdlog::info("coefficient matrix: " + ss.str());
//		}*/
//
//		//根据拟合的函数进行计算
//		auto normalDistribCalulate = [&](int x)->float {
//			float res = coefficient(0); //b0常数项
//			//N个插值点，对应n个高斯函数项
//			for (int i = 1; i < N; ++i) {
//				float xi = data->interploate_xi[i - 1];
//				res += coefficient(i) * std::exp(-((x - xi) * (x - xi)) / (2 * delt * delt));
//			}
//			return res;
//		};
//
//		//采样
//		float start_x = data->points[0][0] - 50.f;
//		for (int i = 0; i < data->sample_num; ++i) {
//			float x = start_x + i * data->stride;
//			float y = normalDistribCalulate(x);
//			data->sample_points.push_back(Ubpa::pointf2{ x,  y });
//		}
//		data->gaussbase_draw = true;
//
//	}
//}
