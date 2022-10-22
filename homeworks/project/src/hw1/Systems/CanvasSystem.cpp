#include "CanvasSystem.h"

#include "../Components/CanvasData.h"

#include <_deps/imgui/imgui.h>
#include<string>
#include <spdlog/spdlog.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<vector>
#include<sstream>

using namespace Ubpa;

void CanvasSystem::OnUpdate(Ubpa::UECS::Schedule& schedule) {
//{
//		schedule.RegisterCommand([](Ubpa::UECS::World* w) {
//			auto data = w->entityMngr.GetSingleton<CanvasData>();
//			if (!data)
//				return;
//
//			if (ImGui::Begin("Canvas")) {
//				ImGui::Checkbox("Enable grid", &data->opt_enable_grid);
//				ImGui::Checkbox("Enable context menu", &data->opt_enable_context_menu);
//				ImGui::Text("Mouse Left: drag to add lines,\nMouse Right: drag to scroll, click for context menu.");
//
//				// Typically you would use a BeginChild()/EndChild() pair to benefit from a clipping region + own scrolling.
//				// Here we demonstrate that this can be replaced by simple offsetting + custom drawing + PushClipRect/PopClipRect() calls.
//				// To use a child window instead we could use, e.g:
//				//      ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));      // Disable padding
//				//      ImGui::PushStyleColor(ImGuiCol_ChildBg, IM_COL32(50, 50, 50, 255));  // Set a background color
//				//      ImGui::BeginChild("canvas", ImVec2(0.0f, 0.0f), true, ImGuiWindowFlags_NoMove);
//				//      ImGui::PopStyleColor();
//				//      ImGui::PopStyleVar();
//				//      [...]
//				//      ImGui::EndChild();
//
//				// Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
//				ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
//				ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
//				if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
//				if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
//				ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);
//
//				// Draw border and background color
//				ImGuiIO& io = ImGui::GetIO();
//				ImDrawList* draw_list = ImGui::GetWindowDrawList(); //我要绘制的图元
//				draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
//				draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));
//
//				// This will catch our interactions
//				ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
//				const bool is_hovered = ImGui::IsItemHovered(); // Hovered
//				const bool is_active = ImGui::IsItemActive();   // Held
//
//				//获取 mouse相对canvas面板左上角的坐标，data->scrolling是 canvas左上角与初始
//				const ImVec2 origin(canvas_p0.x + data->scrolling[0], canvas_p0.y + data->scrolling[1]); // Lock scrolled origin
//				const pointf2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);
//
//				// Add first and second point, 左键点击，并且在hovered，那么加入两个点
//				if (is_hovered && !data->adding_line && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
//				{
//					data->points.push_back(mouse_pos_in_canvas);
//					data->points.push_back(mouse_pos_in_canvas);
//					data->adding_line = true;
//				}
//				if (data->adding_line)
//				{
//					data->points.back() = mouse_pos_in_canvas;
//					if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
//						data->adding_line = false;
//				}
//
//				// Pan (we use a zero mouse threshold when there's no context menu)
//				// You may decide to make that threshold dynamic based on whether the mouse is hovering something etc.
//				const float mouse_threshold_for_pan = data->opt_enable_context_menu ? -1.0f : 1.f;
//				if (is_active && ImGui::IsMouseDragging(ImGuiMouseButton_Right, mouse_threshold_for_pan))
//				{
//					data->scrolling[0] += io.MouseDelta.x;
//					data->scrolling[1] += io.MouseDelta.y;
//				}
//
//				// Context menu (under default mouse threshold)
//				ImVec2 drag_delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
//				if (data->opt_enable_context_menu && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && drag_delta.x == 0.0f && drag_delta.y == 0.0f)
//					ImGui::OpenPopupContextItem("context");
//				if (ImGui::BeginPopup("context"))
//				{
//					if (data->adding_line)
//						data->points.resize(data->points.size() - 2);
//					data->adding_line = false;
//					if (ImGui::MenuItem("Remove one", NULL, false, data->points.size() > 0)) { data->points.resize(data->points.size() - 2); }
//					if (ImGui::MenuItem("Remove all", NULL, false, data->points.size() > 0)) { data->points.clear(); }
//					ImGui::EndPopup();
//				}
//
//				// Draw grid + all lines in the canvas
//				draw_list->PushClipRect(canvas_p0, canvas_p1, true);
//				if (data->opt_enable_grid)
//				{
//					const float GRID_STEP = 64.0f;
//					for (float x = fmodf(data->scrolling[0], GRID_STEP); x < canvas_sz.x; x += GRID_STEP)
//						draw_list->AddLine(ImVec2(canvas_p0.x + x, canvas_p0.y), ImVec2(canvas_p0.x + x, canvas_p1.y), IM_COL32(200, 200, 200, 40));
//					for (float y = fmodf(data->scrolling[1], GRID_STEP); y < canvas_sz.y; y += GRID_STEP)
//						draw_list->AddLine(ImVec2(canvas_p0.x, canvas_p0.y + y), ImVec2(canvas_p1.x, canvas_p0.y + y), IM_COL32(200, 200, 200, 40));
//				}
//				for (int n = 0; n < data->points.size(); n += 2)
//					draw_list->AddLine(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n + 1][0], origin.y + data->points[n + 1][1]), IM_COL32(255, 255, 0, 255), 2.0f);
//				/*for (int n = 0; n < data->points.size(); n += 4)
//					draw_list->AddQuad(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]),
//						ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), IM_COL32(255, 0, 0, 255));*/
//				draw_list->PopClipRect();
//			}
//
//			ImGui::End();
//			});
//}
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
			ImDrawList* draw_list = ImGui::GetWindowDrawList(); //我要绘制的图元
			draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
			draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

			//获取 mouse相对canvas面板左上角的坐标，data->scrolling是 canvas左上角与初始
			const ImVec2 origin(canvas_p0.x + data->scrolling[0], canvas_p0.y + data->scrolling[1]); // Lock scrolled origin
			const pointf2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);
			
			{
				bool is_hoverd, is_clicked;

				{//幂基插值法
					ImGui::Button("fit power base", ImVec2(40.f, 20.f));
					is_hoverd = ImGui::IsItemHovered();
					is_clicked = ImGui::IsItemClicked(ImGuiMouseButton_Left);
					if (is_hoverd && is_clicked) {
						{
							std::string msg{ ": " };
							for (const auto& point : data->points) {
								msg += "(";
								msg += std::to_string(point[0]);
								msg += " ";
								msg += std::to_string(point[1]);
								msg += "), ";
							}
							spdlog::info(std::to_string(data->points.size()) + msg);
						}
						//点拟合函数
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
							std::stringstream ss;
							ss << coefficient;
							spdlog::info("coefficient matrix: " + ss.str());
							auto calculateY = [&](float x)->float {
								float res = 0.f;
								for (int i = 0; i < coefficient.size(); ++i) {
									res += coefficient[i] * std::pow(x, i);
								}
								return res;
							};

							//采样
							pointf2 startPoint = data->points[0];
							data->sample_points.push_back(ImVec2(origin.x + startPoint[0], origin.y + startPoint[1]));
							for (int i = 1; i < data->sample_num; ++i) {
								float x = startPoint[0] + i * data->stride;
								float y = calculateY(x);
								data->sample_points.push_back(ImVec2(origin.x + x, origin.y + y));
							}
							data->fit_draw = true;
						}
					}
				}
				
				{//高斯基插值法
					ImGui::Button("fit Gausss base", ImVec2(40.f, 20.f));
					is_hoverd = ImGui::IsItemHovered();
					is_clicked = ImGui::IsItemClicked(ImGuiMouseButton_Left);
					if (is_hoverd && is_clicked) {
						{
							std::string msg{ ": " };
							for (const auto& point : data->points) {
								msg += "(";
								msg += std::to_string(point[0]);
								msg += " ";
								msg += std::to_string(point[1]);
								msg += "), ";
							}
							spdlog::info(std::to_string(data->points.size()) + msg);
						}
						//点拟合函数
						int N = data->points.size();
						if (N > 0) {
							Eigen::VectorXd coefficient(N);
							Eigen::MatrixXd A(N, N);
							Eigen::VectorXd Y(N);

							float delt = data->delt;
							//计算矩阵中的一个系数
							auto normalDistrb = [&](int x, int xi)->float {
								return std::exp(-((x - xi) * (x - xi)) / (2 * delt * delt));
							};
							//初始化μ（为采样点X值）
							for (int i = 1; i < N; ++i) {
								data->interploate_xi.push_back(data->points[i][0]);
							}
							//gauss基矩阵
							for (int i = 0; i < N; ++i) {
								float x = data->points[i][0];
								float y = data->points[i][1];
								Y(i) = y;
								A(i, 0) = 1;
								for (int j = 1; j < N; ++j) {
									A(i, j) = normalDistrb(x, data->interploate_xi[j - 1]);
								}
							}
							//求解
							coefficient = A.colPivHouseholderQr().solve(Y);
							{
								std::stringstream ss;
								ss << coefficient;
								spdlog::info("coefficient matrix: " + ss.str());
							}

							//根据拟合的函数进行计算
							auto normalDistribCalulate = [&](int x)->float {
								float res = coefficient(0); //b0常数项
								//N个插值点，对应n个高斯函数项
								for (int i = 1; i < N; ++i) {
									float xi = data->interploate_xi[i - 1];
									res += coefficient(i) * std::exp(-((x - xi) * (x - xi)) / (2 * delt * delt));
								}
								return res;
							};

							//采样
							pointf2 startPoint = data->points[0];
							data->sample_points.push_back(ImVec2(origin.x + startPoint[0], origin.y + startPoint[1]));
							for (int i = 1; i < data->sample_num; ++i) {
								float x = startPoint[0] + i * data->stride;
								float y = normalDistribCalulate(x);
								data->sample_points.push_back(ImVec2(origin.x + x, origin.y + y));
							}
							data->fit_draw = true;
						}
					}
				}
				{//最小二乘拟合
					ImGui::Button("ls_fit", ImVec2(40.f, 20.f));
					is_hoverd = ImGui::IsItemHovered();
					is_clicked = ImGui::IsItemClicked();
					if (is_hoverd && is_clicked) {
						{
							std::string msg{ ": " };
							for (const auto& point : data->points) {
								msg += "(";
								msg += std::to_string(point[0]);
								msg += " ";
								msg += std::to_string(point[1]);
								msg += "), ";
							}
							spdlog::info(std::to_string(data->points.size()) + msg);
						}
						int N = data->points.size();//
						if (N > 0) {
							//TransposeA * A * coeffient = TransposeA * Y
							Eigen::VectorXf	coefficient(data->ls_m+1);
							Eigen::MatrixXf A(N, data->ls_m + 1);	//
							Eigen::VectorXf Y(N); //采样点Y值


							//初始化A, Y
							for (int i = 0; i < N; ++i) {
								float x = data->points[i][0];
								float y = data->points[i][1];
								Y(i) = y;
								for (int j = 0; j < data->ls_m + 1; ++j) {
									A(i, j) = std::pow(x, j);
								}
							}

							//计算 coeffient
							Eigen::MatrixXf TA = A.transpose();
							Eigen::VectorXf B =  TA* Y;
							Eigen::MatrixXf A1 = TA * A;
							Eigen::MatrixXf IA1 = A1.inverse();
							coefficient = IA1 * B;
							{
								std::stringstream ss;
								ss << "transpose A = " << TA << "\n";
								ss << "B = " << B << "\n";
								ss << "A1 = " << A1 << "\n";
								ss << "IA1 = " << IA1 << "\n";
								ss << "coefficient matrix: " << coefficient << "\n";
								spdlog::info(ss.str());
							}
							//为什么下面这个方法算出来不对？
							//coefficient = A1.colPivHouseholderQr().solve(B);
							/*{
								std::stringstream ss;
								ss << coefficient;
								spdlog::info("coefficient matrix: " + ss.str());
							}*/
							//采样
							auto calculateY = [&](float x)->float {
								float res = 0.f;
								for (int i = 0; i < coefficient.size(); ++i) {
									res += coefficient[i] * std::pow(x, i);
								}
								return res;
							};
							float startX = data->points[0][0] - 20.f;
							for (int i = 0; i < data->sample_num; ++i) {
								float x = startX + data->stride * i;
								float y = calculateY(x);
								data->ls_sample_points.push_back(ImVec2(origin.x + x, origin.y + y));
							}
							data->ls_drwa = true;
						}
					}
				}
			}
			{//清除拟合函数
				ImGui::Button("clear fit draw", ImVec2(40.f, 20.f));
				bool is_hoverd = ImGui::IsItemHovered();
				bool is_clicked = ImGui::IsItemClicked();
				if (is_hoverd && is_clicked) {
					//清除插值数据
					data->sample_points.clear();
					data->interploate_xi.clear();
					data->fit_draw = false;

					//清除最小二乘数据
					data->ls_sample_points.clear();
					data->ls_drwa = false;
				}
			}
			{//设置高斯函数标准差, 最小二乘最高次数 UI
				ImGui::SliderFloat("delt", &data->delt, 0.1f, 20.f, "%.1f", 1);
				ImGui::InputInt("LSS power", &data->ls_m);
			}

			// This will catch our interactions
			ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
			const bool is_hovered = ImGui::IsItemHovered(); // Hovered
			const bool is_active = ImGui::IsItemActive();   // Held


			// Add first and second point
			if (is_hovered && !data->adding_line && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
			{
				data->points.push_back(mouse_pos_in_canvas);
				//data->points.push_back(mouse_pos_in_canvas);
				data->adding_line = true;
			}

			if (data->adding_line)
			{
				data->points.back() = mouse_pos_in_canvas;
				if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
					data->adding_line = false;
			}

			// Pan (we use a zero mouse threshold when there's no context menu)
			// You may decide to make that threshold dynamic based on whether the mouse is hovering something etc.
			const float mouse_threshold_for_pan = data->opt_enable_context_menu ? -1.0f : 1.f;
			if (is_active && ImGui::IsMouseDragging(ImGuiMouseButton_Right, mouse_threshold_for_pan))
			{
				data->scrolling[0] += io.MouseDelta.x;
				data->scrolling[1] += io.MouseDelta.y;
			}

			// Context menu (under default mouse threshold)
			ImVec2 drag_delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
			if (data->opt_enable_context_menu && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && drag_delta.x == 0.0f && drag_delta.y == 0.0f)
				ImGui::OpenPopupContextItem("context");
			
			if (ImGui::BeginPopup("context"))
			{
				if (data->adding_line)
					data->points.resize(data->points.size() - 2);
				data->adding_line = false;
				if (ImGui::MenuItem("Remove one", NULL, false, data->points.size() > 0)) { data->points.resize(data->points.size() - 2); }
				if (ImGui::MenuItem("Remove all", NULL, false, data->points.size() > 0)) { 
					data->points.clear();
					data->sample_points.clear(); 
					data->fit_draw = false; 
					
					data->ls_drwa = false;
					data->ls_sample_points.clear();	
				}
				ImGui::EndPopup();
			}


			// Draw grid + all lines in the canvas
			draw_list->PushClipRect(canvas_p0, canvas_p1, true);
			if (data->opt_enable_grid)
			{
				const float grid_step = 64.0f;
				for (float x = fmodf(data->scrolling[0], grid_step); x < canvas_sz.x; x += grid_step)
					draw_list->AddLine(ImVec2(canvas_p0.x + x, canvas_p0.y), ImVec2(canvas_p0.x + x, canvas_p1.y), IM_COL32(200, 200, 200, 40));
				for (float y = fmodf(data->scrolling[1], grid_step); y < canvas_sz.y; y += grid_step)
					draw_list->AddLine(ImVec2(canvas_p0.x, canvas_p0.y + y), ImVec2(canvas_p1.x, canvas_p0.y + y), IM_COL32(200, 200, 200, 40));
			}

			for (int n = 0; n < data->points.size(); ++n) {
				draw_list->AddCircleFilled(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), 1.f, IM_COL32(255, 0, 0, 255));
			}

		
			for (int n = 0; n < (int)data->points.size() - 1; n += 1) {
				draw_list->AddLine(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n + 1][0], origin.y + data->points[n + 1][1]), IM_COL32(255, 255, 0, 255), 2.0f);
			}
			/*for (int n = 0; n < data->points.size(); n += 4)
				draw_list->AddQuad(ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]),
					ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), ImVec2(origin.x + data->points[n][0], origin.y + data->points[n][1]), IM_COL32(255, 0, 0, 255));*/
			
			//插值曲线绘制
			if (data->fit_draw) {
				draw_list->AddPolyline(data->sample_points.data(), data->sample_num, IM_COL32(0, 255, 0, 255), false, 1.f);
			}
			//最小二乘曲线绘制
			if (data->ls_drwa) {
				draw_list->AddPolyline(data->ls_sample_points.data(), data->sample_num, ImColor(255, 0, 0, 255), false, 1.f);
			}
			draw_list->PopClipRect();
		}

		ImGui::End();
		});
}
