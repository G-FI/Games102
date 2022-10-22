#pragma once

#include"imgui.h"
#include <vector>

struct CanvasData {
	std::vector<ImVec2> sample_data;

	bool draw_line = { false };
	float draw_data_num = {1000};
	std::vector<ImVec2> draw_data;
};