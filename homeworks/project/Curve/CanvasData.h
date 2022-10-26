#pragma once

#include"imgui.h"
#include <vector>
using namespace std;

struct CanvasData {
	//监控面板开关
	bool monitor_open = { false };
	//控制点
	vector<ImVec2> points;

	//bezier
	vector<ImVec2> bezier_curve;
	int k = 3;
	bool bezier_draw = { false };

	//B-spline
	vector<ImVec2> b_curve;
	bool b_draw = { false };

	//每段采样数量
	int num = 1000;

	//编辑控制点
	bool editing = { false };
	bool is_draging = { false };
	int edit_idx = { -1 };

	void Reset() {
		monitor_open = false;
		points.clear();

		bezier_curve.clear();
		bezier_draw = false;

		b_curve.clear();
		b_draw = false;

		editing = false;
		is_draging = false;
		edit_idx = -1;
	}
};