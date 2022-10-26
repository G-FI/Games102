#pragma once

#include"imgui.h"
#include <vector>
using namespace std;

struct CanvasData {
	//监测面板开关

	bool monitor_open = { false };

	//控制点
	vector<ImVec2> points;

	//细分迭代次数
	int k {0};
	float alpha{0.1};
	bool closure {false};


	//控制开关
	bool draw_line{ false };
	bool chaiukin2_draw{false};
	bool chaiukin3_draw{ false };
	bool fpsubdiv_draw{ false };

	bool is_editing{ false };
	bool is_draging{ false };
	int edit_idx { -1 };

	//细分后的点
	vector<ImVec2> chaiukin2_points;
	vector<ImVec2> chaiukin3_points;
	vector<ImVec2> fpsubdiv_points;



	void Reset() {
		monitor_open = false;
		points.clear();

		chaiukin2_points.clear();
		chaiukin3_points.clear();
		fpsubdiv_points.clear();
	}
};