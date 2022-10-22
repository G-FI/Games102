#pragma once

#include <UGM/UGM.h>
#include<_deps/imgui/imgui.h>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<vector>

struct CanvasData {
	std::vector<Ubpa::pointf2> points;
	Ubpa::valf2 scrolling{ 0.f,0.f };
	bool opt_enable_grid{ true };
	bool opt_enable_context_menu{ true };
	bool adding_line{ false };
	
	//added
	std::vector<Ubpa::pointf2> sample_points;	//使用x,y
	std::vector<Ubpa::pointf2> sample_points1;	//使用均匀采样
	std::vector<Ubpa::pointf2> sample_points2;	//弦长参数化
	
	//三次样条参数,分别对应均匀参数化和弦长参数化
	Eigen::MatrixXf coeffXu;
	Eigen::MatrixXf coeffYu;
	Eigen::MatrixXf coeffXc;
	Eigen::MatrixXf coeffYc;

	//uniform 值型点参数
	std::vector<float> tu;
	//存储导数
	std::vector<float> dx_p;
	std::vector<float> dy_p;
	std::vector<float> dx_n;
	std::vector<float> dy_n;

	bool updatad{ false };

	bool normal{ false };
	bool unifrom{ true };
	bool choral{ false };


	bool draw_curve{ false };//使用x，y
	bool draw_curve1{ false };
	bool draw_curve2{ false };

	//更改值型点
	bool is_darging{ false };
	int draging_idx{-1};
	

	bool is_editing{ false };
	int edit_idx{ -1 };
	bool positive{ true }; //选择的控制点是正向的

	//edit type表示更新顶点的类型0：平滑，2：直线， 3：角部顶点是·
	int edit_type = 0;

	//是否更改顶点平滑属性
	bool changing_point{ false };
	
};


#include "details/CanvasData_AutoRefl.inl"
