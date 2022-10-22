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
	std::vector<Ubpa::pointf2> sample_points;	//ʹ��x,y
	std::vector<Ubpa::pointf2> sample_points1;	//ʹ�þ��Ȳ���
	std::vector<Ubpa::pointf2> sample_points2;	//�ҳ�������
	
	//������������,�ֱ��Ӧ���Ȳ��������ҳ�������
	Eigen::MatrixXf coeffXu;
	Eigen::MatrixXf coeffYu;
	Eigen::MatrixXf coeffXc;
	Eigen::MatrixXf coeffYc;

	//uniform ֵ�͵����
	std::vector<float> tu;
	//�洢����
	std::vector<float> dx_p;
	std::vector<float> dy_p;
	std::vector<float> dx_n;
	std::vector<float> dy_n;

	bool updatad{ false };

	bool normal{ false };
	bool unifrom{ true };
	bool choral{ false };


	bool draw_curve{ false };//ʹ��x��y
	bool draw_curve1{ false };
	bool draw_curve2{ false };

	//����ֵ�͵�
	bool is_darging{ false };
	int draging_idx{-1};
	

	bool is_editing{ false };
	int edit_idx{ -1 };
	bool positive{ true }; //ѡ��Ŀ��Ƶ��������

	//edit type��ʾ���¶��������0��ƽ����2��ֱ�ߣ� 3���ǲ������ǡ�
	int edit_type = 0;

	//�Ƿ���Ķ���ƽ������
	bool changing_point{ false };
	
};


#include "details/CanvasData_AutoRefl.inl"
