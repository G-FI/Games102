#pragma once

#include <UGM/UGM.h>
#include "imgui/imgui.h"
struct CanvasData {
	std::vector<Ubpa::pointf2> points;
	Ubpa::valf2 scrolling{ 0.f,0.f };
	bool opt_enable_grid{ true };
	bool opt_enable_context_menu{ true };
	bool adding_line{ false };

	//added
	std::vector<ImVec2> sample_points;
	bool fit_draw{ false };
	int sample_num = 10000;
	float stride = 1.f;

	//for gauss base
	std::vector<float> interploate_xi;//ʹ��gauss��ʱ���洢�Ħ̣�ȡ�Ĳ�ֵ��xֵ��
	float delt = 1.f;	//��׼��

	//for least square fit
	std::vector<ImVec2> ls_sample_points;
	int ls_m = 0;	//��С���˷�������������ߴ���
	bool ls_drwa{ false };
};

#include "details/CanvasData_AutoRefl.inl"
