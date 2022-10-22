#pragma once

#include <UGM/UGM.h>

struct CanvasData {
	std::vector<Ubpa::pointf2> points;
	Ubpa::valf2 scrolling{ 0.f,0.f };
	bool opt_enable_grid{ true };
	bool opt_enable_context_menu{ true };
	bool adding_line{ false };

	//added
	std::vector<Ubpa::pointf2> sample_points;
	float stride{ 1.2f/1000 };
	int sample_num{ 1000 };
	bool uniform_param{ false };

	bool chordal_param{ false };
	std::vector<Ubpa::pointf2> sample_points2;

};

#include "details/CanvasData_AutoRefl.inl"
