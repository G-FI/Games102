#pragma once

#include"imgui.h"
#include <vector>
#include<random>
#include<Eigen/Core>


using namespace std;
using namespace Eigen;
struct CanvasData {

	bool monitor_open{ false };

	std::vector<ImVec2> points;


	//rbf 参数
	int epochs{ 10 };
	float lr{ 0.01 };
	int rbf_number{ 10 };



	bool draw_line = { false };
	bool is_editing{ false };
	bool is_draging{ false };
	int edit_idx{ -1 };

    bool trained{ false };

	std::vector<ImVec2> curve;

	void Reset() {
		epochs = 10;
		lr = 0.01;
		rbf_number = 10;
		draw_line = false;
		is_editing = false;
		is_draging = false;
		edit_idx = -1;
		points.clear();
		curve.clear();
	}
};

struct RBFNetwork {
    //super parameters
    int epochs_;
    float lr_;

    //forward时缓存的参数，在backward时需要用到
    float x, y, y_pred, loss;

    //Parameters
    ArrayXf W1, W2, B1;
    float b2;

    //对于Parameters的偏导数
    ArrayXf DW1, DW2, DB1;
    float db2;

    //存储forward时变量
    ArrayXf Z, A; //Z = W1 * x + B1,  A = gauss(Z), y_pred = W2 *A + b2
    RBFNetwork(int rbf_size, int epochs, float lr) {
        epochs_ = epochs;
        lr_ = lr;

        random_device e;
        normal_distribution<float> normal_dist(100, 10.f);
        auto random_init = [&](float dummy) {return normal_dist(e); };

        W1 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);
        W2 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);

        B1 = RowVectorXf::Zero(rbf_size).unaryExpr(random_init);

        b2 = random_init(0);

        DW1 = RowVectorXf::Zero(rbf_size);
        DW2 = RowVectorXf::Zero(rbf_size);
        DB1 = RowVectorXf::Zero(rbf_size);
        db2 = 0.f;
    }

    float forward(float x_input, float y_input) {
        //auto gauss_kernel = [](float x)->float {return exp(-0.5 * x * x) / sqrt(2 * 3.1415926); };
        x = x_input;
        y = y_input;

        Z = x * W1 + B1;
        //cout <<"Z = "<< Z << endl;
       // A =  ((-0.5f) * Z.pow(2) / sqrt(2 * 3.1415926)).exp();
        A = (-0.5f * Z.pow(2)).exp() / sqrt(2 * 3.1415926);
        //cout << "A = " << A << endl;

        y_pred = (A * W2).sum() + b2;
        //cout << "y_pred = " << y_pred << endl;
        return y_pred;
    }
    float get_loss() {
        return (y - y_pred) * (y - y_pred);
    }
    void backward() {
        loss = (y - y_pred) * (y - y_pred);

//        cout << "loss = " << loss << endl;

        float DYpred = -2 * (y - y_pred);

        DW2 = DYpred * A;
        //cout << "DW2 = " << DW2 << endl;

        db2 = DYpred * 1;
        //cout << "db2 = " << db2 << endl;
        ArrayXf DA = DYpred * W2;
        ArrayXf DZ = DA * A * (-1 * Z);

        DW1 = DZ * x;
        //cout << "DW1 = " << DW1 << endl;
        DB1 = DZ * 1;
        //cout << "DB1 = " << DB1 << endl;
    }

    void optimial() {
        W1 -= DW1 * lr_;
        B1 -= DB1 * lr_;
        W2 -= DW2 * lr_;
        b2 -= db2 * lr_;
    }

    vector<float> train(const vector<ImVec2>& data) {
        vector<float> ret_loss;
        ret_loss.reserve(data.size());
        for (int epoch = 0; epoch < epochs_; ++epoch) {
            for (const ImVec2& point : data) {
                forward(point.x, point.y);
                ret_loss.push_back(get_loss());
                backward();
                optimial();
            }
        }
        return ret_loss;
    }
};
