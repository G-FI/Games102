#include "DenoiseSystem.h"

#include "../Components/DenoiseData.h"

#include <_deps/imgui/imgui.h>

#include <spdlog/spdlog.h>
#include<sstream>

using namespace Ubpa;

void DenoiseSystem::OnUpdate(Ubpa::UECS::Schedule& schedule) {
	schedule.RegisterCommand([](Ubpa::UECS::World* w) {
		auto data = w->entityMngr.GetSingleton<DenoiseData>();
		if (!data)
			return;

		if (ImGui::Begin("Denoise")) {
			if (ImGui::Button("Mesh to HEMesh")) {
				data->heMesh->Clear();
				[&]() {
					if (!data->mesh) {
						spdlog::warn("mesh is nullptr");
						return;
					}

					if (data->mesh->GetSubMeshes().size() != 1) {
						spdlog::warn("number of submeshes isn't 1");
						return;
					}

					data->copy = *data->mesh;

					std::vector<size_t> indices(data->mesh->GetIndices().begin(), data->mesh->GetIndices().end());
					data->heMesh->Init(indices, 3);
					if (!data->heMesh->IsTriMesh())
						spdlog::warn("HEMesh init fail");
					
					for (size_t i = 0; i < data->mesh->GetPositions().size(); i++)
						data->heMesh->Vertices().at(i)->position = data->mesh->GetPositions().at(i);

					spdlog::info("Mesh to HEMesh success");
				}();
			}

			if (ImGui::Button("Add Noise")) {
				[&]() {
					if (!data->heMesh->IsTriMesh()) {
						spdlog::warn("HEMesh isn't triangle mesh");
						return;
					}

					for (auto* v : data->heMesh->Vertices()) {
						v->position += data->randomScale * (
							2.f * Ubpa::vecf3{ Ubpa::rand01<float>(),Ubpa::rand01<float>() ,Ubpa::rand01<float>() } - Ubpa::vecf3{ 1.f }
						);
					}

					spdlog::info("Add noise success");
				}();
			}

			if (ImGui::Button("Set Normal to Color")) {
				[&]() {
					if (!data->mesh) {
						spdlog::warn("mesh is nullptr");
						return;
					}

					data->mesh->SetToEditable();
					const auto& normals = data->mesh->GetNormals();
					std::vector<rgbf> colors;
					for (const auto& n : normals)
						colors.push_back((n.as<valf3>() + valf3{ 1.f }) / 2.f);
					data->mesh->SetColors(std::move(colors));

					spdlog::info("Set Normal to Color Success");
				}();
			}

			if (ImGui::Button("HEMesh to Mesh")) {
				[&]() {
					if (!data->mesh) {
						spdlog::warn("mesh is nullptr");
						return;
					}

					if (!data->heMesh->IsTriMesh() || data->heMesh->IsEmpty()) {
						spdlog::warn("HEMesh isn't triangle mesh or is empty");
						return;
					}

					data->mesh->SetToEditable();

					const size_t N = data->heMesh->Vertices().size();
					const size_t M = data->heMesh->Polygons().size();
					std::vector<Ubpa::pointf3> positions(N);
					std::vector<uint32_t> indices(M * 3);
					for (size_t i = 0; i < N; i++)
						positions[i] = data->heMesh->Vertices().at(i)->position;
					for (size_t i = 0; i < M; i++) {
						auto tri = data->heMesh->Indices(data->heMesh->Polygons().at(i));
						indices[3 * i + 0] = static_cast<uint32_t>(tri[0]);
						indices[3 * i + 1] = static_cast<uint32_t>(tri[1]);
						indices[3 * i + 2] = static_cast<uint32_t>(tri[2]);
					}
					data->mesh->SetColors({});
					data->mesh->SetUV({});
					data->mesh->SetPositions(std::move(positions));
					data->mesh->SetIndices(std::move(indices));
					data->mesh->SetSubMeshCount(1);
					data->mesh->SetSubMesh(0, { 0, M * 3 });
					data->mesh->GenUV();
					data->mesh->GenNormals();
					data->mesh->GenTangents();

					spdlog::info("HEMesh to Mesh success");
				}();
			}

			if (ImGui::Button("Recover Mesh")) {
				[&]() {
					if (!data->mesh) {
						spdlog::warn("mesh is nullptr");
						return;
					}
					if (data->copy.GetPositions().empty()) {
						spdlog::warn("copied mesh is empty");
						return;
					}

					*data->mesh = data->copy;

					spdlog::info("recover success");
				}();
			}

			if (ImGui::Button("Denoise")) {
				/*
				* 找到边界上的顶点，保存不更新
					for i from 1 to k //迭代次数
						先保存顶点 saved_vertex
						创建更新顶点数组 new_vertex
						for v in every non-boundary v
							计算法线(寻找邻接顶点)
							延法线反向更新
							保存到新顶点数组中

						new_vertex 更新data->heMesh中的顶点
				*/
				/*顶点在data->heMesh->Vertices()的position属性中*/

				/*auto v = data->heMesh->Vertices()[0];
				std::stringstream ss;
				ss << "邻接多边形个数: " << v->AdjPolygons().size() << std::endl;
				ss << "邻接顶点个数: " << v->AdjVertices().size() << std::endl;;
				ss << "邻接边个数: " << v->AdjEdges().size() << std::endl;
				spdlog::info(ss.str());*/
				//
				for (auto* p : data->heMesh->Vertices()) {
					if (p->IsOnBoundary()) {
						auto he = p->HalfEdge();
						p->is_boundary = true;
					}
				}

				//for (int i = 0; i < data->k; ++i) {
				//	auto vertices = data->heMesh->Vertices();	//vertices的元素指向的是底层hemesh的顶点
				//	std::vector<Ubpa::pointf3> new_positions(vertices.size());

				//	for (int j = 0; j < vertices.size(); ++j) {
				//		if (vertices[j]->is_boundary) {
				//			new_positions[j] = vertices[j]->position;
				//			continue;
				//		}

				//		//只是在一个三角形中的邻接点
				//		auto start_he = vertices[j]->HalfEdge();
				//		auto he = start_he;
				//		
				//		Ubpa::vecf3 dir;
				//		float area = 0.f;

				//		do {
				//			auto he1 = he->Pair()->Next();
				//			auto he2 = he->RotateNext();
				//			auto v = he->End();
				//			area = v.area();

				//			auto v1 = he1->End();
				//			auto v2 = he2->End();

				//			dir += (v1->cot() + v2->cot()) * (vertices[j]->position - v->position);
				//			he = he2;
				//		} while (he != start_he);
				//	
				//		new_positions[j] = (vertices[j] + data->lambda * dir)/4*area;
				//	}

				//	for (int j = 0; j < vertices.size(); ++j) {
				//		vertices[j]->position = new_positions[j];
				//	}
				//}
			}
		}
		ImGui::End();
	});
}
