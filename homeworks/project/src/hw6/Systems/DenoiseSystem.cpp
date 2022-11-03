#include "DenoiseSystem.h"

#include "../Components/DenoiseData.h"

#include <_deps/imgui/imgui.h>

#include <spdlog/spdlog.h>
#include<sstream>

using namespace Ubpa;
//
//pointf3 Cross(const pointf3& a, const pointf3& b) {
//	pointf3 res;
//	res[0] = a[1] * b[2] - a[2] * b[1];
//	res[1] = a[2] * b[0] - a[0] * b[2];
//	res[2] = a[0] * b[1] - a[1] * b[0];
//	return res;
//}
//
//float Norm(const pointf3& a) {
//	return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
//}

pointf3 operator /(const pointf3& a, float b) {
	return pointf3(a[0] / b, a[1] / b, a[2] / b);
}

pointf3 operator *(float a, const pointf3& p) {
	return pointf3(a * p[0], a * p[1], a * p[2]);
}

pointf3 operator +(const pointf3& p1, const pointf3& p2) {
	return pointf3( p1[0] + p2[0],  p1[1] * p2[1], p1[2] * p2[2]);
}

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
				
				
				//v1为尖点
				auto Cot = [](const pointf3& v1, const pointf3& v2, const pointf3& v3)->float {
					vecf3 a = v2 - v1;
					vecf3 b = v3 - v1;
					return vecf3::cot_theta(a, b);
				};

				for (auto* p : data->heMesh->Vertices()) {
					if (p->IsOnBoundary()) {
						auto he = p->HalfEdge();
						p->is_boundary = true;
					}
				}

				auto vertices = data->heMesh->Vertices();	//vertices的元素指向的是底层hemesh的顶点
				std::vector<Ubpa::pointf3> new_positions(vertices.size());

				for (int j = 0; j < vertices.size(); ++j) {
					if (vertices[j]->is_boundary) {
						new_positions[j] = vertices[j]->position;
						continue;
					}

					//启示半边
					auto* start_he = vertices[j]->HalfEdge();
					auto* he = start_he;

					Ubpa::vecf3 Hn(0.f, 0.f, 0.f);
					float area = 0.f;

			
					do {
						auto* he1 = he->RotateNext();
						auto* he2 = he->RotatePre();

						assert(he1->Origin() == he->Origin());
						assert(he2->Origin() == he->Origin());

						auto* vi = he->Origin();
						auto* vj = he->End();
						auto* v_alpha = he1->End();
						auto* v_beta = he2->End();



						auto cot_alpha = Cot(v_alpha->position, vi->position, vj->position);
						auto cot_beta = Cot(v_beta->position, vi->position, vj->position);
						vecf3 tmp = vi->position- vj->position;
						area += (float)1 / 8 * (cot_alpha + cot_beta) * tmp.norm2();

						Hn += (cot_alpha + cot_beta) * tmp;
						he = he1;
					} while (he != start_he);
					new_positions[j] = vertices[j]->position + (data->lambda / (4 * area)) * Hn;


				}

				for (int j = 0; j < vertices.size(); ++j) {
					vertices[j]->position = new_positions[j];
				}


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
		}
		ImGui::End();
	});
}
