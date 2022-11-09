#include "DenoiseSystem.h"

#include "../Components/DenoiseData.h"

#include <_deps/imgui/imgui.h>

#include <spdlog/spdlog.h>
#include<sstream>
#define EPSILON 1E-4F

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
float GetTriangleArea(Vertex* v0, Vertex* v1, Vertex* v2) {
	return 0.5f * v0->position.distance(v1->position) * v0->position.distance(v2->position) * (v1->position - v0->position).sin_theta(v2->position - v0->position);
}

float GetAmixed(Vertex* v) {
	float A_mixed = 0.f;
	if (v->IsOnBoundary()) return A_mixed;

	for (auto* adjV : v->AdjVertices()) {
		Vertex* np;
		auto* adjVE = adjV->HalfEdge();
		while (true) {
			if (v == adjVE->Next()->End()) {
				np = adjVE->End();
				break;
			}
			else
				adjVE = adjVE->Pair()->Next();
		}

		if ((adjV->position - v->position).dot(np->position - v->position) >= 0.f && (v->position - adjV->position).dot(np->position - adjV->position) >= 0.f && (v->position - np->position).dot(adjV->position - np->position) >= 0.f) {
			//T is non-obtuse
			if (GetTriangleArea(v, adjV, np) > EPSILON)
				A_mixed += (v->position.distance2(adjV->position) * (v->position - np->position).cot_theta(adjV->position - np->position) + v->position.distance2(np->position) * (v->position - adjV->position).cot_theta(np->position - adjV->position)) / 8.f;
		}
		else {
			if ((adjV->position - v->position).dot(np->position - v->position) < 0.f) {
				//T is obtuse && the angle of T at v is obtuse
				A_mixed += GetTriangleArea(v, adjV, np) / 2.f;
			}
			else {
				//T is obtuse && the angle of T at v is non-obtuse
				A_mixed += GetTriangleArea(v, adjV, np) / 4.f;
			}
		}
	}

	return A_mixed;
}

valf3 GetMeanCurvatureOperator(Vertex* v) {
	valf3 c{ 0.f };
	if (v->IsOnBoundary()) return c;

	float A_mixed = GetAmixed(v);
	if (A_mixed < EPSILON) return c;

	for (auto* adjV : v->AdjVertices()) {
		Vertex* pp;
		Vertex* np;
		auto* adjVE = adjV->HalfEdge();
		while (true) {
			if (v == adjVE->End()) {
				pp = adjVE->Next()->End();
				break;
			}
			else
				adjVE = adjVE->Pair()->Next();
		}
		while (true) {
			if (v == adjVE->Next()->End()) {
				np = adjVE->End();
				break;
			}
			else
				adjVE = adjVE->Pair()->Next();
		}

		if (GetTriangleArea(v, adjV, pp) > EPSILON && GetTriangleArea(v, adjV, np) > EPSILON) {
			float cot_alpha = (adjV->position - pp->position).cot_theta(v->position - pp->position);
			float cot_beta = (adjV->position - np->position).cot_theta(v->position - np->position);
			c += (cot_alpha + cot_beta) * (v->position - adjV->position);
		}
	}

	return c / (A_mixed * 2.f);
}
vecf3 GetMeanCurvatureOperator_me(Vertex* vi) 
{
	//起始半边
	auto* start_he = vi->HalfEdge();
	auto* he = start_he;

	vecf3 K(0.f, 0.f, 0.f);
	float area = 0.f;


	do {
		auto* he1 = he->Pair()->Next();
		auto* he2 = he->Pre()->Pair();

		//assert(he1->Origin() == he->Origin());
		//assert(he2->Origin() == he->Origin());

		//auto* vi = he->Origin();
		auto* vj = he->End();
		auto* v_alpha = he1->End();
		auto* v_beta = he2->End();


		
		auto cot_alpha = (v_alpha->position - vi->position).cot_theta(v_alpha->position - vi->position);
		auto cot_beta = (v_beta->position - vi->position).cot_theta(v_beta->position - vi->position);
		
		area += (float)1 / 8 * (cot_alpha + cot_beta) * (vi->position - vj->position).norm2();

		K += (cot_alpha + cot_beta) * (vi->position - vj->position);
		he = he1;
	} while (he != start_he);

	K = (float)1 / (2 * area) * K;
	return K;
}

void HeMesh2Mesh(DenoiseData* data)
{
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
}
void Mesh2HeMesh(DenoiseData* data)
{
	data->heMesh->Clear();
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
}
void DenoiseSystem::OnUpdate(Ubpa::UECS::Schedule& schedule) {
	schedule.RegisterCommand([](Ubpa::UECS::World* w) {
		auto data = w->entityMngr.GetSingleton<DenoiseData>();
		if (!data)
			return;

		if (ImGui::Begin("Denoise")) {
	
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
				
				Mesh2HeMesh(data);
				//v1为尖点
				/*auto Cot = [](const pointf3& v1, const pointf3& v2, const pointf3& v3)->float {
					vecf3 a = v2 - v1;
					vecf3 b = v3 - v1;
					return vecf3::cot_theta(a, b);
				};*/

				for (auto* p : data->heMesh->Vertices()) {
					p->position = 5 * p->position;
					if (p->IsOnBoundary()) {
						//auto he = p->HalfEdge();
						p->is_boundary = true;
					}
				}

				auto vertices = data->heMesh->Vertices();	//vertices的元素指向的是底层hemesh的顶点
				std::vector<Ubpa::pointf3> new_positions(vertices.size());
				
				for (int k = 0; k < data->k; ++k) {
					for (int j = 0; j < vertices.size(); ++j) {
						if (vertices[j]->is_boundary) {
							new_positions[j] = vertices[j]->position;
							continue;
						}

						vecf3 He = GetMeanCurvatureOperator(vertices[j]) / (2.f);
						new_positions[j] = vertices[j]->position + data->lambda * He;
					}

					for (int j = 0; j < vertices.size(); ++j) {
						vertices[j]->position = new_positions[j];
					}
				}


				HeMesh2Mesh(data);
			}
		}
		ImGui::End();
	});
}
