#pragma once

#include <UHEMesh/HEMesh.h>

#include <UGM/UGM.h>
struct Vertex;
struct Edge;
struct Triangle;
struct HalfEdge;

using HEMeshXTraits = Ubpa::HEMeshTraits<Vertex, Edge, Triangle, HalfEdge>;

struct Vertex : Ubpa::TVertex<HEMeshXTraits> {
	// you can add any attributes and mothods to Vertex
	Ubpa::pointf3 position{ 0.f };

	//added
	bool is_boundary{ false };

	float area() {
		
	}
	
	float cot() {
		auto* he = HalfEdge();
		auto* pre = he->End();
		he = he->RotateNext();
		auto* next = he->End();
		/*Ubpa::vecf3 v1(pre->position - position);
		Ubpa::vecf3 v2(next->position - position);
		return Ubpa::vecf3::cot_theta(v1, v2);*/
		return 0;
	}

};

struct Edge : Ubpa::TEdge<HEMeshXTraits> {
	// you can add any attributes and mothods to Edge

	// [example]

	// Ubpa::pointf3 Midpoint() const {
	//     auto* p = HalfEdge()->Origin();
    //     auto* q = HalfEdge()->End();
	//     return Ubpa::pointf3::combine(std::array{ p,q }, 0.5f);
	// }
};

struct Triangle : Ubpa::TPolygon<HEMeshXTraits> {
	// you can add any attributes and mothods to Triangle

	// [example]
	// 
	//float area{ 0.f };
	
	/*float  GetArea() {

		auto lengthV3 = [](float x, float y, float z) ->float {
			return sqrt(x * x + y * y + z * z);
		};

		auto p0 = HalfEdge()->Origin()->position;
		auto p1 = HalfEdge()->Next()->Origin()->position;
		auto p2 = HalfEdge()->Next()->Next()->Origin()->position;
		auto d01 = p1 - p0;
		auto d02 = p2 - p0;
		auto normal = d01.cross(d02);

		float area = 0.5f * lengthV3(normal[0], normal[1], normal[2]);
		
		return area;
	}*/
	// 
	// bool IsTriangle() const {
	//     return Degree() == 3;
	// }
	// 
	/* void UpdateArea() {
	     assert(IsTriangle());
	     auto* p0 = HalfEdge()->Origin();
	     auto* p1 = HalfEdge()->HalfEdge()->Origin();
	     auto* p2 = HalfEdge()->HalfEdge()->HalfEdge()->Origin();
	     auto d01 = p1 - p0;
	     auto d02 = p2 - p0;
	     area = 0.5f * d02.cross(d01);
	 }*/
};

struct HalfEdge : Ubpa::THalfEdge<HEMeshXTraits> {
	// you can add any attributes and mothods to HalfEdge

};

struct HEMeshX : Ubpa::HEMesh<HEMeshXTraits> {
	// you can add any attributes and mothods to HEMeshX
};
