#pragma once

#pragma warning(push)
#pragma warning(disable: 4127)
#include <Eigen/Dense>
#pragma warning(pop)
#include "half_edges.h"

struct mesh
{
	struct triangle
	{
		uint32_t idx[3];

		triangle() = default;
		triangle(triangle&&) = default;
		triangle(const triangle&) = default;
		triangle& operator=(triangle&&) = default;
		triangle& operator=(const triangle&) = default;
		triangle(const std::array<uint32_t, 3>& t) { set(t); }

		void set(uint32_t i, uint32_t j, uint32_t k) { idx[0] = i; idx[1] = j; idx[2] = k; }
		void set(const std::array<uint32_t, 3>& t) { set(t[0], t[1], t[2]); }
		triangle& operator=(const std::array<uint32_t, 3>& t) { set(t); return *this; }
		uint32_t operator[](unsigned int k) const { return idx[k]; }
		uint32_t& operator[](unsigned int k) { return idx[k]; }
		bool contains(uint32_t v) const { return idx[0] == v || idx[1] == v || idx[2] == v; }
	};

	std::vector<Eigen::Vector3d> vertices;
	std::vector<triangle> triangles;

	mesh() = default;
	mesh(mesh&&) = default;
	mesh(const mesh&) = default;
	mesh& operator=(mesh&&) = default;
	mesh& operator=(const mesh&) = default;
	mesh(const char* filename);

	void save(const std::string& filename) const;

	Eigen::MatrixXd matrix_vertices() const;
	Eigen::MatrixXi matrix_triangles() const;

	void scale(double a)
	{
		for(Eigen::Vector3d& v : vertices)
			v *= a;
	}

	half_edge_connectivity half_edges() const
	{
		const uint32_t* begin = (const uint32_t*)triangles.data();
		return { vertices.size(), begin, begin + 3 * triangles.size() };
	}

	bool empty() const { return vertices.empty(); }
	size_t num_vertices() const { return vertices.size(); }
	size_t num_triangles() const { return triangles.size(); }
};

struct weighted_position
{
	Eigen::Vector3d position;
	double weight;

	weighted_position() = default;
	weighted_position(const Eigen::Vector3d& p, double c) : position(p), weight(c) {}
	weighted_position(const std::pair<Eigen::Vector3d, double>& p) : position(p.first), weight(p.second) {}
	weighted_position(weighted_position&&) = default;
	weighted_position(const weighted_position&) = default;
	weighted_position& operator=(weighted_position&&) = default;
	weighted_position& operator=(const weighted_position&) = default;
	weighted_position& operator=(const std::pair<Eigen::Vector3d, double>& p)
	{
		position = p.first;
		weight = p.second;
		return *this;
	}
};
