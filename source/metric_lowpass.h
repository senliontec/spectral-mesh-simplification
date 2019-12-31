#pragma once

#include "metric.h"
#pragma warning(push)
#pragma warning(disable: 4127)
#include <Eigen/Sparse>
#pragma warning(pop)

// |P Z - M^-1 L P F|^2, where Z = M^-1 L F
struct metric_lowpass : public reduction_metric
{
	unsigned int num_eigenvectors = 0;
	bool optimal_pos = true;

	Eigen::VectorXd norms;
	Eigen::MatrixXd signals;
	Eigen::MatrixXd signals_L;
	Eigen::SparseMatrix<double, Eigen::RowMajor> projection;
	std::vector<uint32_t> mapping;
	double total_cost = 0;

	std::vector<mesh::triangle> triangles;
	std::vector<Eigen::Vector3d> cotangents;
	std::vector<double> areas;
	std::vector<uint32_t> half_edge_triangle;

	using cost_diff = std::vector<std::pair<uint32_t, double>>;
	struct edge_continuation
	{
		cost_diff diff;
		double alpha = std::numeric_limits<double>::quiet_NaN();

		edge_continuation() = default;
	};
	mutable std::vector<edge_continuation> computed_collapse;
	mutable std::vector<cost_diff> computed_flip;
	std::vector<uint32_t> vertex_revision;
	mutable std::vector<uint32_t> collapse_revision;
	mutable std::vector<uint32_t> flip_revision;
	uint32_t mesh_revision = 0;

	const mesh* object = nullptr;
	const half_edge_connectivity* connectivity = nullptr;

	~metric_lowpass();

	void setup(const mesh& obj, const half_edge_connectivity& connec) override;

	bool collapse_still_valid(uint32_t half_edge) const override;
	bool flip_still_valid(uint32_t half_edge) const override;

	std::pair<Eigen::Vector3d, double> cost_collapse(uint32_t h, uint32_t to_keep, uint32_t to_remove) const override;
	unsigned int post_collapse(uint32_t h, uint32_t to_keep, uint32_t to_remove) override;

	double cost_flip(uint32_t h) const override;
	unsigned int post_flip(uint32_t h) override;

	serialized_parameters save() const override;
	void load(const serialized_parameters&) override;
};
