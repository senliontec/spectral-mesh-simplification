#pragma once

#include "metric.h"

struct metric_edge_length : public reduction_metric
{
	const mesh* obj;

	void setup(const mesh& object, const half_edge_connectivity&)
	{
		obj = &object;
	}

	std::pair<Eigen::Vector3d, double> cost_collapse(uint32_t, uint32_t to_keep, uint32_t to_remove) const
	{
		const Eigen::Vector3d pos = obj->vertices[to_keep] * 0.5 + obj->vertices[to_remove] * 0.5;
		const double cost = (obj->vertices[to_keep] - obj->vertices[to_remove]).norm();	
		return { pos, cost };
	}

	unsigned int post_collapse(uint32_t, uint32_t, uint32_t)
	{
		return 1;
	}
};
