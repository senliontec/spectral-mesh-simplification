#pragma once

#include "metric.h"
#include "geometry.h"
#include "traversal.h"

struct metric_quadrics : public reduction_metric
{
	bool recompute = false;
	bool avoid_flips = true;
	bool lock_boundaries = false;

	std::vector<Eigen::Matrix4d> quadrics;
	const mesh* object;
	const half_edge_connectivity* connectivity;

	void setup(const mesh& obj, const half_edge_connectivity& connec)
	{
		object = &obj;
		connectivity = &connec;
		quadrics = geometry_quadrics(obj);
	}

	std::pair<Eigen::Vector3d, double> cost_collapse(uint32_t he, uint32_t to_keep, uint32_t to_remove) const
	{
		if(lock_boundaries)
		{
			const std::pair<uint32_t, uint32_t> v = connectivity->edge_vertices(he);
			if(connectivity->is_boundary_vertex(v.first) || connectivity->is_boundary_vertex(v.second))
				return { {}, std::numeric_limits<double>::infinity() };
		}

		const Eigen::Vector3d fallback = (object->vertices[to_keep] + object->vertices[to_remove]) / 2.0;
		const std::pair<Eigen::Vector3d, double> result = qem_merge(quadrics[to_keep], quadrics[to_remove], fallback);
		
		auto has_triangle_flip = [this, &result](uint32_t vertex, uint32_t exclude) -> bool
		{
			const Eigen::Vector3d& old_pos = object->vertices[vertex];
			const Eigen::Vector3d& new_pos = result.first;

			half_edge h = connectivity->handle(connectivity->vertex_half_edge(vertex));
			uint32_t index = uint32_t(-1);
			while(h.is_valid() && h.index != index)
			{
				if(index == uint32_t(-1)) index = h.index;
				const uint32_t vx = h.vertex();
				const Eigen::Vector3d& v = object->vertices[vx];
				h = h.next();
				const uint32_t wx = h.vertex();
				const Eigen::Vector3d& w = object->vertices[wx];
				h = h.next().opposite();

				if(vx == exclude || wx == exclude) continue;
				const Eigen::Vector3d old_normal = (v - old_pos).cross(w - old_pos);
				const Eigen::Vector3d new_normal = (v - new_pos).cross(w - new_pos);
				if(old_normal.dot(new_normal) <= 0) return true;
			}
			return false;
		};

		if(avoid_flips && (has_triangle_flip(to_keep, to_remove) || has_triangle_flip(to_remove, to_keep)))
		{
			const double nan = std::numeric_limits<double>::quiet_NaN();
			return { { nan, nan, nan }, std::numeric_limits<double>::infinity() };
		}

		return result;
	}

	unsigned int post_collapse(uint32_t, uint32_t to_keep, uint32_t to_remove)
	{
		if(recompute)
		{
			quadrics[to_keep] = vertex_quadric(*object, *connectivity, to_keep);
			for(uint32_t v : connectivity->vertex_ring(to_keep))
				quadrics[v] = vertex_quadric(*object, *connectivity, v);
			return 2;
		}
		else
		{
			quadrics[to_keep] += quadrics[to_remove];
			return 1;
		}
	}

	serialized_parameters save() const
	{
		serialized_parameters p;
		serialize(p, "recompute", recompute);
		serialize(p, "noflips", avoid_flips);
		serialize(p, "keepborders", lock_boundaries);
		return p;
	}

	void load(const serialized_parameters& p)
	{
		deserialize(p, "recompute", recompute);
		deserialize(p, "noflips", avoid_flips);
		deserialize(p, "keepborders", lock_boundaries);
	}

};
