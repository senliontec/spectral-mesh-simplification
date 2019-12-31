#include "reduce.h"
#include "mesh.h"
#include "metric.h"
#include "min_priority_queue.h"
#include "traversal.h"
#include "util.h"
#include "ranges.h"
#include <execution>

void remove_standalone_vertices(mesh& obj, const half_edge_connectivity& connectivity)
{
	std::vector<uint32_t> map(obj.vertices.size(), uint32_t(-1));
	uint32_t u = 0;
	for(uint32_t v = 0; v < obj.vertices.size(); v++)
	{
		if(!connectivity.is_connected(v)) continue;
		map[v] = u;
		if(u != v) obj.vertices[u] = obj.vertices[v];
		u++;
	}
	obj.vertices.resize(u);
	for(size_t t = 0; t < obj.triangles.size(); t++)
		for(uint32_t& idx : obj.triangles[t].idx)
			idx = map[idx];
}

namespace detail
{
	weighted_position invalid_wpos()
	{
		return { Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN()), std::numeric_limits<double>::infinity() };
	}

	struct candidate_operation
	{
		weighted_position target = invalid_wpos();
		uint32_t index = uint32_t(-1);
		bool is_collapse = true;

		candidate_operation() = default;
		candidate_operation(uint32_t h, const weighted_position& p) : target(p), index(h) {}

		bool operator<(const candidate_operation& c) const { return target.weight < c.target.weight; }
		void make_cost_infinite() { target.weight = std::numeric_limits<double>::infinity(); }
	};

#if 0
	struct statistics
	{
		statistics(const size_t&) {}
		void on_operation(bool collapse) {}
		void on_invalid() {}
		void on_outdated() {}
		void on_reduction_end(size_t) {}
	};
#else
	struct statistics
	{
		const size_t& current_size;
		size_t num_total = 0, num_flips = 0;
		size_t num_invalid_collapse = 0, num_outdated_collapse = 0;
		size_t num_invalid_flip = 0, num_outdated_flip = 0;
		bool last_is_collapse = false;

		statistics(const size_t& current_size) : current_size(current_size) {}
		void on_operation(bool collapse)
		{
			last_is_collapse = collapse;
			num_total++;
			if(!collapse) num_flips++;
			print_progress<0>("(reduce) " + std::to_string(current_size));
		}
		void on_invalid() { (last_is_collapse ? num_invalid_collapse : num_invalid_flip)++; }
		void on_outdated() { (last_is_collapse ? num_outdated_collapse : num_outdated_flip)++; }
		void on_reduction_end(size_t remaining)
		{
			const size_t num_invalid = num_invalid_collapse + num_invalid_flip;
			const size_t num_outdated = num_outdated_collapse + num_outdated_flip;
			printf("Edges statistics:\n"
				"\tprocessed    : %zu\n"
				"\t  + flips    : %zu (%.2f%%)\n"
				"\t  + invalid  : %zu (%.2f%%)\n"
				"\t    + flips  : %zu (%.2f%%)\n"
				"\t  + outdated : %zu (%.2f%%)\n"
				"\t    + flips  : %zu (%.2f%%)\n"
				"\tunprocessed  : %zu\n",
				num_total,
				num_flips, 100.0 * double(num_flips) / double(num_total),
				num_invalid, 100.0 * double(num_invalid) / double(num_total),
				num_invalid_flip, 100.0 * double(num_invalid_flip) / double(num_total),
				num_outdated, 100.0 * double(num_outdated) / double(num_total),
				num_outdated_flip, 100.0 * double(num_outdated_flip) / double(num_total),
				remaining);
		}
	};
#endif
}

void reduce_stream(mesh& obj, size_t target_size, reduction_metric& metric, const reduce_options& options)
{
	auto edge_index = [&options](uint32_t h, bool collapse) { return options.flip_edges ? 2 * h + uint32_t(collapse) : h; };
	auto candidate_index = [&edge_index](const detail::candidate_operation& c) { return edge_index(c.index, c.is_collapse); };
	auto discard_candidate = [](detail::candidate_operation& c) { c.make_cost_infinite(); };

	// Setup data structures
	half_edge_connectivity connectivity = obj.half_edges();
	min_priority_queue<detail::candidate_operation, decltype(candidate_index)> candidates(candidate_index);
	metric.setup((const mesh&)obj, (const half_edge_connectivity&)connectivity);

	// Helpers
	auto add_collapse = [&](const half_edge& he)
	{
		if(!he.is_valid()) return;
		const uint32_t ho = he.opposite().index;
		if(!options.directed_edges && !he.is_boundary() && he.index < ho) return;

		const uint32_t to_keep = he.vertex();
		const uint32_t to_remove = he.next().next().vertex();

		detail::candidate_operation c;
		c.index = he.index;
		c.target = metric.cost_collapse(he.index, to_keep, to_remove);
		c.is_collapse = true;
		candidates.push(c);
	};
	auto add_flip = [&](const half_edge& he)
	{
		if(!options.flip_edges || !he.is_valid()) return;
		if(!options.directed_edges && !he.is_boundary() && candidates.contains(edge_index(he.opposite().index, false))) return;

		detail::candidate_operation c;
		c.index = he.index;
		c.target.weight = metric.cost_flip(he.index);
		c.is_collapse = false;
		candidates.push(c);
	};
	auto add_operation = [&](const half_edge& he)
	{
		add_collapse(he);
		add_flip(he);
	};
	auto update_neighbourhood = [&](uint32_t v, unsigned int rings)
	{
		if(rings == 1)
		{
			half_edge h = connectivity.handle(connectivity.vertex_half_edge(v));
			uint32_t index = uint32_t(-1);
			while(h.is_valid() && h.index != index)
			{
				if(index == uint32_t(-1)) index = h.index;
				add_operation(h);
				h = h.next().next();
				add_operation(h);
				h = h.opposite();
			}
		}
		else
			k_ring(connectivity, rings, v, add_operation); /// TODO: optimize this, it is very inefficient (especially for 1-rings)
	};
	auto perform_collapse = [&](const detail::candidate_operation& c)
	{
		const auto [to_remove, to_keep] = connectivity.edge_vertices(c.index);
		for(uint32_t h : connectivity.collapse_edge(c.index))
			candidates.update(edge_index(h, true), discard_candidate);
		obj.vertices[to_keep] = c.target.position;
		const unsigned int rings = metric.post_collapse(c.index, to_keep, to_remove);
		update_neighbourhood(to_keep, rings);
	};
	auto perform_flip = [&](const detail::candidate_operation& c)
	{
		connectivity.flip_edge(c.index);
		const auto [u, v] = connectivity.edge_vertices(c.index);
		const unsigned int rings = metric.post_flip(c.index);
		k_ring(connectivity, rings, u, add_operation).extend(rings, v, add_operation);
	};

	// Fill the candidate queue
	if(options.flip_edges)
		candidates.reserve(connectivity.num_half_edges(), connectivity.num_half_edges() * 2);
	else
		candidates.reserve(connectivity.num_half_edges() / 2, connectivity.num_half_edges());
	for(uint32_t hx = 0; hx < connectivity.num_half_edges(); hx++) /// TODO: parallel
		add_operation(connectivity.handle(hx));

	// Simplify the mesh
	size_t current_size = obj.num_vertices();
	detail::statistics stats(current_size);
	while(!candidates.empty() && current_size > target_size)
	{
		const detail::candidate_operation c = candidates.pop();
		stats.on_operation(c.is_collapse);

		// Check that the operation is still valid
		if(c.is_collapse)
		{
			if(!connectivity.valid_collapse_edge(c.index))
			{
				stats.on_invalid();
				continue;
			}
			if(!metric.collapse_still_valid(c.index))
			{
				stats.on_outdated();
				add_collapse(connectivity.handle(c.index));
				continue;
			}
		}
		else
		{
			if(!connectivity.valid_flip_edge(c.index))
			{
				stats.on_invalid();
				continue;
			}
			if(!metric.flip_still_valid(c.index))
			{
				stats.on_outdated();
				add_flip(connectivity.handle(c.index));
				continue;
			}
		}

		// Perform the operation
		if(c.is_collapse)
		{
			perform_collapse(c);
			current_size--;
		}
		else
			perform_flip(c);
	}

	// Finalize
	stats.on_reduction_end(candidates.size());
	obj.triangles.clear();
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { obj.triangles.push_back(t); });
	remove_standalone_vertices(obj, connectivity);
}

void reduce_stream_noflip(mesh& obj, size_t target_size, reduction_metric& metric, const reduce_options& options)
{
	auto candidate_index = [](const detail::candidate_operation& c) { return c.index; };
	auto discard_candidate = [](detail::candidate_operation& c) { c.make_cost_infinite(); };

	// Setup data structures
	half_edge_connectivity connectivity = obj.half_edges();
	min_priority_queue<detail::candidate_operation, decltype(candidate_index)> candidates(candidate_index);
	metric.setup((const mesh&)obj, (const half_edge_connectivity&)connectivity);

	// Helpers to parallelize the filling of the queue
	std::vector<detail::candidate_operation> batch;
	auto compute_batch = [&]()
	{
		std::for_each(std::execution::par_unseq, batch.begin(), batch.end(), [&](detail::candidate_operation& op)
		{
			const auto [to_remove, to_keep] = connectivity.edge_vertices(op.index);
			op.target = metric.cost_collapse(op.index, to_keep, to_remove);
		});
		for(detail::candidate_operation& op : batch)
			candidates.push(std::move(op));
		batch.clear();
	};

	// Helpers
	auto add_collapse = [&](const half_edge& he)
	{
		if(!he.is_valid()) return;
		const uint32_t ho = he.opposite().index;
		if(!options.directed_edges && !he.is_boundary() && he.index < ho) return;

		const uint32_t to_keep = he.vertex();
		const uint32_t to_remove = he.next().next().vertex();

		detail::candidate_operation c;
		c.index = he.index;
		c.is_collapse = true;
		batch.push_back(std::move(c));
	};
	auto add_neighbourhood = [&](uint32_t v, unsigned int rings)
	{
		if(rings == 1)
		{
			half_edge h = connectivity.handle(connectivity.vertex_half_edge(v));
			uint32_t index = uint32_t(-1);
			while(h.is_valid() && h.index != index)
			{
				if(index == uint32_t(-1)) index = h.index;
				add_collapse(h);
				h = h.next().next();
				add_collapse(h);
				h = h.opposite();
			}
		}
		else
			k_ring(connectivity, rings, v, add_collapse); /// TODO: optimize this, it is very inefficient (especially for 1-rings)
	};
	auto perform_collapse = [&](const detail::candidate_operation& c)
	{
		const auto [to_remove, to_keep] = connectivity.edge_vertices(c.index);
		for(uint32_t h : connectivity.collapse_edge(c.index))
			candidates.update(h, discard_candidate);
		obj.vertices[to_keep] = c.target.position;
		const unsigned int rings = metric.post_collapse(c.index, to_keep, to_remove);
		add_neighbourhood(to_keep, rings);
		compute_batch();
	};

	// Fill the candidate queue
	candidates.reserve(connectivity.num_half_edges() / 2, connectivity.num_half_edges());
	for(uint32_t hx = 0; hx < connectivity.num_half_edges(); hx++) /// TODO: parallel
		add_collapse(connectivity.handle(hx));
	compute_batch();

	// Simplify the mesh
	size_t current_size = obj.num_vertices();
	detail::statistics stats(current_size);
	while(!candidates.empty() && current_size > target_size)
	{
		const detail::candidate_operation c = candidates.pop();
		stats.on_operation(c.is_collapse);

		// Check that the operation is still valid
		if(!connectivity.valid_collapse_edge(c.index))
		{
			stats.on_invalid();
			continue;
		}
		if(!metric.collapse_still_valid(c.index))
		{
			stats.on_outdated();
			add_collapse(connectivity.handle(c.index));
			compute_batch();
			continue;
		}

		// Perform the operation
		perform_collapse(c);
		current_size--;
	}

	// Finalize
	stats.on_reduction_end(candidates.size());
	obj.triangles.clear();
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { obj.triangles.push_back(t); });
	remove_standalone_vertices(obj, connectivity);
}

void optimize_flips(mesh& obj, reduction_metric& metric)
{
	struct candidate_flip
	{
		double cost = std::numeric_limits<double>::quiet_NaN();
		uint32_t index = uint32_t(-1);

		candidate_flip() = default;
		candidate_flip(uint32_t h, double c) : index(h), cost(c) {}
		bool operator<(const candidate_flip& c) const { return cost < c.cost; }
	};

	auto candidate_index = [](const candidate_flip& c) { return c.index; };

	// Setup data structures
	half_edge_connectivity connectivity = obj.half_edges();
	min_priority_queue<candidate_flip, decltype(candidate_index)> candidates(candidate_index);
	metric.setup((const mesh&)obj, (const half_edge_connectivity&)connectivity);

	// Helpers
	auto add_flip = [&](const half_edge& he)
	{
		if(!he.is_valid()) return;
		if(!he.is_boundary() && candidates.contains(he.opposite().index)) return;
		candidates.push(candidate_flip(he.index, metric.cost_flip(he.index)));
	};
	auto perform_flip = [&](const candidate_flip& c)
	{
		connectivity.flip_edge(c.index);
		const auto [u, v] = connectivity.edge_vertices(c.index);
		const unsigned int rings = metric.post_flip(c.index);
		k_ring(connectivity, rings, u, add_flip).extend(rings, v, add_flip);
	};

	// Fill the candidate queue
	candidates.reserve(connectivity.num_half_edges() / 2, connectivity.num_half_edges());
	for(uint32_t hx = 0; hx < connectivity.num_half_edges(); hx++) /// TODO: parallel
		add_flip(connectivity.handle(hx));

	// Simplify the mesh
	while(!candidates.empty())
	{
		const candidate_flip c = candidates.pop();

		// Check that the operation is still valid
		if((c.cost > 0 && std::isinf(c.cost)) || std::isnan(c.cost)) break;
		if(!connectivity.valid_flip_edge(c.index)) continue;
		if(!metric.flip_still_valid(c.index))
		{
			add_flip(connectivity.handle(c.index));
			continue;
		}

		// Perform the operation
		perform_flip(c);
	}

	// Finalize
	obj.triangles.clear();
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { obj.triangles.push_back(t); });
}
