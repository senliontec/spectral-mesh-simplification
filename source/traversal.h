#pragma once

#include "half_edges.h"
#include <unordered_set>

struct k_ring
{
	template<typename F>
	k_ring(const half_edge_connectivity& connec, unsigned int k, uint32_t v, F f) :
		connectivity(connec)
	{
		traverse_k_ring(k, v, f);
	}

	template<typename F>
	k_ring& extend(unsigned int k, uint32_t v, F f)
	{
		visited_vertices.clear();
		traverse_k_ring(k, v, f);
		return *this;
	}

private:
	template<typename F>
	void traverse_k_ring(unsigned int k, uint32_t v, F f)
	{
		if(k == 0) return;
		visited_vertices.insert(v);
		traverse_1_ring(v, [&](const half_edge& h)
		{
			if(visited_edges.count(h.index) == 0) f(h);
			visited_edges.insert(h.index);

			const uint32_t vert = h.vertex();
			if(k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring(k - 1, vert, f);
		});
	}

	template<typename F>
	void traverse_1_ring(uint32_t v, F f)
	{
		half_edge h = connectivity.handle(connectivity.vertex_half_edge(v));
		uint32_t start = uint32_t(-1);
		while(h.index != start && h.is_valid())
		{
			if(start == uint32_t(-1)) start = h.index;
			f(h);
			f(h = h.next());
			f(h = h.next());
			h = h.opposite();
		}
	}

private:
	std::unordered_set<uint32_t> visited_edges;
	std::unordered_set<uint32_t> visited_vertices;
	const half_edge_connectivity& connectivity;
};
