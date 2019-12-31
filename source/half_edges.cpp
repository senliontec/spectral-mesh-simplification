#include "half_edges.h"
#include <execution>
#include <algorithm>
#include "debugbreak.h"

half_edge_connectivity::half_edge_connectivity(size_t num_vertices, const std::vector<uint32_t>& indices)
{
	reset(num_vertices, indices);
}

void half_edge_connectivity::link_edges(std::vector<std::pair<uint32_t, uint32_t>>& vh)
{
	// Link half edges by setting the 'opposite' field
	auto next = [&](uint32_t h) { return half_edges[h].next; };
	auto vertex = [&](uint32_t h) { return half_edges[h].vertex; };
	auto sibling = [&](uint32_t h, uint32_t u) { return vertex(h) == u ? vertex(next(next(h))) : vertex(h); };
	std::sort(std::execution::par_unseq, vh.begin(), vh.end(), [&](const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b)
	{
		if(a.first != b.first) return a.first < b.first;
		return sibling(a.second, a.first) < sibling(b.second, a.first);
	});
	for(size_t t = 0; t < vh.size();)
	{
		const uint32_t u = vh[t].first;
		const uint32_t v = sibling(vh[t].second, u);
		size_t k = t;
		while(k < vh.size() && vh[k].first == u && sibling(vh[k].second, u) == v)
			k++;
		// k - t = number of half-edges for edge (u, v), equal to the number of triangles containing that edge
		// 1 : nothing to do, this is a boundary so no opposite
		// 2 : classical manifold edge
		// 3+: non-manifold edge, treat if as some kind of boundary
		if(k - t > 2) printf("Non-manifold edge detected !\n"); /// TODO: better error reporting
		if(k - t == 2)
		{
			const uint32_t h1 = vh[t].second, h2 = vh[t + 1].second;
			half_edges[h1].opposite = h2;
			half_edges[h2].opposite = h1;
		}
		t = k;
	}
}

void half_edge_connectivity::link_vertices(size_t num_vertices)
{
	// Vertex -> half-edge; always start by a boundary edge
	vertex_to_half_edge.resize(num_vertices, uint32_t(-1));
	for(uint32_t t = 0; t < (uint32_t)half_edges.size(); t++)
	{
		const uint32_t v = handle(t).next().next().vertex();
		if(vertex_to_half_edge[v] == uint32_t(-1)) vertex_to_half_edge[v] = t;
	}
	for(uint32_t v = 0; v < vertex_to_half_edge.size(); v++)
	{
		uint32_t& v2h = vertex_to_half_edge[v];
		const half_edge h = handle(v2h);
		half_edge k = h.opposite().next();
		while(k.is_valid() && k != h)
		{
			v2h = k.index;
			k = k.opposite().next();
		}
	}
}

bool half_edge_connectivity::is_boundary_vertex(uint32_t v) const
{
	const uint32_t h = vertex_to_half_edge[v];
	if(h == uint32_t(-1)) return true;
	return is_boundary_half_edge(h);
}

bool half_edge_connectivity::is_boundary_half_edge(uint32_t h) const
{
	return half_edges[h].opposite == uint32_t(-1);
}

bool half_edge_connectivity::is_connected(uint32_t u, uint32_t v) const
{
	for(uint32_t w : vertex_ring(u))
		if(w == v) return true;
	return false;
}

std::vector<uint32_t> half_edge_connectivity::joint_neighbours(uint32_t u, uint32_t v) const
{
	std::vector<uint32_t> n = vertex_ring(u).to_vector();
	size_t found = 0;
	for(uint32_t w : vertex_ring(v))
	{
		auto it = std::find(n.begin() + found, n.end(), w);
		if(it != n.end()) std::swap(n[found++], *it);
	}
	n.resize(found);
	return n;
}

std::array<uint32_t, 3> half_edge_connectivity::triangle(uint32_t h) const
{
	half_edge hh = handle(h);
	return { hh.vertex(), hh.next().vertex(), hh.next().next().vertex() };
}

std::vector<uint32_t> half_edge_connectivity::create_indices() const
{
	std::vector<uint32_t> indices;
	indices.reserve(half_edges.size());
	on_triangles([&](const std::array<uint32_t, 3>& v)
	{
		for(unsigned int k = 0; k < 3; k++)
			indices.push_back(v[k]);
	});
	return indices;
}

uint32_t half_edge_connectivity::vertex_arity(uint32_t v) const
{
	uint32_t count = 0;
	half_edge he = handle(vertex_to_half_edge[v]);
	uint32_t start = uint32_t(-1);
	while(he.is_valid() && he.index != start)
	{
		if(start == uint32_t(-1)) start = he.index;
		count++;
		he = he.next().next().opposite();
	}
	return count;
}

std::vector<uint32_t> half_edge_connectivity::boundary_edges() const
{
	std::vector<uint32_t> idx;
	for(uint32_t t = 0; t < (uint32_t)half_edges.size(); t++)
		if(half_edges[t].is_valid() && is_boundary_half_edge(t)) idx.push_back(t);
	return idx;
}

std::vector<uint32_t> half_edge_connectivity::interior_edges() const
{
	std::vector<bool> visited(half_edges.size(), false);
	std::vector<uint32_t> idx;
	for(uint32_t t = 0; t < (uint32_t)half_edges.size(); t++)
	{
		if(visited[t]) continue;
		visited[t] = true;
		if(half_edges[t].is_valid() && !is_boundary_half_edge(t))
		{
			idx.push_back(t);
			visited[half_edges[t].opposite] = true;
		}
	}
	return idx;
}

uint32_t half_edge_connectivity::find_edge(uint32_t from, uint32_t to) const
{
	half_edge he = handle(vertex_to_half_edge[from]);
	if(he.vertex() == to) return he.index;
	const uint32_t start = he.vertex();
	do he = he.next().next().opposite();
	while(he.is_valid() && he.vertex() != to && he.vertex() != start);
	return he.vertex() == to ? he.index : uint32_t(-1);
}

bool half_edge_connectivity::valid_split_edge(uint32_t h) const
{
	// Ensure we split an existing edge
	return handle(h).is_valid();
}

bool half_edge_connectivity::valid_split_edge(uint32_t h, uint32_t new_vertex) const
{
	// Vertices of the edge to split must not be connected
	if(is_connected(new_vertex)) return false;

	// Ensure we split an existing edge
	return handle(h).is_valid();
}

bool half_edge_connectivity::valid_flip_edge(uint32_t h) const
{
	// Must be a valid half edge
	half_edge he = handle(h);
	if(!he.is_valid()) return false;

	// An edge must be in 2 triangles to be flipable
	if(!he.opposite().is_valid()) return false;

	// It is impossible to flip an edge if the flipped edge already exists
	return !is_connected(he.next().vertex(), he.opposite().next().vertex());
}

bool half_edge_connectivity::valid_collapse_edge(uint32_t h) const
{
	// Vertices to merge must be connected
	half_edge he = handle(h);
	if(!he.is_valid()) return false;

	// All joint neighbors of the edge are in a face with it ("Mesh Optimization", Hoppe et al 1993)
	const uint32_t u = he.next().next().vertex(), v = he.vertex();
	const size_t n = joint_neighbours(u, v).size();
	if(n != he.arity()) return false;

	// Edge must be boundary if both vertices are boundary ("Mesh Optimization", Hoppe et al 1993)
	return !is_boundary_vertex(u) || !is_boundary_vertex(v) || he.is_boundary();
}

void half_edge_connectivity::split_edge(uint32_t he, uint32_t new_vertex)
{
	if(is_boundary_half_edge(he)) split_edge_1(he, new_vertex);
	else split_edge_2(he, new_vertex);
}

void half_edge_connectivity::split_edge_1(uint32_t he, uint32_t new_vertex)
{
	const std::array<uint32_t, 3> h = triangle_half_edges(he);
	const uint32_t j[3] = { alloc_half_edge(), alloc_half_edge(), alloc_half_edge() };
	const uint32_t v[3] = { half_edges[h[0]].vertex, half_edges[h[1]].vertex, half_edges[h[2]].vertex };
	half_edge_data& j0 = half_edges[j[0]];
	j0.next = h[0];
	j0.opposite = j[1];
	j0.vertex = new_vertex;
	half_edge_data& j1 = half_edges[j[1]];
	j1.next = h[2];
	j1.opposite = j[0];
	j1.vertex = v[1];
	half_edge_data& j2 = half_edges[j[2]];
	j2.next = j[1];
	j2.opposite = uint32_t(-1);
	j2.vertex = new_vertex;
	half_edges[h[1]].next = j[0];
	half_edges[h[2]].next = j[2];
	while(new_vertex > vertex_to_half_edge.size())
		vertex_to_half_edge.push_back(uint32_t(-1));
	vertex_to_half_edge[new_vertex] = h[0];
	vertex_to_half_edge[v[2]] = j[2];
}

void half_edge_connectivity::split_edge_2(uint32_t he, uint32_t new_vertex)
{
	// Inner edges going to/from the new vertex
	const uint32_t e[4][2] =
	{
		{ alloc_half_edge(), alloc_half_edge() },
		{ alloc_half_edge(), alloc_half_edge() },
		{ alloc_half_edge(), alloc_half_edge() },
		{ he, half_edges[he].opposite }
	};

	// Half-edges on the side
	uint32_t p[4] = { uint32_t(-1), half_edges[e[3][0]].next, uint32_t(-1), half_edges[e[3][1]].next };
	p[0] = half_edges[p[1]].next;
	p[2] = half_edges[p[3]].next;

	// Link edges
	for(unsigned t = 0; t < 4; t++)
	{
		half_edges[p[t]].next = e[(t + 3) % 4][0];
		half_edge_data& e0 = half_edges[e[t][0]];
		e0.vertex = new_vertex;
		e0.next = e[(t + 1) % 4][1];
		e0.opposite = e[t][1];
		half_edge_data& e1 = half_edges[e[t][1]];
		e1.vertex = half_edges[p[(t + 1) % 4]].vertex;
		e1.next = p[t];
		e1.opposite = e[t][0];
	}

	// Link to vertices
	while(new_vertex >= vertex_to_half_edge.size())
		vertex_to_half_edge.push_back(uint32_t(-1));
	vertex_to_half_edge[new_vertex] = e[0][1];
	const uint32_t v = half_edges[p[2]].vertex;
	if(vertex_to_half_edge[v] == e[3][1]) vertex_to_half_edge[v] = e[1][0];
}

void half_edge_connectivity::flip_edge(uint32_t he)
{
	//	h1: u -> v
	//
	//	      v
	//	      *
	//	k3  / | \  h2
	//	  /   |   \
	//	*   k1|h1   *
	//	  \   |   /
	//	k2  \ | /  h3
	//	      *
	//		  u
	half_edge_data& h1 = half_edges[he];
	half_edge_data& h2 = half_edges[h1.next];
	half_edge_data& h3 = half_edges[h2.next];
	half_edge_data& k1 = half_edges[h1.opposite];
	half_edge_data& k2 = half_edges[k1.next];
	half_edge_data& k3 = half_edges[k2.next];
	if(vertex_to_half_edge[k1.vertex] == k1.opposite) vertex_to_half_edge[k1.vertex] = h3.opposite;
	if(vertex_to_half_edge[h1.vertex] == h1.opposite) vertex_to_half_edge[h1.vertex] = k3.opposite;
	h1.vertex = h2.vertex;
	k1.vertex = k2.vertex;
	k3.next = h1.next;
	h1.next = h2.next;
	h3.next = k1.next;
	k1.next = k2.next;
	h2.next = h1.opposite;
	k2.next = he;
}

std::array<uint32_t, 6> half_edge_connectivity::collapse_edge(uint32_t he)
{
	// Change vertex of all half-edges from 'to_remove' to 'to_keep'
	half_edge h = handle(he);
	const uint32_t to_keep = h.vertex(), to_remove = h.next().next().vertex();
	h = handle(vertex_to_half_edge[to_remove]);
	const uint32_t loop_h = h.index;
	do
	{
		h = h.next().next();
		half_edges[h.index].vertex = to_keep;
		h = h.opposite();
	} while(h.is_valid() && h.index != loop_h);
	h = handle(he);

	// Edges to delete, adjust the vertex correspondence
	const std::array<uint32_t, 6> to_delete =
	{
		h.index,
		h.opposite().index,
		h.next().next().index,
		h.next().next().opposite().index,
		h.opposite().next().index,
		h.opposite().next().opposite().index
	};
	auto adjust_vertex_pointed_by = [&](half_edge k)
	{
		if(!k.is_valid()) return;
		const uint32_t v = k.vertex();
		k = k.next();
		if(vertex_to_half_edge[v] == k.index) vertex_to_half_edge[v] = k.next().next().opposite().index;
	};
	half_edge i = h.next().next().opposite();
	half_edge j = h.opposite().next().opposite();
	adjust_vertex_pointed_by(h.next());
	if(!i.is_valid()) adjust_vertex_pointed_by(h);
	adjust_vertex_pointed_by(h.opposite().next().next());
	if(!j.is_valid()) adjust_vertex_pointed_by(h.opposite().next());
	adjust_vertex_pointed_by(h.opposite().next().opposite().next().next());
	vertex_to_half_edge[to_remove] = uint32_t(-1);

	// Adjust the 4 triangles adjacent to 'he'
	if(i.is_valid())
	{
		const uint32_t tri[3] = { i.next().index, i.next().next().index, i.opposite().next().next().index };
		half_edges[tri[1]].next = tri[2];
		half_edges[tri[2]].next = tri[0];
	}
	else
	{
		i = h.next();
		half_edges[i.opposite().index].opposite = uint32_t(-1);
		dealloc_half_edge(i.index);
	}
	if(j.is_valid())
	{
		const uint32_t tri[3] = { j.next().index, j.next().next().index, j.opposite().next().index };
		half_edges[tri[1]].next = tri[2];
		half_edges[tri[2]].next = tri[0];
	}
	else
	{
		j = h.opposite().next().next();
		if(j.is_valid())
		{
			half_edges[j.opposite().index].opposite = uint32_t(-1);
			dealloc_half_edge(j.index);
		}
	}

	// Deallocate edges
	for(uint32_t t = 0; t < 6; t++)
		dealloc_half_edge(to_delete[t]);

	// A collapse can change the topology of 'to_keep' from a closed disk to an open disk, hence the following
	h = handle(vertex_to_half_edge[to_keep]);
	const uint32_t start_h = h.index;
	while(h.opposite().next().is_valid())
	{
		h = h.opposite().next();
		if(h.index == start_h) break;
	}
	if(h.index != start_h) vertex_to_half_edge[to_keep] = h.index;
	return to_delete;
}

uint32_t half_edge_connectivity::alloc_half_edge()
{
	if(free_half_edges.empty())
	{
		half_edges.emplace_back();
		return uint32_t(half_edges.size()) - 1;
	}

	std::pop_heap(free_half_edges.begin(), free_half_edges.end(), std::greater<uint32_t>());
	const uint32_t h = free_half_edges.back();
	free_half_edges.pop_back();
	return h;
}

void half_edge_connectivity::dealloc_half_edge(uint32_t h)
{
	if(h == uint32_t(-1)) return;
	half_edges[h] = half_edge_data();
	free_half_edges.push_back(h);
	std::push_heap(free_half_edges.begin(), free_half_edges.end(), std::greater<uint32_t>());
}

void half_edge_connectivity::realign()
{
	enum class flags : uint8_t
	{
		uncharted,
		added,
		added_reorient,
		visited
	};
	std::vector<flags> state(half_edges.size(), flags::uncharted);
	std::vector<uint32_t> to_visit;
	auto push = [&](half_edge he, bool reorient)
	{
		if(!he.is_valid()) return;
		if(state[he.index] == flags::uncharted)
		{
			state[he.index] = reorient ? flags::added_reorient : flags::added;
			to_visit.push_back(he.index);
		}
		else if(state[he.index] == flags::visited && reorient)
		{
			fprintf(stderr, "Non-manifold mesh: non-orientable (cannot reorient half-edge %u)\n", he.index);
			debug_break();
		}
	};

	// Add all edges, but leave them in "uncharted". It allows to process disconnected components; and
	// if the edge is pushed, then it will be in "visited" preventing it to be processed twice.
	to_visit.reserve(half_edges.size());
	for(size_t t = 0; t < half_edges.size(); t++)
		to_visit.push_back((uint32_t)t);

	while(!to_visit.empty())
	{
		const uint32_t he = to_visit.back();
		to_visit.pop_back();

		if(state[he] == flags::visited) continue;
		const std::array<uint32_t, 3> h = triangle_half_edges(he);
		if(state[he] == flags::added_reorient)
		{
			const uint32_t v[3] = { half_edges[h[0]].vertex, half_edges[h[1]].vertex, half_edges[h[2]].vertex };
			for(unsigned int t = 0; t < 3; t++)
			{
				half_edges[h[t]].next = h[(t + 2) % 3];
				half_edges[h[t]].vertex = v[(t + 2) % 3];
			}
		}
		for(unsigned int t = 0; t < 3; t++)
		{
			state[h[t]] = flags::visited;
			half_edge hh = handle(h[t]);
			push(hh.opposite(), hh.vertex() == hh.opposite().vertex());
		}
	}
}

void half_edge_connectivity::check_invariants()
{
	for(size_t t = 0; t < vertex_to_half_edge.size(); t++)
	{
		const uint32_t h = vertex_to_half_edge[t];
		if(h == uint32_t(-1)) continue;
		if(!half_edges[h].is_valid())
		{
			fprintf(stderr, "Invalid vertex_to_half_edge[%zu] !\n", t);
			debug_break();
		}
		if(uint32_t(t) != handle(h).next().next().vertex())
		{
			fprintf(stderr, "Invalid vertex triangle (for vertex %zu)\n", t);
			debug_break();
		}
		half_edge hh = handle(h);
		if(hh.is_valid())
		{
			while(hh.opposite().next().is_valid())
			{
				hh = hh.opposite().next();
				if(hh.index == h) break;
			}
			if(hh.index != h)
			{
				fprintf(stderr, "Invalid starting half-edge (for vertex %zu), should be %u\n", t, hh.index);
				debug_break();
			}
		}
	}

	std::vector<uint8_t> num_boundaries(vertex_to_half_edge.size(), 0);
	for(size_t h = 0; h < half_edges.size(); h++)
	{
		if(!half_edges[h].is_valid()) continue;
		const uint32_t v = half_edges[h].vertex;
		if(vertex_to_half_edge[v] == uint32_t(-1))
		{
			fprintf(stderr, "Invalid reference (v2h invalid)\n");
			debug_break();
		}
		const uint32_t n = half_edges[h].next;
		const uint32_t nn = half_edges[n].next;
		const uint32_t nnn = half_edges[nn].next;
		const bool correct_next = n != h && nn != h && nn != n && nnn == h;
		if(!correct_next)
		{
			fprintf(stderr, "Invalid triangle, half-edges %u -> %u -> %u -> %u\n", uint32_t(h), n, nn, nnn);
			debug_break();
		}
		const uint32_t vv = half_edges[nn].vertex;
		if(half_edges[h].opposite == uint32_t(-1) && ++num_boundaries[vv] > 2)
		{
			fprintf(stderr, "Invalid vertex %u, too many boundary edges\n", vv);
			debug_break();
		}
		const uint32_t o = half_edges[h].opposite;
		const uint32_t oo = o != uint32_t(-1) ? half_edges[o].opposite : uint32_t(-1);
		if(o != uint32_t(-1) && oo != h)
		{
			fprintf(stderr, "Invalid opposite, %u -> %u -> %u\n", uint32_t(h), o, oo);
			debug_break();
		}
		if(o != uint32_t(-1) && half_edges[o].vertex == v)
		{
			fprintf(stderr, "Invalid orientation, %zu and its opposite have the same vertex (%u)\n", h, v);
			debug_break();
		}
	}
}
