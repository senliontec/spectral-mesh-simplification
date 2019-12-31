#pragma once

#include <array>
#include <vector>

class half_edge_connectivity;

struct half_edge_data
{
	uint32_t next;
	uint32_t opposite;
	uint32_t vertex;

	half_edge_data(uint32_t v = uint32_t(-1), uint32_t n = uint32_t(-1), uint32_t o = uint32_t(-1)) : next(n), opposite(o), vertex(v) {}
	bool is_valid() const { return vertex != uint32_t(-1) && next != uint32_t(-1); }
};

struct half_edge /// TODO: allow const AND non-const usage
{
	const half_edge_connectivity* connec;
	uint32_t index;

	bool is_valid() const;
	uint32_t vertex() const;
	half_edge next() const;
	half_edge opposite() const;
	size_t arity() const;
	bool is_boundary() const;

	bool operator==(const half_edge& he) const;
	bool operator!=(const half_edge& he) const;
};

/// TODO: structure of array
/// TODO: put opposites adjacent in the arrays (so that (h, h_opp) are always at (i, i+1))
/// TODO: face indices
class half_edge_connectivity
{
private:
	struct range_vertex_ring
	{
		struct sentinel {};
		struct iterator
		{
			half_edge he;
			uint32_t start, current; // 'start' needed to break the cycle

			iterator(half_edge h) : he(h), start(h.vertex()), current(h.vertex()) {}
			uint32_t operator*() const { return current; }
			bool operator!=(const sentinel&) const { return current != uint32_t(-1); }
			iterator& operator++()
			{
				he = he.next();
				current = he.vertex();
				he = he.next().opposite();
				if(current == start) current = uint32_t(-1);
				return *this;
			}
		};

		half_edge start;

		iterator begin() const { return start; }
		sentinel end() const { return {}; };
		std::vector<uint32_t> to_vector() const
		{
			std::vector<uint32_t> vec;
			for(uint32_t v : *this)
				vec.push_back(v);
			return vec;
		}
	};

public:
	/// Convert a list of indices to half-edge connectivity
	half_edge_connectivity(size_t num_vertices, const std::vector<uint32_t>& indices);
	template<typename It>
	half_edge_connectivity(size_t num_vertices, It begin, It end);

	/// Remove everything; memory is not freed
	void clear();

	/// Resets the half-edges to represent the specified mesh
	// Note: It is expected that \p num_vertices >= max(indices)
	void reset(size_t num_vertices, const std::vector<uint32_t>& indices);
	template<typename It>
	void reset(size_t num_vertices, It begin, It end);

	/// Return the number of half-edges
	uint32_t num_half_edges() const;

	/// Return a handle to the specified edge
	half_edge handle(uint32_t h) const;

	/// Return the vertices of the specified half-edge
	std::pair<uint32_t, uint32_t> edge_vertices(uint32_t h) const;

	/// Return whether the specified vertex has at least 1 boundary edge
	bool is_boundary_vertex(uint32_t v) const;

	/// Return whether the specified half edge is boundary
	bool is_boundary_half_edge(uint32_t h) const;

	/// Return whether \p u and \p v are connected
	bool is_connected(uint32_t u, uint32_t v) const;

	/// Return whether \p v is connected to another vertex
	bool is_connected(uint32_t v) const;

	/// Return joint neighbours of \p u and \p v
	std::vector<uint32_t> joint_neighbours(uint32_t u, uint32_t v) const;

	/// Return the triangle containing the half-edge
	std::array<uint32_t, 3> triangle(uint32_t h) const;

	/// Return the half-edges of the triangle containing \p h
	std::array<uint32_t, 3> triangle_half_edges(uint32_t h) const;

	/// Create the array of indices for this mesh
	std::vector<uint32_t> create_indices() const;

	/// Return the range of all vertices connected to \p v; any modification to the half edges of \p v will invalidate the range
	range_vertex_ring vertex_ring(uint32_t v) const;

	/// Return the first half-edge associated with the specified vertex
	uint32_t vertex_half_edge(uint32_t v) const;

	/// Return the number of vertices connected to \p v
	uint32_t vertex_arity(uint32_t v) const;

	/// Return the indices of boundary half edges
	std::vector<uint32_t> boundary_edges() const;

	/// Return the indices of interior half edges
	std::vector<uint32_t> interior_edges() const;

	/// Return the index of the half edge (\p from, \p to)
	uint32_t find_edge(uint32_t from, uint32_t to) const;

	/// Return whether the specified half edge is valid
	bool valid_half_edge(uint32_t he) const;

	/// Return whether splitting \p h with \p new_vertex is valid
	bool valid_split_edge(uint32_t h, uint32_t new_vertex) const;

	/// Return whether splitting \p h is valid, without consideration for the new vertex to use
	bool valid_split_edge(uint32_t h) const;

	/// Return  whether flipping \p h is valid
	bool valid_flip_edge(uint32_t h) const;

	/// Return whether collapsing \p h is valid; note that there must be at least 5 vertices for it to be valid
	// The removed vertices will be the one at the start of \p h (i.e., h->next->next->vertex)
	bool valid_collapse_edge(uint32_t h) const;

	/// Split the half edge; no verification made, undefined behaviour if !valid_split_edge
	void split_edge(uint32_t he, uint32_t new_vertex);

	/// Flip the half edge \p he; no verification made, undefined behaviour if !valid_flip_edge
	void flip_edge(uint32_t he);

	/// Collapse the half edge \p he; no verification made, undefined behaviour if !valid_collapse_edge
	// Return the index of removed half edges (which will all be invalid following this operation).
	std::array<uint32_t, 6> collapse_edge(uint32_t he);

	/// Call \p f on all triangles (note that they are in the right orientation)
	// Note: f shall take a std::array<uint32_t, 3> as parameter (w/o const&)
	template<typename F>
	void on_triangles(F f) const;

private:
	// Add half edge loops for each triangle
	template<typename It>
	std::vector<std::pair<uint32_t, uint32_t>> create_loops(It begin, It end);
	// Link half edges by setting the 'opposite' field
	void link_edges(std::vector<std::pair<uint32_t, uint32_t>>& vh);
	void link_vertices(size_t num_vertices);

	uint32_t alloc_half_edge();
	void dealloc_half_edge(uint32_t h);
	void realign();
	void split_edge_1(uint32_t he, uint32_t new_vertex);
	void split_edge_2(uint32_t he, uint32_t new_vertex);
	void check_invariants();

private:
	friend struct half_edge;
	std::vector<half_edge_data> half_edges;

	// Invariant: this always point to a boundary edge going *from* the vertex
	std::vector<uint32_t> vertex_to_half_edge;

	// Store invalid half edges in a heap; this allow to find them fast for re-use
	std::vector<uint32_t> free_half_edges;
};

// Half edge
inline bool half_edge::is_valid() const      { return index != uint32_t(-1) && connec->half_edges[index].is_valid(); }
inline uint32_t half_edge::vertex() const    { return index != uint32_t(-1) ? connec->half_edges[index].vertex : uint32_t(-1); }
inline half_edge half_edge::next() const     { return index != uint32_t(-1) ? half_edge { connec, connec->half_edges[index].next } : *this; }
inline half_edge half_edge::opposite() const { return index != uint32_t(-1) ? half_edge { connec, connec->half_edges[index].opposite } : *this; }
inline size_t half_edge::arity() const       { return index != uint32_t(-1) ? (!is_boundary() ? 2 : 1) : 0; }
inline bool half_edge::is_boundary() const   { return index != uint32_t(-1) && !opposite().is_valid(); }
inline bool half_edge::operator==(const half_edge& he) const { return connec == he.connec && index == he.index; }
inline bool half_edge::operator!=(const half_edge& he) const { return !operator==(he); }

// Half edge connectivity
template<typename It>
inline half_edge_connectivity::half_edge_connectivity(size_t num_vertices, It begin, It end) { reset(num_vertices, begin, end); }
inline uint32_t half_edge_connectivity::num_half_edges() const { return (uint32_t)half_edges.size(); }
inline half_edge half_edge_connectivity::handle(uint32_t h) const { return { this, h }; }
inline half_edge_connectivity::range_vertex_ring half_edge_connectivity::vertex_ring(uint32_t v) const { return { handle(vertex_to_half_edge[v]) }; }
inline uint32_t half_edge_connectivity::vertex_half_edge(uint32_t v) const { return vertex_to_half_edge[v]; }
inline bool half_edge_connectivity::valid_half_edge(uint32_t he) const { return half_edges[he].is_valid(); }
inline std::array<uint32_t, 3> half_edge_connectivity::triangle_half_edges(uint32_t h) const
{
	const uint32_t k = half_edges[h].next;
	return { h, k, half_edges[k].next };
}
inline bool half_edge_connectivity::is_connected(uint32_t v) const { return v < vertex_to_half_edge.size() && vertex_to_half_edge[v] != uint32_t(-1); }
inline std::pair<uint32_t, uint32_t> half_edge_connectivity::edge_vertices(uint32_t he) const
{
	half_edge h = handle(he);
	return { h.next().next().vertex(), h.vertex() };
}
template<typename F>
inline void half_edge_connectivity::on_triangles(F f) const
{
	for(uint32_t h = 0; h < (uint32_t)half_edges.size(); h++)
	{
		const std::array<uint32_t, 3> v = triangle(h);
		if(v[2] <= v[1] || v[2] <= v[0]) continue; // Avoids duplicates (bc 3 permutations per triangle)
		f(v);
	}
}
inline void half_edge_connectivity::clear()
{
	half_edges.clear();
	vertex_to_half_edge.clear();
	free_half_edges.clear();
}
template<typename It>
inline std::vector<std::pair<uint32_t, uint32_t>> half_edge_connectivity::create_loops(It begin, It end)
{
	std::vector<std::pair<uint32_t, uint32_t>> vh;
	while(begin != end)
	{
		const uint32_t i0 = *begin++;
		const uint32_t i1 = *begin++;
		const uint32_t i2 = *begin++;
		const uint32_t h = (uint32_t)half_edges.size();
		half_edges.emplace_back(i1, h + 1);
		half_edges.emplace_back(i2, h + 2);
		half_edges.emplace_back(i0, h + 0);
		vh.emplace_back(i0, h + 0);
		vh.emplace_back(i0, h + 2);
		vh.emplace_back(i1, h + 0);
		vh.emplace_back(i1, h + 1);
		vh.emplace_back(i2, h + 1);
		vh.emplace_back(i2, h + 2);
	}
	return vh;
}
inline void half_edge_connectivity::reset(size_t num_vertices, const std::vector<uint32_t>& indices) { reset(num_vertices, indices.begin(), indices.end()); }
template<typename It>
inline void half_edge_connectivity::reset(size_t num_vertices, It begin, It end)
{
	clear();
	std::vector<std::pair<uint32_t, uint32_t>> vh = create_loops(begin, end);
	link_edges(vh);
	realign();
	link_vertices(num_vertices);
}

/// Flip selected edges, and return the number of flipped edges
template<typename F>
inline unsigned int flip_edges(half_edge_connectivity& connec, F select_edge = [](uint32_t){ return true; })
{
	unsigned int num_flips = 0;
	for(uint32_t he : connec.interior_edges())
	{
		if(!connec.valid_flip_edge(he) || !select_edge(he)) continue;
		num_flips++;
		connec.flip_edge(he);
	}
	return num_flips;
}

/// Flip selected edges only if it brings the vertex arity toward the specified objective function
template<typename F, typename G>
inline unsigned int step_flip_regularize(half_edge_connectivity& connec, F target_vertex_arity, G select_edge = [](uint32_t){ return true; })
{
	return flip_edges(connec, [&](uint32_t he)
	{
		if(!select_edge(he)) return false;
		half_edge h = connec.handle(he);
		const uint32_t u = h.next().next().vertex(), v = h.vertex();
		const int tgt[2] = { (int)target_vertex_arity(u), (int)target_vertex_arity(v) };
		const int a[2] = { (int)connec.vertex_arity(u), (int)connec.vertex_arity(v) };
		const int delta_pre = std::abs(a[0] - tgt[0]) + std::abs(a[1] - tgt[1]);
		const int delta_post = std::abs(a[0] - 1 - tgt[0]) + std::abs(a[1] - 1 - tgt[1]);
		return delta_post < delta_pre;
	});
}
