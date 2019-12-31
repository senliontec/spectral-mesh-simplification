#include "metric_lowpass.h"
#include "mesh.h"
#include "geometry.h"
#include "traversal.h"
#include "ranges.h"
#include "util.h"
#include <random>
#include <chrono>

void metric_lowpass::setup(const mesh& obj, const half_edge_connectivity& connec)
{
	using clock = std::chrono::high_resolution_clock;
	const clock::time_point ts = clock::now();

	object = &obj;
	connectivity = &connec;
	mapping.resize(obj.vertices.size(), uint32_t(-1));

	triangles = object->triangles;
	cotangents = triangle_cotangents(obj);
	areas = triangle_areas(obj);

	computed_collapse.resize(connec.num_half_edges());
	computed_flip.resize(connec.num_half_edges());
	vertex_revision.resize(obj.vertices.size(), mesh_revision);
	collapse_revision.resize(connec.num_half_edges(), mesh_revision);
	flip_revision.resize(connec.num_half_edges(), mesh_revision);

	half_edge_triangle.resize(connec.num_half_edges(), uint32_t(-1));
	for(uint32_t t = 0; t < (uint32_t)obj.triangles.size(); t++)
	{
		const uint32_t he = connec.find_edge(obj.triangles[t][0], obj.triangles[t][1]);
		assert(he != uint32_t(-1));
		half_edge h = connec.handle(he);
		half_edge_triangle[h.index] = t;
		half_edge_triangle[(h = h.next()).index] = t;
		half_edge_triangle[(h = h.next()).index] = t;
	}

	norms.setZero(obj.vertices.size());
	projection.resize(obj.vertices.size(), obj.vertices.size());
	projection.setIdentity();

	Eigen::VectorXd mass = mass_barycentric(obj.vertices.size(), triangles, areas);
	Eigen::SparseMatrix<double> L = laplacian_geometric(obj.vertices.size(), triangles, cotangents);

	smallest_eigenvectors(L, mass, num_eigenvectors, nullptr, &signals);
	signals.colwise().normalize();

	signals_L = mass.cwiseInverse().asDiagonal() * L * signals;

	/// Note: here, signals are pretty big and this cause bad perfs (for example, 100 eigenvectors means
	/// having PF being on average 30x100 = 3000 coefficients ~= 24Kb, *per candidate collapse*). Tested
	/// strategies for dimensionality reduction:
	///  - PCA: sadly, the variance is uniformly dispatched between eigenvectors, so for x% cumulative
	///    variance about x eigenvectors are needed (ie, 90 eigenvectors for 90% cumulative variance...)
	///    This holds for both a PCA on signals, and (signals^T | signals_L^T)^T.
	///    A PCA on (signals | signals_L) will need 50% of the columns, since half of them are linear
	///    combinations of the previous ones.

	signals.transposeInPlace();
	signals_L.transposeInPlace();

	const clock::time_point te = clock::now();
	printf("Setup time: %.6f s\n", std::chrono::duration_cast<std::chrono::milliseconds>(te - ts).count() / 1000.0);
}

bool metric_lowpass::collapse_still_valid(uint32_t h) const
{
	const uint32_t rev = collapse_revision[h];
	const auto [to_remove, to_keep] = connectivity->edge_vertices(h);
	for(uint32_t v : connectivity->vertex_ring(to_keep))
		if(vertex_revision[v] > rev) return false;
	for(uint32_t v : connectivity->vertex_ring(to_remove))
		if(vertex_revision[v] > rev) return false;
	return true;
}

bool metric_lowpass::flip_still_valid(uint32_t h) const
{
	// Vertices: [old edge, new edge]
	half_edge hx = connectivity->handle(h);
	const uint32_t vertices[4] =
	{
		hx.vertex(),
		hx.opposite().vertex(),
		hx.next().vertex(),
		hx.opposite().next().vertex()
	};

	const uint32_t rev = flip_revision[h];
	for(uint32_t v : vertices)
		if(vertex_revision[v] > rev) return false;
	return true;
}

std::pair<Eigen::Vector3d, double> metric_lowpass::cost_collapse(uint32_t h, uint32_t to_keep, uint32_t to_remove) const
{
	// Get 2-rings
	thread_local std::vector<uint32_t> V, F;
	std::unordered_map<uint32_t, uint32_t> V_inv;
	{
		V.clear();
		F.clear();
		get_edge_2_ring(*connectivity, half_edge_triangle, h, V, F);

		for(size_t t = 0; t < V.size(); t++)
			V_inv[(uint32_t)V[t]] = (uint32_t)t;
	}

	// Inner vertices
	std::vector<bool> inner(V.size(), false);
	for(uint32_t f : F)
	{
		const mesh::triangle& tri = triangles[f];
		if(tri.contains(to_keep) != tri.contains(to_remove))
			for(uint32_t v : tri.idx)
				inner[V_inv.at(v)] = true;
	}

	// Previous local cost
	double prev_cost_local = 0;
	for(size_t v = 0; v < V.size(); v++)
		if(inner[v]) prev_cost_local += norms[V[v]];

	// Precompute signals restriction
	Eigen::MatrixXd PF(signals.rows(), V.size());
	for(size_t v = 0; v < V.size(); v++)
	{
		if(V[v] == to_remove)
			PF.col(v).setConstant(std::numeric_limits<double>::quiet_NaN());
		else if(V[v] != to_keep)
			PF.col(v) = signals.col(V[v]);
	}
	const uint32_t local_vK = V_inv[to_keep];
	const uint32_t local_vR = V_inv[to_remove];

	// Preallocate structures
	Eigen::SparseMatrix<double> L(V.size(), V.size());
	Eigen::MatrixXd LPF(signals.rows(), V.size());

	auto eval = [&](double alpha, std::vector<std::pair<uint32_t, double>>& diff) -> std::pair<Eigen::Vector3d, double>
	{
		const Eigen::Vector3d pos = object->vertices[to_keep] * (1 - alpha) + object->vertices[to_remove] * alpha;
		auto new_pos = [&](uint32_t i) { return i == to_keep || i == to_remove ? pos : object->vertices[i]; };
		auto new_index = [&](uint32_t i) { return i == to_remove ? to_keep : i; };

		// Restrict signals
		PF.col(local_vK) = (signals.col(to_keep) * (1 - alpha) + signals.col(to_remove) * alpha);

		// Compute mass and laplacian
		Eigen::VectorXd M = Eigen::VectorXd::Zero(V.size());
		thread_local std::vector<Eigen::Triplet<double>> coeffs;
		coeffs.clear();
		for(uint32_t f : F)
		{
			// Exclude the 2 triangles that will get removed by the edge collapse
			const mesh::triangle& tri = triangles[f];
			const bool has_keep = tri.contains(to_keep);
			const bool has_remove = tri.contains(to_remove);
			if(has_keep && has_remove) continue;

			// Cotangents and area
			Eigen::Vector3d c = cotangents[f];
			double a = areas[f];
			if(has_keep || has_remove)
			{
				c = triangle_cotangents(new_pos(tri[0]), new_pos(tri[1]), new_pos(tri[2]));
				a = triangle_area(new_pos(tri[0]), new_pos(tri[1]), new_pos(tri[2]));
			}

			// Local indices, after collapse (=> re-index to_remove)
			const uint32_t local_tri[3] = { V_inv.at(new_index(tri[0])), V_inv.at(new_index(tri[1])), V_inv.at(new_index(tri[2])) };

			// Mass
			for(uint32_t v : local_tri)
				M[v] += a / 3;

			// Laplacian
			auto add_coeff = [&inner, &to_remove, &local_vR](uint32_t row, uint32_t col, double val)
			{
				if(col == local_vR || !inner[col]) return;
				coeffs.emplace_back(row, col, val);
			};
			for(unsigned int j = 0; j < 3; j++)
			{
				add_coeff(local_tri[j], local_tri[j], c[(j + 1) % 3] + c[(j + 2) % 3]);
				add_coeff(local_tri[j], local_tri[(j + 1) % 3], -c[(j + 2) % 3]);
				add_coeff(local_tri[j], local_tri[(j + 2) % 3], -c[(j + 1) % 3]);
			}
		}
		L.setFromTriplets(coeffs.begin(), coeffs.end());

		// Cost (local)
		double cost_local = 0.0;
		LPF = PF * L;
		for(size_t v = 0; v < V.size(); v++)
		{
			if(!inner[v] || V[v] == to_remove) continue;
			const double Mv = M[v];

			double cost_v = std::numeric_limits<double>::quiet_NaN();
			if(V[v] == to_keep)
			{
				const Eigen::VectorXd PZv = signals_L.col(to_keep) * (1 - alpha) + signals_L.col(to_remove) * alpha;
				cost_v = Mv * (PZv - (1.0 / Mv) * LPF.col(v)).squaredNorm();
			}
			else
				cost_v = Mv * (signals_L.col(V[v]) - (1.0 / Mv) * LPF.col(v)).squaredNorm();

			cost_local += cost_v;
			diff.emplace_back(V[v], cost_v);
		}

		// Finalize
		const double cost = cost_local - prev_cost_local;
		return { pos, cost };
	};

	if(optimal_pos)
	{
		std::vector<std::pair<uint32_t, double>> diff[3];
		std::pair<Eigen::Vector3d, double> cost[3] = { eval(0.0, diff[0]), eval(0.5, diff[1]), eval(1.0, diff[2]) };
		const Eigen::Vector3d c = fit_polynom2({ 0.0, 0.5, 1.0 }, { cost[0].second, cost[1].second, cost[2].second });
		const double beta = min_polynom2(c);

		computed_collapse[h].diff.clear();
		std::pair<Eigen::Vector3d, double> q = { {}, std::numeric_limits<double>::infinity() };
		const double tol = 0.005;
		if(beta > tol && std::abs(beta - 0.5) > tol && beta < 1 - tol)
		{
			computed_collapse[h].alpha = beta;
			q = eval(beta, computed_collapse[h].diff);
		}

		for(unsigned int k = 0; k < 3; k++)
			if(cost[k].second < q.second)
			{
				q = cost[k];
				computed_collapse[h].diff = std::move(diff[k]);
				computed_collapse[h].alpha = 0.5 * double(k);
			}

		collapse_revision[h] = mesh_revision;
		return q;
	}
	else
	{
		computed_collapse[h].alpha = 0.5;
		computed_collapse[h].diff.clear();
		collapse_revision[h] = mesh_revision;
		return eval(0.5, computed_collapse[h].diff);
	}
}

unsigned int metric_lowpass::post_collapse(uint32_t h, uint32_t to_keep, uint32_t to_remove)
{
	const double alpha = computed_collapse[h].alpha;
	const std::vector<std::pair<uint32_t, double>>& diff = computed_collapse[h].diff;

	// Update signals
	signals.col(to_keep) = signals.col(to_keep) * (1 - alpha) + signals.col(to_remove) * alpha;
	signals_L.col(to_keep) = signals_L.col(to_keep) * (1 - alpha) + signals_L.col(to_remove) * alpha;
	signals.col(to_remove).setConstant(std::numeric_limits<double>::quiet_NaN());
	signals_L.col(to_remove).setConstant(std::numeric_limits<double>::quiet_NaN());

	// Update projection
	projection.row(to_keep) = projection.row(to_keep) * (1 - alpha) + projection.row(to_remove) * alpha;
	projection.prune([&](int row, int /*col*/, double /*val*/){ return row != (int)to_remove; });
	mapping[to_remove] = to_keep;

	// Update faces
	k_ring(*connectivity, 1, to_keep, [&](const half_edge& h)
	{
		uint32_t& f = half_edge_triangle[h.index];
		mesh::triangle& t = triangles[f];
		unsigned int num_keep = 0;
		for(unsigned int j = 0; j < 3; j++)
		{
			if(t[j] == to_remove) t[j] = to_keep;
			if(t[j] == to_keep) num_keep++;
		}
		if(num_keep > 1)
		{
			t[0] = t[1] = t[2] = uint32_t(-1);
			f = half_edge_triangle[h.next().index];
		}
		else
		{
			cotangents[f] = triangle_cotangents(object->vertices[t[0]], object->vertices[t[1]], object->vertices[t[2]]);
			areas[f] = triangle_area(object->vertices[t[0]], object->vertices[t[1]], object->vertices[t[2]]);
		}
	});

	// Update costs & revisions
	mesh_revision++;
	for(const std::pair<uint32_t, double>& d : diff)
	{
		norms[d.first] = d.second;
		vertex_revision[d.first] = mesh_revision;
	}
	vertex_revision[to_remove] = uint32_t(-1);
	total_cost = norms.sum();

	// We need to update the 3-ring to be sure all modified collapses are updated
	computed_collapse[h] = edge_continuation();
	return 3;
}

double metric_lowpass::cost_flip(uint32_t h) const
{
	// Prepare continuation
	computed_flip[h].clear();
	flip_revision[h] = mesh_revision;

	// Vertices: [old edge, new edge]
	half_edge hx = connectivity->handle(h);
	uint32_t vertices[4] = { uint32_t(-1), uint32_t(-1), uint32_t(-1), uint32_t(-1) };
	vertices[0] = hx.vertex();
	vertices[1] = hx.opposite().vertex();
	vertices[2] = hx.next().vertex();
	vertices[3] = hx.opposite().next().vertex();
	const uint32_t tri1 = half_edge_triangle[hx.index];
	const uint32_t tri2 = half_edge_triangle[hx.opposite().index];

	// Anti-cycle criterion
#if 0
	const double len_old = (object->vertices[vertices[0]] - object->vertices[vertices[1]]).norm();
	const double len_new = (object->vertices[vertices[2]] - object->vertices[vertices[3]]).norm();
	if(len_new >= len_old) return std::numeric_limits<double>::infinity();
#endif
#if 0
	auto circumscribed_sphere = [&](uint32_t iA, uint32_t iB, uint32_t iC) -> weighted_position
	{
		const Eigen::Vector3d& vA = object->vertices[iA];
		const Eigen::Vector3d& vB = object->vertices[iB];
		const Eigen::Vector3d& vC = object->vertices[iC];
		const Eigen::Vector3d ea = vA - vC;
		const Eigen::Vector3d eb = vB - vC;
		const Eigen::Vector3d nn = ea.cross(eb);
		const Eigen::Vector3d p = vC + (ea.squaredNorm() * eb - eb.squaredNorm() * ea).cross(nn) / (2 * nn.squaredNorm());
		return { p, (vA - p).norm() };
	};
	auto in_sphere = [&](uint32_t v, const weighted_position& sph)
	{
		return (object->vertices[v] - sph.position).norm() <= sph.weight;
	};
	const weighted_position sph1 = circumscribed_sphere(vertices[0], vertices[2], vertices[3]);
	const weighted_position sph2 = circumscribed_sphere(vertices[3], vertices[2], vertices[1]);
	if(in_sphere(vertices[1], sph1) || in_sphere(vertices[0], sph2))
		return std::numeric_limits<double>::infinity();
#endif

	// Get local mesh
	std::vector<uint32_t> V, F;
	std::unordered_map<uint32_t, uint32_t> V_inv;
	{
		auto get_VF = [&](const half_edge& h)
		{
			V.push_back(h.vertex());
			F.push_back(half_edge_triangle[h.index]);
		};
		auto on_1_ring = [this](uint32_t v, auto f)
		{
			half_edge hh = connectivity->handle(connectivity->vertex_half_edge(v));
			uint32_t start = uint32_t(-1);
			while(hh.index != start && hh.is_valid())
			{
				if(start == uint32_t(-1)) start = hh.index;
				f(hh);
				f(hh = hh.next());
				f(hh = hh.next());
				hh = hh.opposite();
			}
		};
		for(uint32_t v : vertices)
			on_1_ring(v, get_VF);
		sort_unique_inplace(V);
		sort_unique_inplace(F);

		for(size_t t = 0; t < V.size(); t++)
			V_inv[V[t]] = (uint32_t)t;
	}

	// Previous local cost
	double prev_cost_local = 0;
	for(uint32_t v : vertices)
		prev_cost_local += norms[v];

	// Signals
	Eigen::MatrixXd PF(signals.rows(), V.size());
	for(size_t v = 0; v < V.size(); v++)
		PF.col(v) = signals.col(V[v]);

	// Compute mass and laplacian
	Eigen::VectorXd M = Eigen::VectorXd::Zero(V.size());
	std::vector<Eigen::Triplet<double>> coeffs;
	for(uint32_t f : F)
	{
		mesh::triangle tri = triangles[f];
		if(f == tri1)
		{
			for(uint32_t& w : tri.idx)
				if(w == vertices[1])
					w = vertices[3];
		}
		else if(f == tri2)
		{
			for(uint32_t& w : tri.idx)
				if(w == vertices[0])
					w = vertices[2];
		}
		const bool modified = f == tri1 || f == tri2;

		// Cotangents and area
		Eigen::Vector3d c = cotangents[f];
		double a = areas[f];
		if(modified)
		{
			const Eigen::Vector3d& v0 = object->vertices[tri[0]];
			const Eigen::Vector3d& v1 = object->vertices[tri[1]];
			const Eigen::Vector3d& v2 = object->vertices[tri[2]];
			c = triangle_cotangents(v0, v1, v2);
			a = triangle_area(v0, v1, v2);
		}

		// Local indices
		const uint32_t local_tri[3] = { V_inv.at(tri[0]), V_inv.at(tri[1]), V_inv.at(tri[2]) };

		// Mass
		for(uint32_t v : local_tri)
			M[v] += a / 3;

		// Laplacian
		for(unsigned int j = 0; j < 3; j++)
		{
			coeffs.emplace_back(local_tri[j], local_tri[j], c[(j + 1) % 3] + c[(j + 2) % 3]);
			coeffs.emplace_back(local_tri[j], local_tri[(j + 1) % 3], -c[(j + 2) % 3]);
			coeffs.emplace_back(local_tri[j], local_tri[(j + 2) % 3], -c[(j + 1) % 3]);
		}
	}
	Eigen::SparseMatrix<double> L(V.size(), V.size());
	L.setFromTriplets(coeffs.begin(), coeffs.end());

	// Cost (local)
	double cost_local = 0.0;
	const Eigen::MatrixXd LF = PF * L;
	for(uint32_t v : vertices)
	{
		const uint32_t w = (uint32_t)std::distance(V.begin(), std::find(V.begin(), V.end(), v));
		const double Mv = M[w];
		const double cost_v = Mv * (signals_L.col(v) - LF.col(w) * (1.0 / Mv)).squaredNorm();
		cost_local += cost_v;
		computed_flip[h].emplace_back(v, cost_v);
	}

	// If there are degenerate triangles, cotangents will be either +inf or -inf, which means
	// that when summed, there will be NaNs in the laplacian matrix. Since we don't want flips
	// that create degenerate triangles, simply make the cost +inf
	if(std::isnan(cost_local)) return std::numeric_limits<double>::infinity();

	// Return delta cost, and only allow flips that will decrease the global cost
	const double cost = cost_local - prev_cost_local;
	return cost < 0 ? cost : std::numeric_limits<double>::infinity();
}

unsigned int metric_lowpass::post_flip(uint32_t h)
{
	// Update faces
	half_edge hx = connectivity->handle(h);
	half_edge ho = hx.opposite();
	const uint32_t t1 = half_edge_triangle[hx.index];
	const uint32_t t2 = half_edge_triangle[ho.index];
	half_edge_triangle[hx.next().index] = t1;
	half_edge_triangle[hx.next().next().index] = t1;
	half_edge_triangle[ho.next().index] = t2;
	half_edge_triangle[ho.next().next().index] = t2;
	triangles[t1].set(hx.vertex(), hx.next().vertex(), hx.next().next().vertex());
	triangles[t2].set(ho.vertex(), ho.next().vertex(), ho.next().next().vertex());

	// Update costs & revisions
	mesh_revision++;
	for(const std::pair<uint32_t, double>& d : computed_flip[h])
	{
		norms[d.first] = d.second;
		vertex_revision[d.first] = mesh_revision;
	}
	total_cost = norms.sum();

	// While only 1-ring is needed to update other flips, we need to update the 3-ring to be sure all modified collapses are updated
	computed_flip[h].clear();
	return 3;
}

metric_lowpass::~metric_lowpass()
{
	// Remove all empty rows, keep the order
	static_assert(decltype(projection)::IsRowMajor, "Projection needs to be row major");
	projection.prune(0.0);
	assert(!projection.innerNonZeroPtr() && "Need compressed mode");
	const Eigen::Index count = projection.rows();
	int* outer_start = projection.outerIndexPtr();
	auto is_empty = [&](int row) { return outer_start[row + 1] == outer_start[row]; };

	int empty_row = 0;
	while(empty_row < count && !is_empty(empty_row))
		empty_row++;
	for(int row = empty_row; row < count; row++)
	{
		if(is_empty(row)) continue;
		outer_start[++empty_row] = outer_start[row + 1];
	}

	projection.conservativeResize(empty_row, projection.cols());
	projection.makeCompressed();
	save_matrix_ascii("projection.txt", projection);

	// Final mapping
	for(size_t v = 0; v < mapping.size(); v++)
	{
		if(mapping[v] == uint32_t(-1)) continue;
		while(mapping[mapping[v]] != uint32_t(-1))
			mapping[v] = mapping[mapping[v]];
	}
	std::vector<uint32_t> new_index(mapping.size(), uint32_t(-1));
	uint32_t idx = 0;
	for(size_t v = 0; v < mapping.size(); v++)
		if(mapping[v] == uint32_t(-1))
			new_index[v] = idx++;
	for(size_t v = 0; v < mapping.size(); v++)
		mapping[v] = new_index[mapping[v] != uint32_t(-1) ? mapping[v] : v];
	std::ofstream fm("mapping.txt");
	for(uint32_t w : mapping)
		fm << w << '\n';
}

reduction_metric::serialized_parameters metric_lowpass::save() const
{
	serialized_parameters p;
	serialize(p, "eigenvectors",  num_eigenvectors);
	serialize(p, "optimal_pos",   optimal_pos);
	return p;
}

void metric_lowpass::load(const serialized_parameters& p)
{
	deserialize(p, "eigenvectors",  num_eigenvectors);
	deserialize(p, "optimal_pos",   optimal_pos);
}
