#include "geometry.h"
#include <execution>
#include "debugbreak.h"
#include "ranges.h"

#pragma warning(push)
#pragma warning(disable: 4127)
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#pragma warning(pop)

Eigen::Vector3d triangle_qem_weights(qem_face_weight type, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
{
	switch(type)
	{
	case qem_face_weight::constant: return Eigen::Vector3d::Constant(1.0);
	case qem_face_weight::area:     return Eigen::Vector3d::Constant(triangle_area(v0, v1, v2));
	case qem_face_weight::angle:    return triangle_angles(v0, v1, v2);
	default: printf("Invalid weight type\n"); debug_break(); return Eigen::Vector3d::Zero();
	}
}

Eigen::Matrix4d vertex_quadric(const mesh& obj, const half_edge_connectivity& connectivity, uint32_t v0, qem_face_weight type)
{
	Eigen::Matrix4d quadric = Eigen::Matrix4d::Zero();
	half_edge h = connectivity.handle(connectivity.vertex_half_edge(v0));
	uint32_t v1 = uint32_t(-1), v2 = uint32_t(-1);
	uint32_t end = uint32_t(-1);
	while(h.index != end && h.is_valid())
	{
		if(end == uint32_t(-1)) end = h.index;
		v1 = h.vertex();
		h = h.next();
		v2 = h.vertex();
		h = h.next().opposite();
		// Compute quadric
		const Eigen::Vector3d& p0 = obj.vertices[v0];
		const Eigen::Vector3d& p1 = obj.vertices[v1];
		const Eigen::Vector3d& p2 = obj.vertices[v2];
		const Eigen::Vector3d n = (p1 - p0).cross(p2 - p0).normalized();
		const Eigen::Vector4d p(n[0], n[1], n[2], -n.dot(p0));
		const Eigen::Matrix4d Q = p * p.transpose();
		const Eigen::Vector3d w = triangle_qem_weights(type, p0, p1, p2);
		quadric += Q * w[0];
	}
	return quadric;
}

std::vector<Eigen::Matrix4d> geometry_quadrics(const mesh& obj, qem_face_weight type)
{
	std::vector<Eigen::Matrix4d> quadrics(obj.vertices.size(), Eigen::Matrix4d::Zero());
	for(uint32_t t = 0; t < obj.num_triangles(); t++)
	{
		const mesh::triangle& tri = obj.triangles[t];
		const Eigen::Vector3d& p0 = obj.vertices[tri[0]];
		const Eigen::Vector3d& p1 = obj.vertices[tri[1]];
		const Eigen::Vector3d& p2 = obj.vertices[tri[2]];
		const Eigen::Vector3d n = (p1 - p0).cross(p2 - p0).normalized();
		const Eigen::Vector4d p(n[0], n[1], n[2], -n.dot(p0));
		const Eigen::Matrix4d Q = p * p.transpose();
		const Eigen::Vector3d w = triangle_qem_weights(type, p0, p1, p2);
		quadrics[tri[0]] += Q * w[0];
		quadrics[tri[1]] += Q * w[1];
		quadrics[tri[2]] += Q * w[2];
	}
	return quadrics;
}

std::pair<Eigen::Vector3d, double> qem_merge(const Eigen::Matrix4d& Q1, const Eigen::Matrix4d& Q2, const Eigen::Vector3d& fallback)
{
	const Eigen::Matrix4d Q = Q1 + Q2;
	Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
	bool invertible = true;
	Eigen::Matrix3d(Q.block(0, 0, 3, 3)).computeInverseWithCheck(A, invertible);
	Eigen::Vector3d target;
	if(invertible)
		target = -A * Q.block(0, 3, 3, 1);
	else
		target = fallback;
	const Eigen::Vector4d hg_tgt = { target[0], target[1], target[2], 1.0 };
	const double qem = hg_tgt.transpose() * Q * hg_tgt;
	return { target, qem };
}

double total_area(const mesh& obj)
{
	double area = 0.0;
	for(const mesh::triangle& t : obj.triangles)
		area += triangle_area(obj.vertices[t[0]], obj.vertices[t[1]], obj.vertices[t[2]]);
	return area;
}

double triangle_area(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
	const double e[3] = { (a - b).norm(), (a - c).norm(), (b - c).norm() };
	const double p = (e[0] + e[1] + e[2]) / 2;
	return std::sqrt(std::max(p * (p - e[0]) * (p - e[1]) * (p - e[2]), 0.0));
}

std::vector<double> triangle_areas(const mesh& obj)
{
	std::vector<double> areas;
	areas.reserve(obj.triangles.size());
	for(const mesh::triangle& t : obj.triangles)
		areas.push_back(triangle_area(obj.vertices[t[0]], obj.vertices[t[1]], obj.vertices[t[2]]));
	return areas;
}

Eigen::Vector3d triangle_cosines(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
	const Eigen::Vector3d e0 = (b - a).normalized();
	const Eigen::Vector3d e1 = (c - b).normalized();
	const Eigen::Vector3d e2 = (a - c).normalized();
	return Eigen::Vector3d(e0.dot(-e2), e1.dot(-e0), e2.dot(-e1)).cwiseMax(-1.0).cwiseMin(1.0);
}

Eigen::Vector3d triangle_angles(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
	const Eigen::Vector3d d = triangle_cosines(a, b, c);
	return { std::acos(d[0]), std::acos(d[1]), std::acos(d[2]) };
}

Eigen::Vector3d triangle_cotangents(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
	const Eigen::Vector3d u = b - a;
	const Eigen::Vector3d v = c - a;
	const Eigen::Vector3d w = c - b;
	return { u.dot(v) / u.cross(v).norm(), w.dot(-u) / w.cross(-u).norm(), v.dot(w) / v.cross(w).norm() };
}

std::vector<Eigen::Vector3d> triangle_cotangents(const mesh& obj)
{
	std::vector<Eigen::Vector3d> cotangents;
	cotangents.reserve(obj.triangles.size());
	for(const mesh::triangle& t : obj.triangles)
		cotangents.push_back(triangle_cotangents(obj.vertices[t[0]], obj.vertices[t[1]], obj.vertices[t[2]]));
	return cotangents;
}

Eigen::Vector3d fit_polynom2(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
{
	Eigen::Matrix3d S;
	for(unsigned int k = 0; k < 3; k++)
	{
		double t = 1.0;
		for(unsigned int d = 0; d < 3; d++)
		{
			S(k, d) = t;
			t *= x(k);
		}
	}
	return S.colPivHouseholderQr().solve(y);
}

double eval_polynom2(const Eigen::Vector3d& c, double x)
{
	return c[0] + c[1] * x + c[2] * x * x;
}

double min_polynom2(const Eigen::Vector3d& c)
{
	return -c[1] / (2 * c[2]);
}

Eigen::SparseMatrix<double> laplacian_geometric(size_t num_vertices, const std::vector<mesh::triangle>& triangles, const std::vector<Eigen::Vector3d>& cotangents)
{
	std::vector<Eigen::Triplet<double>> coeffs;
	for(size_t t = 0; t < triangles.size(); t++)
	{
		const Eigen::Vector3d& c = cotangents[t];
		const mesh::triangle& f = triangles[t];
		for(unsigned int j = 0; j < 3; j++)
		{
			coeffs.emplace_back(f[j], f[j], c[(j + 1) % 3] + c[(j + 2) % 3]);
			coeffs.emplace_back(f[j], f[(j + 1) % 3], -c[(j + 2) % 3]);
			coeffs.emplace_back(f[j], f[(j + 2) % 3], -c[(j + 1) % 3]);
		}
	}

	Eigen::SparseMatrix<double> L(num_vertices, num_vertices);
	L.setFromTriplets(coeffs.begin(), coeffs.end());
	return L;
}

Eigen::VectorXd mass_barycentric(size_t num_vertices, const std::vector<mesh::triangle>& triangles, const std::vector<double>& areas)
{
	Eigen::VectorXd M = Eigen::VectorXd::Zero(num_vertices);
	for(size_t t = 0; t < triangles.size(); t++)
		for(unsigned int k = 0; k < 3; k++)
			M[triangles[t][k]] += areas[t] / 3;
	return M;
}

std::vector<std::pair<uint32_t, uint32_t>> undirected_edges(const mesh& obj)
{
	using edge_t = std::pair<uint32_t, uint32_t>;
	std::vector<edge_t> edges;
	auto add_edge = [&](uint32_t a, uint32_t b) { edges.emplace_back(std::min(a, b), std::max(a, b)); };
	for(size_t t = 0; t < obj.triangles.size(); t++)
	{
		add_edge(obj.triangles[t][0], obj.triangles[t][1]);
		add_edge(obj.triangles[t][1], obj.triangles[t][2]);
		add_edge(obj.triangles[t][2], obj.triangles[t][0]);
	}
	sort_unique_inplace(std::execution::par_unseq, edges, [](const edge_t& a, const edge_t& b)
	{
		return a.first != b.first ? a.first < b.first : a.second < b.second;
	});
	return edges;
}

size_t num_edges(const mesh& obj)
{
	return undirected_edges(obj).size();
}

bool smallest_eigenvectors(const Eigen::SparseMatrix<double>& L, const Eigen::VectorXd& M, unsigned int num,
	Eigen::VectorXd* out_eigenvalues, Eigen::MatrixXd* out_eigenvectors)
{
	if(!out_eigenvalues && !out_eigenvectors) return false;

	/// 1. Instead of solving the generalized problem (L X = M X D), solve the simpler problem (H Y = Y D), where
	///     - H = N L N
	///     - N = M^(-1/2)
	///     - X = N Y
	///    This allow for faster solve.
	///
	/// 2. Since L is symmetric positive definite, all its eigenvalues are in R+. So, instead of specifying we want
	///    the smallest eigenvalues, we can ask for the largest eigenvalues that are close to 0, or any negative number.
	///    Using a negative epsilon (-1e-6) instead of 0 makes the convergence more robust, and also faster.
	const Eigen::VectorXd N = M.cwiseSqrt().cwiseInverse();
	const Eigen::SparseMatrix<double> H = N.asDiagonal() * L * N.asDiagonal();

	const unsigned int ncv = std::min(3 * num, (unsigned int)L.rows());
	Spectra::SparseSymShiftSolve<double> op(H);
	Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, decltype(op)> solver(&op, num, ncv, -1e-6);
	solver.init();
	solver.compute(1000, 1e-10, Spectra::SMALLEST_ALGE);
	if(solver.info() != Spectra::SUCCESSFUL)
	{
		printf("Failed to compute eigenvectors\n");
		return false;
	}

	if(out_eigenvalues) *out_eigenvalues = solver.eigenvalues();
	if(out_eigenvectors) *out_eigenvectors = N.asDiagonal() * solver.eigenvectors();
	return true;
}

void get_edge_2_ring(const half_edge_connectivity& connectivity, const std::vector<uint32_t>& half_edge_triangle, uint32_t h,
	std::vector<uint32_t>& V, std::vector<uint32_t>& F)
{
	const auto [v1, v2] = connectivity.edge_vertices(h);
	V.push_back(v1);
	V.push_back(v2);

	auto add_vert_ring1 = [&](uint32_t vert)
	{
		half_edge he = connectivity.handle(connectivity.vertex_half_edge(vert));
		uint32_t index = uint32_t(-1);
		while(he.is_valid() && he.index != index)
		{
			if(index == uint32_t(-1)) index = he.index;
			const uint32_t w = he.vertex();
			if(w != v1 && w != v2) V.push_back(w);
			F.push_back(half_edge_triangle[he.index]);
			half_edge hx = he.next().next().opposite();
			if(!hx.is_valid()) V.push_back(he.next().vertex());
			he = hx;
		}
	};

	add_vert_ring1(v1);
	add_vert_ring1(v2);
	const size_t end_ring1 = V.size();
	for(size_t v = 2; v < end_ring1; v++)
		add_vert_ring1(V[v]);

	sort_unique_inplace(V);
	sort_unique_inplace(F);
}
