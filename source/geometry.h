#pragma once

#include "mesh.h"
#pragma warning(push)
#pragma warning(disable: 4127)
#include <Eigen/Sparse>
#pragma warning(pop)

enum class qem_face_weight
{
	constant,
	area,
	angle
};

/// Return the weights of the vertices of the triangle
Eigen::Vector3d triangle_qem_weights(qem_face_weight type, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

/// Compute the vertex quadric, by summing plane quadrics
Eigen::Matrix4d vertex_quadric(const mesh& obj, const half_edge_connectivity& connectivity, uint32_t vertex, qem_face_weight weight = qem_face_weight::area);

/// Compute 1 quadric per vertex, by summing plane quadrics
std::vector<Eigen::Matrix4d> geometry_quadrics(const mesh& obj, qem_face_weight weight = qem_face_weight::area);

/// Merge 2 quadrics, returning the optimal point and the value of Q1+Q2 at this point
std::pair<Eigen::Vector3d, double> qem_merge(const Eigen::Matrix4d& Q1, const Eigen::Matrix4d& Q2, const Eigen::Vector3d& fallback);

/// Return the area of the mesh (sum of triangle areas)
double total_area(const mesh& obj);

/// Return the area of the specified triangle
double triangle_area(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

/// Return the area of all the triangles
std::vector<double> triangle_areas(const mesh& obj);

/// Return (cos a, cos b, cos c) for a triangle (a, b, c)
Eigen::Vector3d triangle_cosines(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

/// Return angles for a triangle (a, b, c)
Eigen::Vector3d triangle_angles(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

/// Return (cot a, cot b, cot c) for a triangle (a, b, c)
Eigen::Vector3d triangle_cotangents(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

/// Return triangle_cotangents for each triangle
std::vector<Eigen::Vector3d> triangle_cotangents(const mesh& obj);

/// Return a polynom of degree 2 passing though the given points (P(x) = y)
// result(i) = coeff of x^i
Eigen::Vector3d fit_polynom2(const Eigen::Vector3d& x, const Eigen::Vector3d& y);

/// Return P(x)
double eval_polynom2(const Eigen::Vector3d& c, double x);

/// Return the x such that P(x) minimum
double min_polynom2(const Eigen::Vector3d& c);

/// Convert a std::vector to a Eigen dynamic vector
inline Eigen::VectorXd to_eigen(const std::vector<double>& v)
{
	Eigen::VectorXd x(v.size());
	for(size_t t = 0; t < v.size(); t++) x[t] = v[t];
	return x;
}

/// Convert a std::vector of vectors to a Eigen dynamic matrix
inline Eigen::MatrixXd to_eigen(const std::vector<Eigen::Vector3d>& v)
{
	Eigen::MatrixXd x(v.size(), 3);
	for(size_t t = 0; t < v.size(); t++) x.row(t) = v[t];
	return x;
}

/// Return the laplacian matrix from cotangents, and without mass
Eigen::SparseMatrix<double> laplacian_geometric(size_t num_vertices, const std::vector<mesh::triangle>& triangles, const std::vector<Eigen::Vector3d>& cotangents);

/// Return the vector of mass from triangles and areas; mass(vertex) = 1/3 * areas of neighbouring triangles
Eigen::VectorXd mass_barycentric(size_t num_vertices, const std::vector<mesh::triangle>& triangles, const std::vector<double>& areas);

/// Return edges of the specified mesh, undirected
std::vector<std::pair<uint32_t, uint32_t>> undirected_edges(const mesh& obj);

/// Return the number of edges in the specified mesh
size_t num_edges(const mesh& obj);

/// Solve the generalized eigenvalue problem: L X = M X D:
// L is a real symmetric positive definite matrix
// M is a diagonal matrix (only the diagonal is given to this function)
// D is the diagonal matrix of eigenvalues (returned by its diagonal, if not nullptr)
// X is the eigenvectors matrix (returned if not nullptr)
// This solves for the smallest eigenvalues, and returns whether the convergence was successful
bool smallest_eigenvectors(const Eigen::SparseMatrix<double>& L, const Eigen::VectorXd& M, unsigned int num,
	Eigen::VectorXd* out_eigenvalues, Eigen::MatrixXd* out_eigenvectors);

/// Fill V and F with the 2-ring of edge h (Vertices and Faces indices, respectively)
void get_edge_2_ring(const half_edge_connectivity& connectivity, const std::vector<uint32_t>& half_edge_triangle, uint32_t h,
	std::vector<uint32_t>& V, std::vector<uint32_t>& F);
