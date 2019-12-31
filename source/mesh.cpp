#include "mesh.h"
#pragma warning(push, 0)
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#pragma warning(pop)

mesh::mesh(const char* filename)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	if(!igl::read_triangle_mesh(filename, V, F)) return;

	vertices.resize(V.rows());
	for(size_t v = 0; v < vertices.size(); v++)
		vertices[v] = V.row(v);

	triangles.resize(F.rows());
	for(size_t t = 0; t < triangles.size(); t++)
		triangles[t].set(F(t, 0), F(t, 1), F(t, 2));
}

void mesh::save(const std::string& filename) const
{
	igl::write_triangle_mesh(filename, matrix_vertices(), matrix_triangles());
}

Eigen::MatrixXd mesh::matrix_vertices() const
{
	Eigen::MatrixXd V(vertices.size(), 3);
	for(size_t v = 0; v < vertices.size(); v++)
		V.row(v) = vertices[v];
	return V;
}

Eigen::MatrixXi mesh::matrix_triangles() const
{
	Eigen::MatrixXi F(triangles.size(), 3);
	for(size_t t = 0; t < triangles.size(); t++)
		for(unsigned int k = 0; k < 3; k++)
			F(t, k) = triangles[t][k];
	return F;
}
