#pragma once

#pragma warning(push)
#pragma warning(disable: 4127)
#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma warning(pop)
#include <string>
#include <fstream>
#include <cctype>

template<typename Derived>
void save_matrix_ascii(const std::string& filename, const Eigen::DenseBase<Derived>& mat, char sep = ',')
{
	std::ofstream out(filename);
	out.precision(14);
	for(Eigen::Index r = 0; r < mat.rows(); r++)
	{
		for(Eigen::Index c = 0; c < mat.cols(); c++)
			if(c == 0) out << mat(r, c);
			else out << sep << mat(r, c);
		out << '\n';
	}
}

template<typename Scalar, int Options, typename StorageIndex>
void save_matrix_ascii(const std::string& filename, const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat, char sep = '\t')
{
	using SpMat = Eigen::SparseMatrix<Scalar, Options, StorageIndex>;
	std::ofstream out(filename);
	out.precision(14);
	bool max_pos = false;
	for(int outer = 0; outer < mat.outerSize(); outer++)
		for(SpMat::InnerIterator it(mat, outer); it; ++it)
		{
			out << (it.row() + 1) << sep << (it.col() + 1) << sep << it.value() << '\n';
			if(it.row() + 1 == mat.rows() && it.col() + 1 == mat.cols()) max_pos = true;
		}
	if(!max_pos) out << mat.rows() << sep << mat.cols() << sep << "0.0\n";
}

template<size_t ID>
void print_progress(const std::string& label, double interval = 1.0)
{
	const std::chrono::milliseconds dt(int(1000 * interval));
	using clock = std::chrono::high_resolution_clock;
	static clock::time_point ttt = clock::now() - 2 * dt;
	const clock::time_point tn = clock::now();
	if(tn - ttt > dt)
	{
		ttt = tn;
		printf("%s                        \r", label.c_str());
	}
}

inline bool ci_equal(const char* a, const char* b)
{
	if(!a || !b) return false;
	while(*a && *b && std::tolower(*a) == std::tolower(*b))
		a++, b++;
	return !*a && !*b;
}
