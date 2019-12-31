#include <charconv>
#include <memory>
#include <chrono>
#include "mesh.h"

#include "reduce.h"
#include "metric_lowpass.h"
#include "metric_quadrics.h"
#include "metric_edge_length.h"

template<typename C, typename D>
double elapsed(const std::chrono::time_point<C, D>& start, const std::chrono::time_point<C, D>& end)
{
	return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;
}

using metric_create_f = std::unique_ptr<reduction_metric>();

template<typename T>
std::unique_ptr<reduction_metric> make_metric() { return std::make_unique<T>(); }

constexpr const std::pair<const char*, metric_create_f*>  metric_register[] =
{
	{ "edge_length", make_metric<metric_edge_length> },
	{ "qem",         make_metric<metric_quadrics> },
	{ "lowpass",     make_metric<metric_lowpass> }
};

std::unique_ptr<reduction_metric> make_metric(const char* name, const reduction_metric::serialized_parameters& params)
{
	for(const std::pair<const char*, metric_create_f*>& m : metric_register)
		if(ci_equal(name, m.first))
		{
			std::unique_ptr<reduction_metric> metric = (*m.second)();
			metric->load(params);
			return metric;
		}
	return nullptr;
}

struct command_line
{
	// Files
	const char* filename_in = nullptr;
	const char* filename_out = nullptr;

	// General options
	const char* metric_type = nullptr;
	unsigned int num_vertices = 0;
	bool normalize_area = true;
	bool enable_flips = false;

	// Metric
	reduction_metric::serialized_parameters metric_params;

	command_line(unsigned int argc, char** argv)
	{
		if(argc == 2)
		{
			magic(argc, argv);
			return;
		}

		if(argc < 6)
		{
			printf("ERROR: not enough arguments\n");
			print_help();
			exit(0);
		}

		filename_in = argv[1];
		filename_out = argv[2];
		std::from_chars(argv[3], argv[3] + strlen(argv[3]), num_vertices);

		enable_flips = false;
		if(ci_equal(argv[4], "stream+flips"))
			enable_flips = true;
		else if(!ci_equal(argv[4], "stream"))
		{
			printf("ERROR: unknown reduction algorithm: %s\n", argv[4]);
			print_help();
			exit(0);
		}

		metric_type = argv[5];

		for(unsigned int a = 6; a < argc; a++)
		{
			std::string pv = argv[a];
			size_t sep = pv.find('=');
			if(sep < pv.size())
				metric_params.emplace_back(pv.substr(0, sep), pv.substr(sep + 1));
			else
				printf(
					"WARNING: unknown parameter \"%s\",\n"
					"         parameters should have the form <name>=<value> (no spaces)\n", argv[a]);
		}
	}

	void magic(unsigned int, char** argv)
	{
		printf("Automagically filling missing parameters\n");
		filename_in = argv[1];
		filename_out = "result.ply";
		num_vertices = 600;
		normalize_area = true;

		auto stream = [&](const char* name, reduction_metric::serialized_parameters p)
		{
			metric_type = name;
			metric_params = std::move(p);
		};

		switch(9)
		{
		case 0: stream("qem", {}); break;
		case 2: stream("local", { { "eigenvectors", "4" } }); break;
		case 4: stream("global", { { "eigenvectors", "100" } }); break;
		case 6: stream("edge_length", {}); break;
		case 7: stream("noise", {}); break;
		case 9: stream("lowpass", { { "optimal_pos", "true" }, { "eigenvectors", "100" } }); break;
		case 10: stream("lowpass", { { "optimal_pos", "true" }, { "eigenvectors", "100" } }); enable_flips = true; break;
		}
	}

	void print_help() const
	{
		printf("Usage:\n"
			"spectral-collapsing <file_in> <file_out> <num_vertices> <algo> <metric> [<param>=<value>]*\n"
			"\n"
			"\t<file_in>       : mesh file to load and process\n"
			"\t<file_out>      : mesh file to save\n"
			"\t<num_vertices>  : target number of vertices in the final mesh\n"
			"\t<algo>          : reduction algorithm, can be \"stream\" or \"batch\"\n"
			"\t<metric>        : reduction metric, can be:\n");
		for(const std::pair<const char*, metric_create_f*>& m : metric_register)
			printf("\t                  - %s\n", m.first);
		printf("\t<param>=<value> : a metric parameter (no spaces)\n");
	}

	void print_parameters() const
	{
		printf("Parameters:\n"
			"\tFilename (in)  : %s\n"
			"\tFilename (out) : %s\n"
			"\tNormalize area : %s\n"
			"\tNum vertices   : %u\n"
			"\tReduction      : %s\n"
			"\tMetric         : %s\n",
			filename_in, filename_out, normalize_area ? "yes" : "no",
			num_vertices, enable_flips ? "stream [collapses + flips]" : "stream [collapses]",
			metric_type);
	}
};

mesh load_mesh(const char* filename)
{
	using clock = std::chrono::high_resolution_clock;

	const clock::time_point ts = clock::now();
	mesh m(filename);
	const clock::time_point te = clock::now();

	if(m.empty())
	{
		printf("Cannot read triangle mesh '%s'\n", filename);
		exit(1);
	}

	printf("Input mesh:\n"
		"\tvertices  : %zu\n"
		"\tedges     : %zu\n"
		"\ttriangles : %zu\n",
		m.num_vertices(), num_edges(m), m.num_triangles());

	printf("\tload time : %.6f s\n", elapsed(ts, te));
	return m;
}

double reduce(mesh& m, unsigned int num_vertices, bool enable_flips, reduction_metric& metric)
{
	using clock = std::chrono::high_resolution_clock;

	reduce_options opts;
	opts.directed_edges = false;
	opts.flip_edges = enable_flips;

	const clock::time_point ts = clock::now();
	if(!enable_flips) reduce_stream_noflip(m, num_vertices, metric, opts);
	else reduce_stream(m, num_vertices, metric, opts);
	const clock::time_point te = clock::now();
	const double duration = elapsed(ts, te);
	printf("Reduction time: %.6f s\n", duration);
	return duration;
}

void print_metric_parameters(const reduction_metric& metric)
{
	const reduction_metric::serialized_parameters params = metric.save();
	size_t len = 0;
	for(const std::pair<std::string, std::string>& p : params)
		len = std::max(len, p.first.size());

	printf("Metric parameters:\n");
	for(const std::pair<std::string, std::string>& p : params)
		printf("\t%-*s: %s\n", int(len + 1), p.first.c_str(), p.second.c_str());
	if(params.empty()) printf("\t<none>\n");
}

int main(unsigned int argc, char** argv)
{
	const command_line cmd(argc, argv);
	cmd.print_parameters();

	std::unique_ptr<reduction_metric> metric = make_metric(cmd.metric_type, cmd.metric_params);
	if(metric) print_metric_parameters(*metric);
	else
	{
		printf("Invalid metric\n");
		exit(1);
	}

	mesh m = load_mesh(cmd.filename_in);
	const double area = total_area(m);
	printf("\tarea      : %.12f%s\n", area, cmd.normalize_area ? " (normalized)" : "");
	if(cmd.normalize_area) m.scale(1.0 / sqrt(area));

	reduce(m, cmd.num_vertices, cmd.enable_flips, *metric);

	if(cmd.normalize_area) m.scale(sqrt(area));
	m.save(cmd.filename_out);
}
