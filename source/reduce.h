#pragma once

struct mesh;
class half_edge_connectivity;
struct reduction_metric;

void remove_standalone_vertices(mesh& obj, const half_edge_connectivity& h);

struct reduce_options
{
	bool directed_edges = false;
	bool flip_edges = false;
};

void reduce_stream(mesh& obj, size_t target_size, reduction_metric& metric, const reduce_options& options = {});

void reduce_stream_noflip(mesh& obj, size_t target_size, reduction_metric& metric, const reduce_options& options = {});

void optimize_flips(mesh& obj, reduction_metric& metric);
