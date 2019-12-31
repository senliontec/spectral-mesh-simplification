Spectral mesh simplification
============================

The spectrum of the Laplace-Beltrami operator is instrumental for a number of geometric modeling applications, from processing to analysis. Recently, multiple methods were developed to retrieve an approximation of a shape that preserves its eigenvectors as much as possible, but these techniques output a subset of input points with no connectivity, which limits their potential applications. Furthermore, the obtained Laplacian results from an optimization procedure, implying its storage alongside the output selected points: focusing on keeping a mesh instead of an operator would allow to retrieve the latter using the standard cotangent formulation, enabling easier processing afterwards.

**Instead, we propose to simplify the input mesh using a spectrum-preserving mesh decimation scheme, so that the Laplacian computed on the simplified mesh is spectrally close to the one of the input mesh.**

Run
---

The tool whose code is in this repository is command-line only. The arguments are as follow:

    spectral-collapsing.exe <INPUT> <OUTPUT> <SIZE> <ALGO> <METRIC> <PARAMETERS>

where:
 - **INPUT** is the input mesh, in a format readable by libigl (that is, *ply*, *obj*, *off*, ...)
 - **OUTPUT** is the output mesh, in a format writable by libigl
 - **SIZE** is the number of vertices in the output mesh
 - **ALGO** is the type of algorithm to reduce the input mesh, it is either `stream` (the standard decimation algorithm) or `stream+flips` (the usual decimation but also with edge flips)
 - **METRIC** is the metric to use, it can be `qem` for the Quadric Error Metric of [Garland & Heckbert 1997], `lowpass` for our metric, or `edge_length` for a very simple metric in which shorter edges are collapsed first
 - **PARAMETERS** are the parameters of the metric, each one given following the format `name=value`, *without spaces*.

There are no parameters for the metric `edge_length`. For the other metrics, here are the parameters:

| Metric    | Parameter                  | Explanation                                                                                                                  |
|-----------|----------------------------|------------------------------------------------------------------------------------------------------------------------------|
| `lowpass` | `eigenvectors=<integer>`   | number of eigenvectors to preserve                                                                                           |
| `lowpass` | `optimal_pos=<true/false>` | enable the approximation of the minimizer on collapse edges (if false, the collapse position will be the center of the edge) |
| `qem`     | `recompute=<true/false>`   | fully recompute quadrics instead of interpolating them on collapse                                                           |
| `qem`     | `noflips=<true/false>`     | disallow collapses which would flip the triangle normals                                                                     |
| `qem`     | `keepborders=<true/false>` | avoid collapsing boundary edges                                                                                              |

Build
-----

This project requires C++17 to compile, with a working implementation of `<execution>` (such as VS2017). Dependencies are provided as submodules and should compile straightforwardly. The rest is the usual CMake setup:

    mkdir build
    cd build
    cmake .. -A x64
    make

The output will be in the folder `/bin`.
