#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace fem3d {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

struct Tet4 {
    std::array<std::size_t, 4> node{};
};

struct Mesh {
    std::vector<Vec3> points;
    std::vector<Tet4> tetrahedra;
};

struct ElementData {
    double volume = 0.0;
    std::array<Vec3, 4> shape_gradients{};
    std::array<std::array<double, 4>, 4> stiffness{};
};

struct CSRMatrix {
    std::size_t rows = 0;
    std::vector<std::size_t> row_offsets;
    std::vector<std::size_t> column_indices;
    std::vector<double> values;

    std::vector<double> multiply(const std::vector<double>& x) const;
    double diagonal(std::size_t row) const;
};

struct SolverOptions {
    double relative_tolerance = 1e-12;
    std::size_t max_iterations = 10000;
};

struct SolverResult {
    std::vector<double> solution;
    std::size_t iterations = 0;
    double residual_norm = 0.0;
    bool converged = false;
};

// Computes the absolute physical volume and constant P1 shape gradients.
ElementData tet4_element(const std::array<Vec3, 4>& points, double conductivity = 1.0);

// Assemble K_ij = integral(k grad(N_i) . grad(N_j)) dV into deterministic CSR.
CSRMatrix assemble(const Mesh& mesh, double conductivity = 1.0);

// Apply prescribed nodal temperatures while preserving symmetry (strong elimination).
void apply_dirichlet(const CSRMatrix& matrix,
                     const std::vector<double>& load,
                     const std::vector<std::pair<std::size_t, double>>& prescribed,
                     CSRMatrix& constrained_matrix,
                     std::vector<double>& constrained_load);

SolverResult solve_pcg(const CSRMatrix& matrix,
                       const std::vector<double>& rhs,
                       const SolverOptions& options = {});

void write_legacy_vtk(const std::string& path, const Mesh& mesh,
                      const std::vector<double>& nodal_temperature);

// Eight unit-cube corners and a center node, split into 12 tetrahedra.
Mesh make_unit_cube_sample();

} // namespace fem3d
