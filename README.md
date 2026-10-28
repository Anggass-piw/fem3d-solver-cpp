# fem3d-solver-cpp

A small, dependency-free C++17 educational MVP for steady heat conduction on linear tetrahedral (Tet4) meshes. The included executable solves a deterministic unit-cube example and writes results as legacy VTK unstructured-grid data.

## What it implements

- Tet4 volume, constant linear shape-function gradients, and element conduction stiffness.
- Deterministic CSR sparse matrix assembly.
- Strong Dirichlet temperature constraints using symmetric row/column elimination and right-hand-side correction.
- Jacobi-preconditioned conjugate gradient (PCG) for the resulting symmetric positive-definite system.
- Legacy ASCII VTK output with nodal `temperature` scalar data.
- A small CTest suite; no third-party libraries are required.

The steady equation is `-div(k grad(T)) = q`. This MVP assembles the conduction stiffness for constant, positive, scalar conductivity `k`; the CLI uses `k = 1` and zero volumetric source. Natural boundary conditions are homogeneous Neumann (zero heat flux) where no boundary condition is prescribed.

## Build

Requirements: CMake 3.16 or newer and a C++17 compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Run the sample

```sh
./build/fem3d_solver
# or choose a destination:
./build/fem3d_solver --output cube.vtk
```

The mesh is the unit cube with its eight corners and one center point, divided into 12 tetrahedra. The solver imposes `T=0` on the `x=0` face and `T=1` on the `x=1` face; all other faces are insulated. The exact linear solution is `T=x`, so the center temperature is `0.5`. Open the generated `.vtk` file in ParaView or another legacy-VTK reader. The CLI currently uses this built-in example; mesh-file input is not implemented.

## API sketch

Include `fem3d/fem3d.hpp`, create a `fem3d::Mesh`, call `fem3d::assemble(mesh, conductivity)`, apply constraints with `fem3d::apply_dirichlet`, then call `fem3d::solve_pcg`. `fem3d::write_legacy_vtk` exports the nodal result. `make_unit_cube_sample()` provides the built-in demonstration mesh.

## Assumptions and limitations

This is an **educational small-problem MVP**, not a validated general-purpose simulation package. It is **not validated for engineering, design, regulatory, or safety use**. Use independent verification and qualified review before relying on any numerical result.

The implementation supports only linear four-node tetrahedra, isotropic constant scalar conductivity, steady-state conduction, nodal temperature Dirichlet constraints, and zero-flux natural boundaries. It has no mesh reader, source/load integration, transient or nonlinear physics, anisotropic/material fields, contact/interface conditions, adaptive refinement, preconditioning beyond Jacobi, parallel assembly, or rigorous mesh-quality diagnostics. A valid constrained problem must remove the constant-temperature nullspace; otherwise PCG reports breakdown or a nonpositive diagonal. Near-degenerate tetrahedra may be rejected as singular. The sparse representation is assembled from in-memory element contributions and is intended for small problems.

## Roadmap

1. Add a documented mesh input format and input validation/quality diagnostics.
2. Add volumetric heat-source integration and nonzero Neumann boundary loads.
3. Add richer material assignment, anisotropic conductivity, and transient heat conduction.
4. Improve solver diagnostics and add stronger preconditioners and larger-mesh benchmarks.
5. Expand verification with manufactured solutions and external reference comparisons.

## Repository contents

- `include/fem3d/fem3d.hpp` — public data structures and API.
- `src/fem3d.cpp` — element, sparse assembly, solver, and VTK implementation.
- `app/main.cpp` — sample command-line executable.
- `tests/test_fem3d.cpp` — dependency-free tests.
