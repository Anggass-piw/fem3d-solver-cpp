#include "fem3d/fem3d.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string& message) { if (!condition) throw std::runtime_error(message); }
void near(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual-expected)>tolerance) throw std::runtime_error(message+": got "+std::to_string(actual)+", expected "+std::to_string(expected));
}
template <typename Function>
void require_invalid_argument(Function&& function, const std::string& message) {
    bool threw=false;
    try { function(); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,message);
}
void test_tet_geometry() {
    const std::array<fem3d::Vec3,4> p{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    const auto e=fem3d::tet4_element(p,2.0);
    near(e.volume,1.0/6.0,1e-14,"unit right tetra volume");
    near(e.shape_gradients[0].x,-1,1e-14,"N0 gradient x");
    near(e.shape_gradients[1].x,1,1e-14,"N1 gradient x");
    for (int i=0;i<4;++i) {
        double row=0;
        for (int j=0;j<4;++j) { near(e.stiffness[i][j],e.stiffness[j][i],1e-14,"stiffness symmetry"); row+=e.stiffness[i][j]; }
        near(row,0,1e-14,"stiffness row sum");
    }
    const std::array<fem3d::Vec3,4> reversed{{p[0],p[2],p[1],p[3]}};
    near(fem3d::tet4_element(reversed).volume,1.0/6.0,1e-14,"orientation-independent volume");
}
void test_invalid_tet() {
    require_invalid_argument(
        [] { (void)fem3d::tet4_element({{{0,0,0},{1,0,0},{2,0,0},{0,1,0}}}); },
        "degenerate element must be rejected");
}
void test_csr_invariants() {
    const auto expect_rejected_by_all_operations=[](const fem3d::CSRMatrix& malformed) {
        require_invalid_argument(
            [&] { (void)malformed.multiply(std::vector<double>(malformed.rows,0.0)); },
            "multiply rejects malformed CSR");
        require_invalid_argument(
            [&] { (void)malformed.diagonal(0); },
            "diagonal rejects malformed CSR");
        require_invalid_argument(
            [&] { (void)fem3d::solve_pcg(malformed,std::vector<double>(malformed.rows,0.0)); },
            "PCG rejects malformed CSR");
        fem3d::CSRMatrix output; std::vector<double> rhs;
        require_invalid_argument(
            [&] { fem3d::apply_dirichlet(malformed,std::vector<double>(malformed.rows,0.0),{},output,rhs); },
            "Dirichlet elimination rejects malformed CSR");
    };

    fem3d::CSRMatrix nonmonotone_offsets;
    nonmonotone_offsets.rows=3;
    nonmonotone_offsets.row_offsets={0,2,1,2};
    nonmonotone_offsets.column_indices={0,1};
    nonmonotone_offsets.values={1.0,1.0};
    expect_rejected_by_all_operations(nonmonotone_offsets);

    fem3d::CSRMatrix nonzero_start;
    nonzero_start.rows=1;
    nonzero_start.row_offsets={1,1};
    nonzero_start.column_indices={0};
    nonzero_start.values={1.0};
    expect_rejected_by_all_operations(nonzero_start);

    fem3d::CSRMatrix bad_column;
    bad_column.rows=1;
    bad_column.row_offsets={0,1};
    bad_column.column_indices={1};
    bad_column.values={1.0};
    expect_rejected_by_all_operations(bad_column);

    fem3d::CSRMatrix missing_offset;
    missing_offset.rows=1;
    missing_offset.row_offsets={0};
    expect_rejected_by_all_operations(missing_offset);

    fem3d::CSRMatrix mismatched_entries;
    mismatched_entries.rows=1;
    mismatched_entries.row_offsets={0,1};
    mismatched_entries.column_indices={0};
    expect_rejected_by_all_operations(mismatched_entries);

    fem3d::CSRMatrix bad_terminal_offset;
    bad_terminal_offset.rows=1;
    bad_terminal_offset.row_offsets={0,0};
    bad_terminal_offset.column_indices={0};
    bad_terminal_offset.values={1.0};
    expect_rejected_by_all_operations(bad_terminal_offset);
}
void test_apply_dirichlet_aliasing() {
    fem3d::CSRMatrix matrix;
    matrix.rows=2;
    matrix.row_offsets={0,2,4};
    matrix.column_indices={0,1,0,1};
    matrix.values={4.0,1.0,1.0,3.0};
    std::vector<double> load={5.0,5.0};

    // Exercise both supported aliases: matrix == constrained_matrix and load == constrained_load.
    fem3d::apply_dirichlet(matrix,load,{{0,2.0}},matrix,load);

    require(matrix.rows==2,"in-place constrained matrix row count");
    require(matrix.row_offsets==std::vector<std::size_t>({0,1,2}),"in-place constrained matrix offsets");
    require(matrix.column_indices==std::vector<std::size_t>({0,1}),"in-place constrained matrix columns");
    require(matrix.values==std::vector<double>({1.0,3.0}),"in-place constrained matrix values");
    require(load==std::vector<double>({2.0,3.0}),"in-place constrained load values");
}
void test_cube_solve_and_vtk() {
    const auto mesh=fem3d::make_unit_cube_sample();
    require(mesh.points.size()==9 && mesh.tetrahedra.size()==12,"sample mesh dimensions");
    const auto a=fem3d::assemble(mesh);
    std::vector<double> ones(mesh.points.size(),1.0);
    const auto product=a.multiply(ones);
    for (double value:product) near(value,0.0,1e-12,"assembled stiffness annihilates constants");
    std::vector<std::pair<std::size_t,double>> bc;
    for (std::size_t i=0;i<mesh.points.size();++i) {
        if (mesh.points[i].x==0) bc.emplace_back(i,0.0);
        else if (mesh.points[i].x==1) bc.emplace_back(i,1.0);
    }
    fem3d::CSRMatrix constrained; std::vector<double> rhs;
    fem3d::apply_dirichlet(a,std::vector<double>(mesh.points.size(),0.0),bc,constrained,rhs);
    const auto result=fem3d::solve_pcg(constrained,rhs);
    require(result.converged,"PCG converged");
    near(result.solution[8],0.5,1e-12,"linear-field center temperature");
    for (std::size_t i=0;i<mesh.points.size();++i) near(result.solution[i],mesh.points[i].x,1e-12,"linear field nodal temperature");
    const std::string path="fem3d_test_output.vtk";
    fem3d::write_legacy_vtk(path,mesh,result.solution);
    std::ifstream file(path); require(static_cast<bool>(file),"VTK file created");
    std::string content((std::istreambuf_iterator<char>(file)),{});
    require(content.find("DATASET UNSTRUCTURED_GRID")!=std::string::npos,"VTK dataset header");
    require(content.find("SCALARS temperature double 1")!=std::string::npos,"VTK temperature field");
    std::remove(path.c_str());
}
void test_dirichlet_conflict() {
    fem3d::CSRMatrix a; a.rows=1; a.row_offsets={0,1}; a.column_indices={0}; a.values={1};
    fem3d::CSRMatrix out; std::vector<double> rhs; bool threw=false;
    try { fem3d::apply_dirichlet(a,{0},{{0,1},{0,2}},out,rhs); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,"conflicting boundary values rejected");
}
}
int main() {
    try {
        test_tet_geometry(); test_invalid_tet(); test_csr_invariants();
        test_apply_dirichlet_aliasing(); test_cube_solve_and_vtk(); test_dirichlet_conflict();
        std::cout << "All fem3d tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failure: " << e.what() << '\n';
        return 1;
    }
}
