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
    bool threw=false;
    try { (void)fem3d::tet4_element({{{0,0,0},{1,0,0},{2,0,0},{0,1,0}}}); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,"degenerate element must be rejected");
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
        test_tet_geometry(); test_invalid_tet(); test_cube_solve_and_vtk(); test_dirichlet_conflict();
        std::cout << "All fem3d tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failure: " << e.what() << '\n';
        return 1;
    }
}
