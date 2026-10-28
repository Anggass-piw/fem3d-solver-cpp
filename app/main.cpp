#include "fem3d/fem3d.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    std::string output="sample.vtk";
    for (int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if (arg=="--help" || arg=="-h") {
            std::cout << "Usage: fem3d_solver [--output FILE]\nSolves the unit-cube Tet4 steady heat-conduction demo (T=0 at x=0, T=1 at x=1).\n";
            return 0;
        }
        if (arg=="--output" && i+1<argc) { output=argv[++i]; continue; }
        std::cerr << "Unknown or incomplete argument: " << arg << "\nUse --help for usage.\n";
        return 2;
    }
    try {
        const auto mesh=fem3d::make_unit_cube_sample();
        const auto stiffness=fem3d::assemble(mesh,1.0);
        std::vector<double> load(mesh.points.size(),0.0);
        std::vector<std::pair<std::size_t,double>> bc;
        for (std::size_t i=0;i<mesh.points.size();++i) {
            if (mesh.points[i].x==0.0) bc.emplace_back(i,0.0);
            else if (mesh.points[i].x==1.0) bc.emplace_back(i,1.0);
        }
        fem3d::CSRMatrix constrained; std::vector<double> rhs;
        fem3d::apply_dirichlet(stiffness,load,bc,constrained,rhs);
        const auto result=fem3d::solve_pcg(constrained,rhs);
        if (!result.converged) {
            std::cerr << "PCG did not converge in " << result.iterations << " iterations (residual " << result.residual_norm << ")\n";
            return 1;
        }
        fem3d::write_legacy_vtk(output,mesh,result.solution);
        std::cout << "Solved " << mesh.points.size() << " nodes / " << mesh.tetrahedra.size()
                  << " tetrahedra in " << result.iterations << " PCG iterations; residual="
                  << result.residual_norm << "\nWrote " << output << "\nInterior center temperature: "
                  << result.solution[8] << " (expected 0.5)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
