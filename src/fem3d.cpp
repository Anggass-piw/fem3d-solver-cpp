#include "fem3d/fem3d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>

namespace fem3d {
namespace {

Vec3 subtract(const Vec3& a, const Vec3& b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

struct Inverse3 {
    double a[3][3];
    double determinant;
};

Inverse3 inverse3(const Vec3& c1, const Vec3& c2, const Vec3& c3) {
    // Matrix columns are the three edges from vertex zero.
    const double a=c1.x, b=c2.x, c=c3.x;
    const double d=c1.y, e=c2.y, f=c3.y;
    const double g=c1.z, h=c2.z, i=c3.z;
    const double det=a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
    const double scale=std::max({std::abs(a),std::abs(b),std::abs(c),std::abs(d),std::abs(e),std::abs(f),std::abs(g),std::abs(h),std::abs(i)});
    if (!(scale>0.0) || !std::isfinite(det) || std::abs(det) <= 1e-14*scale*scale*scale)
        throw std::invalid_argument("degenerate or numerically singular tetrahedron");
    Inverse3 out{};
    out.determinant=det;
    out.a[0][0]=(e*i-f*h)/det; out.a[0][1]=(c*h-b*i)/det; out.a[0][2]=(b*f-c*e)/det;
    out.a[1][0]=(f*g-d*i)/det; out.a[1][1]=(a*i-c*g)/det; out.a[1][2]=(c*d-a*f)/det;
    out.a[2][0]=(d*h-e*g)/det; out.a[2][1]=(b*g-a*h)/det; out.a[2][2]=(a*e-b*d)/det;
    return out;
}

void validate_csr(const CSRMatrix& a) {
    if (a.row_offsets.size()!=a.rows+1 || a.column_indices.size()!=a.values.size() || a.row_offsets.back()!=a.values.size())
        throw std::invalid_argument("malformed CSR matrix");
}

} // namespace

ElementData tet4_element(const std::array<Vec3,4>& p, double conductivity) {
    if (!(conductivity>0.0) || !std::isfinite(conductivity)) throw std::invalid_argument("conductivity must be finite and positive");
    const auto inv=inverse3(subtract(p[1],p[0]),subtract(p[2],p[0]),subtract(p[3],p[0]));
    ElementData e;
    e.volume=std::abs(inv.determinant)/6.0;
    // Rows of B^-1 are gradients of reference coordinates, i.e. shape functions 1..3.
    e.shape_gradients[1]={inv.a[0][0],inv.a[0][1],inv.a[0][2]};
    e.shape_gradients[2]={inv.a[1][0],inv.a[1][1],inv.a[1][2]};
    e.shape_gradients[3]={inv.a[2][0],inv.a[2][1],inv.a[2][2]};
    e.shape_gradients[0]={-(e.shape_gradients[1].x+e.shape_gradients[2].x+e.shape_gradients[3].x),
                          -(e.shape_gradients[1].y+e.shape_gradients[2].y+e.shape_gradients[3].y),
                          -(e.shape_gradients[1].z+e.shape_gradients[2].z+e.shape_gradients[3].z)};
    for (std::size_t r=0;r<4;++r) for (std::size_t c=0;c<4;++c)
        e.stiffness[r][c]=conductivity*e.volume*dot(e.shape_gradients[r],e.shape_gradients[c]);
    return e;
}

CSRMatrix assemble(const Mesh& mesh, double conductivity) {
    if (mesh.points.empty()) throw std::invalid_argument("mesh has no points");
    std::vector<std::map<std::size_t,double>> rows(mesh.points.size());
    for (const auto& tet: mesh.tetrahedra) {
        std::array<Vec3,4> p;
        for (std::size_t i=0;i<4;++i) {
            if (tet.node[i]>=mesh.points.size()) throw std::out_of_range("tetrahedron references a missing point");
            p[i]=mesh.points[tet.node[i]];
        }
        const auto e=tet4_element(p,conductivity);
        for (std::size_t i=0;i<4;++i) for (std::size_t j=0;j<4;++j)
            rows[tet.node[i]][tet.node[j]]+=e.stiffness[i][j];
    }
    CSRMatrix a; a.rows=mesh.points.size(); a.row_offsets.reserve(a.rows+1); a.row_offsets.push_back(0);
    for (const auto& row: rows) {
        for (const auto& entry: row) { a.column_indices.push_back(entry.first); a.values.push_back(entry.second); }
        a.row_offsets.push_back(a.values.size());
    }
    return a;
}

std::vector<double> CSRMatrix::multiply(const std::vector<double>& x) const {
    validate_csr(*this);
    if (x.size()!=rows) throw std::invalid_argument("matrix/vector size mismatch");
    std::vector<double> y(rows,0.0);
    for (std::size_t r=0;r<rows;++r) for (std::size_t k=row_offsets[r];k<row_offsets[r+1];++k)
        y[r]+=values[k]*x[column_indices[k]];
    return y;
}

double CSRMatrix::diagonal(std::size_t row) const {
    validate_csr(*this);
    if (row>=rows) throw std::out_of_range("diagonal row out of range");
    for (std::size_t k=row_offsets[row];k<row_offsets[row+1];++k) if (column_indices[k]==row) return values[k];
    return 0.0;
}

void apply_dirichlet(const CSRMatrix& a, const std::vector<double>& load,
                     const std::vector<std::pair<std::size_t,double>>& prescribed,
                     CSRMatrix& out, std::vector<double>& rhs) {
    validate_csr(a);
    if (load.size()!=a.rows) throw std::invalid_argument("load/matrix size mismatch");
    std::vector<bool> fixed(a.rows,false); std::vector<double> values(a.rows,0.0);
    for (const auto& bc: prescribed) {
        if (bc.first>=a.rows) throw std::out_of_range("Dirichlet node out of range");
        if (!std::isfinite(bc.second)) throw std::invalid_argument("Dirichlet temperature must be finite");
        if (fixed[bc.first] && values[bc.first]!=bc.second) throw std::invalid_argument("conflicting Dirichlet values");
        fixed[bc.first]=true; values[bc.first]=bc.second;
    }
    out.rows=a.rows; out.row_offsets.clear(); out.column_indices.clear(); out.values.clear();
    out.row_offsets.push_back(0); rhs=load;
    for (std::size_t r=0;r<a.rows;++r) {
        if (fixed[r]) {
            out.column_indices.push_back(r); out.values.push_back(1.0); rhs[r]=values[r];
        } else {
            for (std::size_t k=a.row_offsets[r];k<a.row_offsets[r+1];++k) {
                const auto c=a.column_indices[k];
                if (fixed[c]) rhs[r]-=a.values[k]*values[c];
                else { out.column_indices.push_back(c); out.values.push_back(a.values[k]); }
            }
        }
        out.row_offsets.push_back(out.values.size());
    }
}

SolverResult solve_pcg(const CSRMatrix& a, const std::vector<double>& b, const SolverOptions& options) {
    validate_csr(a);
    if (b.size()!=a.rows) throw std::invalid_argument("rhs/matrix size mismatch");
    if (!(options.relative_tolerance>0.0) || !std::isfinite(options.relative_tolerance) || options.max_iterations==0)
        throw std::invalid_argument("invalid PCG options");
    SolverResult result; result.solution.assign(a.rows,0.0);
    if (a.rows==0) { result.converged=true; return result; }
    auto ax=a.multiply(result.solution); std::vector<double> r(b.size());
    double bnorm2=0.0;
    for (std::size_t i=0;i<b.size();++i) { r[i]=b[i]-ax[i]; bnorm2+=b[i]*b[i]; }
    const double target=options.relative_tolerance*std::sqrt(bnorm2);
    result.residual_norm=std::sqrt(bnorm2);
    if (result.residual_norm<=target) { result.converged=true; return result; }
    std::vector<double> z(a.rows),p(a.rows);
    double rz=0.0;
    for (std::size_t i=0;i<a.rows;++i) {
        const double d=a.diagonal(i);
        if (!(d>0.0) || !std::isfinite(d)) throw std::runtime_error("PCG requires a positive diagonal; check mesh and boundary conditions");
        z[i]=r[i]/d; p[i]=z[i]; rz+=r[i]*z[i];
    }
    for (std::size_t it=0;it<options.max_iterations;++it) {
        const auto ap=a.multiply(p); double pap=0.0;
        for (std::size_t i=0;i<a.rows;++i) pap+=p[i]*ap[i];
        if (!(pap>0.0) || !std::isfinite(pap)) throw std::runtime_error("PCG breakdown: matrix is not positive definite (unconstrained nullspace?)");
        const double alpha=rz/pap; double residual2=0.0;
        for (std::size_t i=0;i<a.rows;++i) { result.solution[i]+=alpha*p[i]; r[i]-=alpha*ap[i]; residual2+=r[i]*r[i]; }
        result.iterations=it+1; result.residual_norm=std::sqrt(residual2);
        if (result.residual_norm<=target) { result.converged=true; break; }
        double next_rz=0.0;
        for (std::size_t i=0;i<a.rows;++i) { z[i]=r[i]/a.diagonal(i); next_rz+=r[i]*z[i]; }
        const double beta=next_rz/rz;
        for (std::size_t i=0;i<a.rows;++i) p[i]=z[i]+beta*p[i];
        rz=next_rz;
    }
    return result;
}

void write_legacy_vtk(const std::string& path, const Mesh& mesh, const std::vector<double>& t) {
    if (t.size()!=mesh.points.size()) throw std::invalid_argument("temperature/mesh size mismatch");
    for (const auto& tet: mesh.tetrahedra) for (auto n:tet.node) if (n>=mesh.points.size()) throw std::out_of_range("tetrahedron references a missing point");
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open VTK output: "+path);
    out << "# vtk DataFile Version 3.0\nTet4 steady heat conduction\nASCII\nDATASET UNSTRUCTURED_GRID\n" << std::setprecision(17);
    out << "POINTS " << mesh.points.size() << " double\n";
    for (const auto& p:mesh.points) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
    out << "CELLS " << mesh.tetrahedra.size() << ' ' << mesh.tetrahedra.size()*5 << "\n";
    for (const auto& tet:mesh.tetrahedra) out << "4 " << tet.node[0] << ' ' << tet.node[1] << ' ' << tet.node[2] << ' ' << tet.node[3] << '\n';
    out << "CELL_TYPES " << mesh.tetrahedra.size() << "\n";
    for (std::size_t i=0;i<mesh.tetrahedra.size();++i) out << "10\n";
    out << "POINT_DATA " << t.size() << "\nSCALARS temperature double 1\nLOOKUP_TABLE default\n";
    for (double value:t) out << value << '\n';
    if (!out) throw std::runtime_error("failed while writing VTK output: "+path);
}

Mesh make_unit_cube_sample() {
    Mesh m;
    m.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1},{0.5,0.5,0.5}};
    // Each square boundary face is split into two triangles; center-to-face tetrahedra tile the cube.
    const std::array<std::array<std::size_t,4>,6> faces{{
        {{0,3,2,1}}, {{4,5,6,7}}, {{0,1,5,4}},
        {{3,7,6,2}}, {{0,4,7,3}}, {{1,2,6,5}}
    }};
    for (const auto& f:faces) {
        m.tetrahedra.push_back({{8,f[0],f[1],f[2]}});
        m.tetrahedra.push_back({{8,f[0],f[2],f[3]}});
    }
    return m;
}

} // namespace fem3d
