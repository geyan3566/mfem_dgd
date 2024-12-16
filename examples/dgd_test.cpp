#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

/// \brief Defines the exact solution for the steady isentropic vortex
/// \param[in] x - coordinate of the point at which the state is needed
/// \param[out] u - state variables stored as a 4-vector
void upoly(const Vector &x, Vector &u);
void utest(const Vector &x, Vector &u);
void usingle(const mfem::Vector &x, mfem::Vector &u);

/// Generate quarter annulus mesh
/// \param[in] degree - polynomial degree of the mapping
/// \param[in] num_rad - number of nodes in the radial direction
/// \param[in] num_ang - number of nodes in the angular direction
Mesh buildQuarterAnnulusMesh(int degree, int num_rad, int num_ang);

/// Build basis centers for the dgd space
/// \param[in] num_rad - number of nodes in the radial direction
/// \param[in] num_ang - number of nodes in the angular direction
mfem::Vector buildBasisCentersWithPerturbation(Mesh &mesh);

int main(int argc, char *argv[]) {
  int myid = 0;
  int degree = 1;
  int nx = 1;
  int ny = 1;
  int numRad = 10;
  int numTheta = 10;
  int extra = 1;
  // Parse command-line options
  OptionsParser args(argc, argv);
  args.AddOption(&degree, "-d", "--degree", "poly. degree of mesh mapping");
  args.AddOption(&nx, "-nr", "--num-rad", "number of radial segments");
  args.AddOption(&ny, "-nt", "--num-theta", "number of angular segments");
  args.AddOption(&numRad, "-br", "--numrad", "number of radius points");
  args.AddOption(&numTheta, "-bt", "--numtheta", "number of anglular points");
  args.AddOption(&extra, "-e", "--extra", "number of anglular points");
  args.ParseCheck();

  // working mesh
  Mesh smesh = buildQuarterAnnulusMesh(degree + 1, nx, ny);
  std::cout << "Number of elements " << smesh.GetNE() << '\n';
  int dim = smesh.Dimension();
  int num_state = dim + 2;

  // initialize the basis centers
  Vector center = buildBasisCentersWithPerturbation(smesh);

  // initialize the dgd space
  DG_FECollection fec(degree, dim, BasisType::GaussLobatto);
  DGDSpace dgdSpace(&smesh, &fec, center, degree, extra, num_state,
                    Ordering::byVDIM);
  FiniteElementSpace fes(&smesh, &fec, num_state, Ordering::byVDIM);

  // initialize the solution
  mfem::VectorFunctionCoefficient u0_fun(num_state, upoly);
  mfem::DgdGridFunction x_dgd(&dgdSpace);
  x_dgd.ProjectCoefficient(u0_fun);
  std::ofstream dgd_write("x_dgd.txt");
  x_dgd.Print(dgd_write, num_state);
  dgd_write.close();

  mfem::GridFunction x_exact(&fes);
  x_exact.ProjectCoefficient(u0_fun);
  std::ofstream fes_write("x_exact.txt");
  x_exact.Print(fes_write, num_state);
  fes_write.close();

  // prolong the solution
  mfem::GridFunction x_prolong(&fes);
  x_prolong = 0.0;
  dgdSpace.GetProlongationMatrix()->Mult(x_dgd, x_prolong);
  std::ofstream prolong_write("x_prolong.txt");
  x_prolong.Print(prolong_write, num_state);
  prolong_write.close();

  double max_diff = 0.0;
  int location = 0;
  for (int i = 0; i < x_prolong.Size(); ++i) {
    double diff = std::abs(x_prolong[i] - x_exact[i]);
    if (diff > max_diff) {
      max_diff = diff;
      location = i;
    }
  }

  cout << "max difference is " << max_diff << " location " << location << '\n';

  x_prolong -= x_exact;
  std::ofstream diff_write("x_diff.txt");
  x_prolong.Print(diff_write, num_state);
  diff_write.close();

  double l1 = x_prolong.Norml1();
  double l2 = x_prolong.Norml2();
  double linf = x_prolong.Normlinf();
  cout << "Check the projection l1 norm: " << l1 << '\n';
  cout << "Check the projection l2 error: " << l2 << '\n';
  cout << "Check the projection l2 error: " << linf << '\n';

  return 0;
}

Mesh buildQuarterAnnulusMesh(int degree, int num_rad, int num_ang) {
  auto mesh = Mesh::MakeCartesian2D(num_rad, num_ang, Element::TRIANGLE, true,
                                    2.0, M_PI * 0.5, true);
  // strategy:
  // 1) generate a fes for Lagrange elements of desired degree
  // 2) create a Grid Function using a VectorFunctionCoefficient
  // 4) use mesh_ptr->NewNodes(nodes, true) to set the mesh nodes

  // Problem: fes does not own fec, which is generated in this function's scope
  // Solution: the grid function can own both the fec and fes
  H1_FECollection *fec = new H1_FECollection(degree, 2 /* = dim */);
  FiniteElementSpace *fes =
      new FiniteElementSpace(&mesh, fec, 2, Ordering::byVDIM);

  // This lambda function transforms from (r,\theta) space to (x,y) space
  auto xy_fun = [](const Vector &rt, Vector &xy) {
    xy(0) =
        (rt(0) + 1.0) * cos(rt(1)); // need + 1.0 to shift r away from origin
    xy(1) = (rt(0) + 1.0) * sin(rt(1));
  };
  VectorFunctionCoefficient xy_coeff(2, xy_fun);
  GridFunction *xy = new GridFunction(fes);
  xy->MakeOwner(fec);
  xy->ProjectCoefficient(xy_coeff);

  mesh.NewNodes(*xy, true);
  return mesh;
}

void upoly(const mfem::Vector &x, mfem::Vector &u) {
  u.SetSize(4);
  u = 0.0;
  for (int p = 0; p < 4; p++) {
    u(p) = pow(x(0), p);
  }
}

void utest(const mfem::Vector &x, mfem::Vector &u) {
  u.SetSize(4);
  u(0) = 1.0;
  u(1) = 2.0;
  u(2) = 3.0;
  u(3) = 4.0;
}

void usingle(const mfem::Vector &x, mfem::Vector &u) {
  u.SetSize(1);
  u(0) = 2.0;
}

mfem::Vector buildBasisCentersWithPerturbation(Mesh &mesh) {
  // Setup random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(-1e-2, 1e-2);

  int numBasis = mesh.GetNE();
  Vector center(2 * numBasis);
  Vector loc(2);
  for (int k = 0; k < numBasis; k++) {
    mesh.GetElementCenter(k, loc);
    // Add random perturbation in [-0.01, 0.01]
    center(k * 2) = loc(0) + dist(gen);
    center(k * 2 + 1) = loc(1) + dist(gen);
  }
  return center;
}
