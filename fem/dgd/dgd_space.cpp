// DGD space class definition

#include "dgd_space.hpp"

#include <algorithm>
#include <memory>
#include <numeric>

namespace mfem {

DGDSpace::DGDSpace(mfem::Mesh *m, const mfem::FiniteElementCollection *fec,
                   mfem::Vector center, int degree, int extra, int vdim,
                   int ordering)
    : mfem::FiniteElementSpace(m, fec, vdim, ordering) {
  int dim = m->Dimension();
  numBasis = center.Size() / dim;
  switch (dim) {
  case 1:
    numPolyBasis = polyOrder + 1;
    break;
  case 2:
    numPolyBasis = (polyOrder + 1) * (polyOrder + 2) / 2;
    break;
  case 3:
    numPolyBasis = (polyOrder + 1) * (polyOrder + 2) * (polyOrder + 3) / 6;
    break;
  }
  numLocalBasis = numPolyBasis + extra;
  out << "Number of total basis center is " << center.Size() / dim << '\n';
  out << "Number of required polynomial basis is " << numPolyBasis << '\n';
  out << "Number of element local basis is " << numLocalBasis << '\n';

  // initialize the stencil/patch
  InitializeStencil();

  // build the initial prolongation matrix
  cP = std::make_unique<SparseMatrix>(GetVSize(), GetTrueVSize());
  BuildProlongationMatrix(center);
}

void DGDSpace::InitializeStencil() {
  int i, j, k;
  int num_el = mesh->GetNE();
  int dim = mesh->Dimension();

  // initialize the all element centers for later used
  elementBasisDist.assign(num_el, {});
  selectedBasis.assign(num_el, {});
  selectedElement.assign(numBasis, {});
  coef.SetSize(GetMesh()->GetNE());

  // aux data
  double dist;
  Vector elemCenter(dim);
  Vector center(dim);
  // loop over all the elements to construct the stencil
  std::vector<size_t> temp;
  for (i = 0; i < num_el; i++) {
    coef[i] = new DenseMatrix(numPolyBasis, numLocalBasis);
    mesh->GetElementCenter(i, elemCenter);
    // loop over all basis
    for (j = 0; j < numBasis; j++) {
      GetBasisCenter(j, center, basisCenter);
      center -= elemCenter;
      dist = center.Norml2();
      elementBasisDist[i].push_back(dist);
    }
    // build element/basis stencil based on distance
    temp = sort_indexes(elementBasisDist[i]);
    for (k = 0; k < numLocalBasis; k++) {
      selectedBasis[i].push_back(temp[k]);
      selectedElement[temp[k]].push_back(i);
    }
  }

  out << "------Check the stencil------\n";
  out << "------Basis center loca------\n";
  for (int i = 0; i < numBasis; i++) {
    out << "basis " << i << ": ";
    out << basisCenter(2 * i) << ' ' << basisCenter(2 * i + 1) << '\n';
  }
  out << '\n';
  out << "------Elem's  stencil------\n";
  for (int i = 0; i < num_el; i++) {
    out << "Element " << i << ": ";
    for (int j = 0; j < selectedBasis[i].size(); j++) {
      out << selectedBasis[i][j] << ' ';
    }
    out << '\n';
  }
  out << '\n';
  out << "------Basis's  element------\n";
  for (int k = 0; k < numBasis; k++) {
    out << "basis " << k << ": ";
    for (int l = 0; l < selectedElement[k].size(); l++) {
      out << selectedElement[k][l] << ' ';
    }
    out << '\n';
  }
}

void DGDSpace::GetBasisCenter(const int b_id, Vector &center,
                              const Vector &basisCenter) const {
  int dim = basisCenter.Size();
  for (int i = 0; i < dim; i++) {
    center(i) = basisCenter(b_id * dim + i);
  }
}

std::vector<size_t> DGDSpace::sort_indexes(const std::vector<double> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

void DGDSpace::BuildProlongationMatrix(const Vector &x) {
  DenseMatrix V, Vn;
  DenseMatrix localMat;
  int num_el = mesh->GetNE();
  for (int i = 0; i < num_el; ++i) {
    // 1. build basis matrix
    BuildElementDataMat(i, x, V, Vn);

    // 2. build the interpolation matrix
    SolveLocalProlongationMat(i, V, Vn, localMat);

    // 3. Assemble prolongation matrix
    AssembleProlongationMatrix(i, localMat);
  }
  cP->Finalize();
  // ofstream cp_save("prolong.txt");
  // cP->PrintMatlab(cp_save);
  // cp_save.close();
}

void DGDSpace::BuildElementDataMat(int el_id, const Vector &x, DenseMatrix &V,
                                   DenseMatrix &Vn) const {
  // get element related data
  const Element *el = mesh->GetElement(el_id);
  const FiniteElement *fe =
      fec->FiniteElementForGeometry(el->GetGeometryType());
  const int numDofs = fe->GetDof();
  ElementTransformation *eltransf = mesh->GetElementTransformation(el_id);

  // get the dofs coord
  Array<Vector> dofs_coord;
  // std::vector<std::vector<double>> dofs_coord;
  dofs_coord.SetSize(numDofs);
  int dim = mesh->Dimension();
  Vector coord(dim);
  for (int k = 0; k < numDofs; k++) {
    dofs_coord[k].SetSize(dim);
    eltransf->Transform(fe->GetNodes().IntPoint(k), dofs_coord[k]);
  }

  V.SetSize(numLocalBasis, numPolyBasis);
  Vn.SetSize(numDofs, numPolyBasis);
  // build the data matrix
  buildElementPolyBasisMat(el_id, x, numDofs, dofs_coord, V, Vn);

  // free the aux variable
  // for (int k = 0; k < numDofs; k++)
  // {
  //    delete dofs_coord[k];
  // }
}

} // namespace mfem