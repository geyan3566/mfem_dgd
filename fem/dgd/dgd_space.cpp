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
  cP.reset(new SparseMatrix(GetVSize(), GetTrueVSize()));
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
  coef.resize(mesh->GetNE());

  // aux data
  double dist;
  Vector elemCenter(dim);
  Vector center(dim);
  // loop over all the elements to construct the stencil
  std::vector<size_t> temp;
  for (i = 0; i < num_el; i++) {
    coef[i] = DenseMatrix(numPolyBasis, numLocalBasis);
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
    // 1. build basis
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
  std::vector<Vector> dofs_coord;
  // std::vector<std::vector<double>> dofs_coord;
  dofs_coord.reserve(numDofs);
  int dim = mesh->Dimension();
  Vector coord(dim);
  for (int k = 0; k < numDofs; k++) {
    dofs_coord[k].SetSize(dim);
    eltransf->Transform(fe->GetNodes().IntPoint(k), dofs_coord[k]);
  }

  V.SetSize(numLocalBasis, numPolyBasis);
  Vn.SetSize(numDofs, numPolyBasis);
  // build the data matrix
  BuildElementPolyBasisMat(el_id, x, numDofs, dofs_coord, V, Vn);
}

void DGDSpace::BuildElementPolyBasisMat(const int el_id,
                                        const Vector &basisCenter,
                                        const int numDofs,
                                        const std::vector<Vector> &dofs_coord,
                                        DenseMatrix &V, DenseMatrix &Vn) const {
  int i, j, k, l;
  double dx, dy, dz;
  int b_id;
  int dim = mesh->Dimension();
  Vector loc_coord(dim);
  Vector el_center(dim);
  mesh->GetElementCenter(el_id, el_center);
  if (1 == dim) {
    // form the V matrix
    for (i = 0; i < numLocalBasis; i++) {
      b_id = selectedBasis[el_id][i];
      GetBasisCenter(b_id, loc_coord, basisCenter);
      dx = loc_coord[0] - el_center[0];
      for (j = 0; j <= polyOrder; j++) {
        V(i, j) = pow(dx, j);
      }
    }

    // form the Vn matrix
    for (i = 0; i < numDofs; i++) {
      loc_coord = dofs_coord[i];
      dx = loc_coord[0] - el_center[0];
      for (j = 0; j <= polyOrder; j++) {
        Vn(i, j) = pow(dx, j);
      }
    }
  } else if (2 == dim) {
    // form the V matrix
    for (i = 0; i < numLocalBasis; i++) {
      b_id = selectedBasis[el_id][i];
      GetBasisCenter(b_id, loc_coord, basisCenter);
      dx = loc_coord[0] - el_center[0];
      dy = loc_coord[1] - el_center[1];
      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          V(i, col) = pow(dx, j - k) * pow(dy, k);
          col++;
        }
      }
    }
    // form the Vn matrix
    for (i = 0; i < numDofs; i++) {
      loc_coord = dofs_coord[i];
      dx = loc_coord[0] - el_center[0];
      dy = loc_coord[1] - el_center[1];
      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          Vn(i, col) = pow(dx, j - k) * pow(dy, k);
          col++;
        }
      }
    }
  } else if (3 == dim) {
    // form the V matrix
    for (i = 0; i < numLocalBasis; i++) {
      b_id = selectedBasis[el_id][i];
      GetBasisCenter(b_id, loc_coord, basisCenter);
      dx = loc_coord[0] - el_center[0];
      dy = loc_coord[1] - el_center[1];
      dz = loc_coord[2] - el_center[2];
      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          for (l = 0; l <= j - k; l++) {
            V(i, col) = pow(dx, j - k - l) * pow(dy, l) * pow(dz, k);
            col++;
          }
        }
      }
    }

    // form the Vn matrix
    for (i = 0; i < numDofs; i++) {
      loc_coord = dofs_coord[i];
      dx = loc_coord[0] - el_center[0];
      dy = loc_coord[1] - el_center[1];
      dz = loc_coord[2] - el_center[2];
      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          for (l = 0; l <= j - k; l++) {
            Vn(i, col) = pow(dx, j - k - l) * pow(dy, l) * pow(dz, k);
            col++;
          }
        }
      }
    }
  }
}

void DGDSpace::SolveLocalProlongationMat(const int el_id, const DenseMatrix &V,
                                         const DenseMatrix &Vn,
                                         DenseMatrix &localMat) {
  int numDofs = Vn.Height();
  DenseMatrix b(numLocalBasis, numLocalBasis);
  b = 0.0;
  for (int i = 0; i < numLocalBasis; i++) {
    b(i, i) = 1.0;
  }

  if (numPolyBasis == numLocalBasis) {
    DenseMatrixInverse Vinv(V);
    Vinv.Mult(b, coef[el_id]);
  } else {
    DenseMatrix Vt(V);
    Vt.Transpose();
    DenseMatrix VtV(numPolyBasis, numPolyBasis);
    Mult(Vt, V, VtV);

    DenseMatrixInverse Vinv(VtV);
    DenseMatrix Vtb(numPolyBasis, numLocalBasis);
    Mult(Vt, b, Vtb);
    Vinv.Mult(Vtb, coef[el_id]);
  }
  localMat.SetSize(numDofs, numLocalBasis);
  Mult(Vn, coef[el_id], localMat);
}

void DGDSpace::AssembleProlongationMatrix(const int el_id,
                                          const DenseMatrix &localMat) {
  // element id coresponds to the column indices
  // dofs id coresponds to the row indices
  // the local reconstruction matrix needs to be assembled `vdim` times
  // assume the mesh only contains only 1 type of element
  const Element *el = mesh->GetElement(el_id);
  const FiniteElement *fe =
      fec->FiniteElementForGeometry(el->GetGeometryType());
  const int numDofs = fe->GetDof();

  int numLocalBasis = selectedBasis[el_id].size();
  Array<int> el_dofs;
  Array<int> col_index(numLocalBasis);
  Array<int> row_index(numDofs);

  GetElementVDofs(el_id, el_dofs);
  for (int e = 0; e < numLocalBasis; e++) {
    col_index[e] = vdim * selectedBasis[el_id][e];
  }
  for (int v = 0; v < vdim; v++) {
    el_dofs.GetSubArray(v * numDofs, numDofs, row_index);
    cP->SetSubMatrix(row_index, col_index, localMat, 1);
    // row_index.LoseData();
    //  elements id also need to be shift accordingly
    for (int e = 0; e < numLocalBasis; e++) {
      col_index[e]++;
    }
  }
}

void DGDSpace::GetdPdc(const int id, const Vector &basisCenter,
                       SparseMatrix &dpdc) {
  int dim = mesh->Dimension();
  int xyz = id % dim;  // determine whether it is x, y, or z
  int b_id = id / dim; // determine the basis id
  int numLocalElem = selectedElement[b_id].size();
  int el_id;
  DenseMatrix V;
  DenseMatrix dV;
  DenseMatrix Vn;
  DenseMatrix dpdc_block;
  for (int i = 0; i < numLocalElem; i++) {
    el_id = selectedElement[b_id][i];
    buildDerivDataMat(el_id, b_id, xyz, basisCenter, V, dV, Vn);
    dpdc_block.SetSize(Vn.Height(), numLocalBasis);

    // V is a square matrix
    if (numPolyBasis == numLocalBasis) {
      DenseMatrix temp_mat1(numLocalBasis);
      DenseMatrix temp_mat2(numLocalBasis);
      Mult(dV, coef[el_id], temp_mat1);
      Mult(coef[el_id], temp_mat1, temp_mat2);
      temp_mat2.Neg();                 // -V^-1 * dV * V^-1
      Mult(Vn, temp_mat2, dpdc_block); //  dpdc = Vn * temp2
    }
    // V is overdetermined
    else {
      DenseMatrix Vt(V);
      Vt.Transpose(); // get V^t

      DenseMatrix vtv(numPolyBasis);
      Mult(Vt, V, vtv);
      DenseMatrixInverse vtvinv(vtv); // get (V^t V)^-1

      DenseMatrix dVt(dV);
      dVt.Transpose(); // get dV^t

      DenseMatrix dvtv(numPolyBasis);
      Mult(Vt, dV, dvtv);
      AddMult(dVt, V, dvtv); // compute d V^tV / dc

      DenseMatrix temp_mat1(numPolyBasis);
      DenseMatrix deriv_p1(numPolyBasis, numLocalBasis);
      vtvinv.Mult(dvtv, temp_mat1);
      Mult(temp_mat1, coef[el_id], deriv_p1);
      deriv_p1.Neg(); // first part of the derivatve

      DenseMatrix deriv_p2(numPolyBasis, numLocalBasis);
      vtvinv.Mult(dVt, deriv_p2);

      deriv_p1 += deriv_p2;
      Mult(Vn, deriv_p1, dpdc_block);
    }
    // assemble is back to the derivative matrix
    AssembleDerivMatrix(el_id, dpdc_block, dpdc);
  }
  dpdc.Finalize();
}

void DGDSpace::buildDerivDataMat(const int el_id, const int b_id, const int xyz,
                                 const mfem::Vector &center,
                                 mfem::DenseMatrix &V, mfem::DenseMatrix &dV,
                                 mfem::DenseMatrix &Vn) const {
  int i, j, k, l;
  double dx, dy, dz;
  int dim = mesh->Dimension();
  Vector loc_coord(dim);
  Vector el_center(dim);
  mesh->GetElementCenter(el_id, el_center);

  // Get element related data
  const Element *el = mesh->GetElement(el_id);
  const FiniteElement *fe =
      fec->FiniteElementForGeometry(el->GetGeometryType());
  const int numDofs = fe->GetDof();

  // Initialize matrices
  V.SetSize(numLocalBasis, numPolyBasis);
  dV.SetSize(numLocalBasis, numPolyBasis);
  Vn.SetSize(numDofs, numPolyBasis);

  if (1 == dim) {
    // form the V and dV matrices
    for (i = 0; i < numLocalBasis; i++) {
      int basis_id = selectedBasis[el_id][i];
      GetBasisCenter(basis_id, loc_coord, center);
      dx = loc_coord[0] - el_center[0];

      for (j = 0; j <= polyOrder; j++) {
        V(i, j) = pow(dx, j);
        dV(i, j) = (j == 0)
                       ? 0.0
                       : j * pow(dx, j - 1) * (basis_id == b_id ? 1.0 : 0.0);
      }
    }

    // form the Vn matrix (unchanged from BuildElementPolyBasisMat)
    for (i = 0; i < numDofs; i++) {
      ElementTransformation *T = mesh->GetElementTransformation(el_id);
      Vector dof_coord(dim);
      T->Transform(fe->GetNodes().IntPoint(i), dof_coord);
      dx = dof_coord[0] - el_center[0];

      for (j = 0; j <= polyOrder; j++) {
        Vn(i, j) = pow(dx, j);
      }
    }
  } else if (2 == dim) {
    for (i = 0; i < numLocalBasis; i++) {
      int basis_id = selectedBasis[el_id][i];
      GetBasisCenter(basis_id, loc_coord, center);
      dx = loc_coord[0] - el_center[0];
      dy = loc_coord[1] - el_center[1];

      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          V(i, col) = pow(dx, j - k) * pow(dy, k);
          if (basis_id == b_id) {
            if (xyz == 0) { // x derivative
              dV(i, col) = (j - k == 0)
                               ? 0.0
                               : (j - k) * pow(dx, j - k - 1) * pow(dy, k);
            } else { // y derivative
              dV(i, col) = (k == 0) ? 0.0 : k * pow(dx, j - k) * pow(dy, k - 1);
            }
          } else {
            dV(i, col) = 0.0;
          }
          col++;
        }
      }
    }

    // form the Vn matrix
    for (i = 0; i < numDofs; i++) {
      ElementTransformation *T = mesh->GetElementTransformation(el_id);
      Vector dof_coord(dim);
      T->Transform(fe->GetNodes().IntPoint(i), dof_coord);
      dx = dof_coord[0] - el_center[0];
      dy = dof_coord[1] - el_center[1];

      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          Vn(i, col) = pow(dx, j - k) * pow(dy, k);
          col++;
        }
      }
    }
  } else if (3 == dim) {
    for (i = 0; i < numLocalBasis; i++) {
      int basis_id = selectedBasis[el_id][i];
      GetBasisCenter(basis_id, loc_coord, center);
      dx = loc_coord[0] - el_center[0];
      dy = loc_coord[1] - el_center[1];
      dz = loc_coord[2] - el_center[2];

      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          for (l = 0; l <= j - k; l++) {
            V(i, col) = pow(dx, j - k - l) * pow(dy, l) * pow(dz, k);
            if (basis_id == b_id) {
              if (xyz == 0) { // x derivative
                dV(i, col) = (j - k - l == 0)
                                 ? 0.0
                                 : (j - k - l) * pow(dx, j - k - l - 1) *
                                       pow(dy, l) * pow(dz, k);
              } else if (xyz == 1) { // y derivative
                dV(i, col) = (l == 0) ? 0.0
                                      : l * pow(dx, j - k - l) *
                                            pow(dy, l - 1) * pow(dz, k);
              } else { // z derivative
                dV(i, col) = (k == 0) ? 0.0
                                      : k * pow(dx, j - k - l) * pow(dy, l) *
                                            pow(dz, k - 1);
              }
            } else {
              dV(i, col) = 0.0;
            }
            col++;
          }
        }
      }
    }

    // form the Vn matrix
    for (i = 0; i < numDofs; i++) {
      ElementTransformation *T = mesh->GetElementTransformation(el_id);
      Vector dof_coord(dim);
      T->Transform(fe->GetNodes().IntPoint(i), dof_coord);
      dx = dof_coord[0] - el_center[0];
      dy = dof_coord[1] - el_center[1];
      dz = dof_coord[2] - el_center[2];

      int col = 0;
      for (j = 0; j <= polyOrder; j++) {
        for (k = 0; k <= j; k++) {
          for (l = 0; l <= j - k; l++) {
            Vn(i, col) = pow(dx, j - k - l) * pow(dy, l) * pow(dz, k);
            col++;
          }
        }
      }
    }
  }
}

void DGDSpace::AssembleDerivMatrix(const int el_id, const DenseMatrix &localMat,
                                   SparseMatrix &dpdc) const {
  // element id coresponds to the column indices
  // dofs id coresponds to the row indices
  // the local reconstruction matrix needs to be assembled `vdim` times
  // assume the mesh only contains only 1 type of element
  const Element *el = mesh->GetElement(el_id);
  const FiniteElement *fe =
      fec->FiniteElementForGeometry(el->GetGeometryType());
  const int numDofs = fe->GetDof();

  int numLocalBasis = selectedBasis[el_id].size();
  Array<int> el_dofs;
  Array<int> col_index(numLocalBasis);
  Array<int> row_index(numDofs);

  GetElementVDofs(el_id, el_dofs);
  for (int e = 0; e < numLocalBasis; e++) {
    col_index[e] = vdim * selectedBasis[el_id][e];
  }

  for (int v = 0; v < vdim; v++) {
    el_dofs.GetSubArray(v * numDofs, numDofs, row_index);
    dpdc.SetSubMatrix(row_index, col_index, localMat, 1);
    // elements id also need to be shift accordingly
    for (int e = 0; e < numLocalBasis; e++) {
      col_index[e]++;
    }
  }
}

} // namespace mfem