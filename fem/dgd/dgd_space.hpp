#ifndef MFEM_DGDSpace
#define MFEM_DGDSpace

#include "../../general/array.hpp"
#include "../../mesh/mesh.hpp"
#include "../fe_coll.hpp"
#include "../fespace.hpp"

namespace mfem {
class DGDSpace : public FiniteElementSpace {
public:
  /// class constructor
  DGDSpace(mfem::Mesh *m, const mfem::FiniteElementCollection *fec,
           mfem::Vector center, int degree, int extra, int vdim = 1,
           int ordering = mfem::Ordering::byVDIM);

  /// class destructor
  ~DGDSpace() = default;

  /// @brief initialize the dgd space given the basis center distribution
  void InitializeStencil();

  void GetBasisCenter(const int b_id, mfem::Vector &center,
                      const mfem::Vector &basisCenter) const;

  const Vector &GetBasisCenter() { return basisCenter; }

  void BuildProlongationMatrix(const mfem::Vector &centers);

  std::vector<std::size_t> sort_indexes(const std::vector<double> &v);

  /// Get the true number of dofs in the DGD space
  int GetTrueVSize() const override { return vdim * numBasis; }

  const std::vector<int> &GetSelectedBasis(int el_id) {
    return selectedBasis[el_id];
  }

  const std::vector<int> &GetSelectedElement(int b_id) {
    return selectedElement[b_id];
  }

protected:
  /// @brief Get data for constructing a dgd operator on an element
  void BuildElementDataMat(const int el_id, const mfem::Vector &basisCenter,
                           mfem::DenseMatrix &V, mfem::DenseMatrix &Vn) const;

  /// build the element-wise polynomial basis matrix
  void BuildElementPolyBasisMat(const int el_id,
                                const mfem::Vector &basisCenter,
                                const int numDofs,
                                const std::vector<mfem::Vector> &dofs_coord,
                                mfem::DenseMatrix &V,
                                mfem::DenseMatrix &Vn) const;

  /// @brief Solve local prolongation matrix
  void SolveLocalProlongationMat(const int el_id, const mfem::DenseMatrix &V,
                                 const mfem::DenseMatrix &Vn,
                                 mfem::DenseMatrix &localMat);

  /// @brief Assemble prolongation matrix
  void AssembleProlongationMatrix(const int el_id,
                                  const mfem::DenseMatrix &localMat);

protected:
  /// compute the derivative of prolongation matrix w.r.t the ith basis center
  void GetdPdc(const int i, const mfem::Vector &basisCenter,
               mfem::SparseMatrix &dpdc);

  void buildDerivDataMat(const int el_id, const int b_id, const int xyz,
                         const mfem::Vector &center, mfem::DenseMatrix &V,
                         mfem::DenseMatrix &dV, mfem::DenseMatrix &Vn) const;

  void AssembleDerivMatrix(const int el_id, const DenseMatrix &dpdc_block,
                           mfem::SparseMatrix &dpdc) const;

protected:
  /// number of radial basis function
  int numBasis;
  /// number of polynomial basis
  int numPolyBasis;
  /// polynomial order
  int polyOrder;
  /// number of required basis to constructe certain order polynomial
  int numLocalBasis;
  /// number of extra basis
  int extra;

  /// basis center
  mfem::Vector basisCenter;

  /// selected basis for each element
  std::vector<std::vector<int>> selectedBasis;

  /// selected element for each basis
  std::vector<std::vector<int>> selectedElement;

  /// element basis centers distances
  std::vector<std::vector<double>> elementBasisDist;

  /// @brief coefficient matrix for radial basis function
  std::vector<mfem::DenseMatrix> coef;
};
} // namespace mfem

#endif