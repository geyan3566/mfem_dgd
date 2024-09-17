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

  /// @brief
  void InitializeStencil();

  void GetBasisCenter(const int b_id, mfem::Vector &center,
                      const mfem::Vector &basisCenter) const;

  std::vector<std::size_t> sort_indexes(const std::vector<double> &v);

  /// Get the true number of dofs in the DGD space
  int GetTrueVSize() const override { return vdim * numBasis; }

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
  mutable mfem::Array<mfem::DenseMatrix *> coef;
};
} // namespace mfem

#endif