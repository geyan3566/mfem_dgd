#ifndef MFEM_DGDGRIDFUNC
#define MFEM_DGDGRIDFUNC

#include "dgd_space.hpp"
#include "mfem.hpp"

namespace mfem {
/// A derived grid function class used for store information on the element
/// center
class DgdGridFunction : public mfem::GridFunction {
public:
  DgdGridFunction() {}
  DgdGridFunction(mfem::FiniteElementSpace *f);
  DgdGridFunction(mfem::FiniteElementSpace *f, mfem::Vector center);

  virtual void ProjectCoefficient(mfem::VectorCoefficient &coeff);
  virtual void ProjectConstVecCoefficient(mfem::VectorCoefficient &coeff);

  DgdGridFunction &operator=(const Vector &v);
  DgdGridFunction &operator=(double value);

private:
  int dim;
  int numBasis;
  mfem::Vector basisCenter;
};

} // end of namespace mfem

#endif