#ifndef MFEM_CENTGRIDFUNC
#define MFEM_CENTGRIDFUNC

#include "gridfunc.hpp"

namespace mfem {
/// A derived grid function class used for store information on the element
/// center
class DGDGridFunction : public mfem::GridFunction {
public:
  DGDGridFunction() {}
  DGDGridFunction(mfem::FiniteElementSpace *f);
  DGDGridFunction(mfem::FiniteElementSpace *f, mfem::Vector center);

  virtual void ProjectCoefficient(mfem::VectorCoefficient &coeff);
  virtual void ProjectConstVecCoefficient(mfem::VectorCoefficient &coeff);

  DGDGridFunction &operator=(const Vector &v);
  DGDGridFunction &operator=(double value);

private:
  int dim;
  int numBasis;
  mfem::Vector basisCenter;
};

} // end of namespace mfem

#endif