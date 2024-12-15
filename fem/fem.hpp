// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_FEM_HPP
#define MFEM_FEM_HPP

#include "bilinearform.hpp"
#include "bilininteg.hpp"
#include "ceed/solvers/algebraic.hpp"
#include "coefficient.hpp"
#include "complex_fem.hpp"
#include "convergence.hpp"
#include "datacollection.hpp"
#include "dgd/dgd_gridfunc.hpp"
#include "dgd/dgd_space.hpp"
#include "dgmassinv.hpp"
#include "doftrans.hpp"
#include "eltrans.hpp"
#include "estimators.hpp"
#include "fe.hpp"
#include "fe_coll.hpp"
#include "fespace.hpp"
#include "fespacehierarchy.hpp"
#include "geom.hpp"
#include "gridfunc.hpp"
#include "gslib.hpp"
#include "hybridization.hpp"
#include "hyperbolic.hpp"
#include "intrules.hpp"
#include "intrules_cut.hpp"
#include "kdtree.hpp"
#include "linearform.hpp"
#include "lininteg.hpp"
#include "lor/lor.hpp"
#include "multigrid.hpp"
#include "nonlinearform.hpp"
#include "nonlininteg.hpp"
#include "quadinterpolator.hpp"
#include "quadinterpolator_face.hpp"
#include "restriction.hpp"
#include "staticcond.hpp"
#include "tmop.hpp"
#include "tmop_amr.hpp"
#include "tmop_tools.hpp"
#include "transfer.hpp"

#ifdef MFEM_USE_MPI
#include "pbilinearform.hpp"
#include "pfespace.hpp"
#include "pgridfunc.hpp"
#include "plinearform.hpp"
#include "pnonlinearform.hpp"
#endif

#ifdef MFEM_USE_SIDRE
#include "sidredatacollection.hpp"
#endif

#ifdef MFEM_USE_CONDUIT
#include "conduitdatacollection.hpp"
#endif

#ifdef MFEM_USE_ADIOS2
#include "adios2datacollection.hpp"
#endif

#ifdef MFEM_USE_FMS
#include "fmsconvert.hpp"
#include "fmsdatacollection.hpp"
#endif

#endif
