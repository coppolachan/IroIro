//---------------------------------------------------------------------
/*! @file bfm_storage.hpp
  @brief Routines to convert IroIro fields to BFM (Bagel) type

  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef BFM_STORAGE_INCLUDED
#define BFM_STORAGE_INCLUDED

#include "include/common_fields.hpp"

#define DIRAC_TO_CHIRAL 0
#define CHIRAL_TO_DIRAC 1

//workaround
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "bfm.h"

//redefine some macros
#include "iroiro_config.h"

template <class Float>
class BFM_Storage{
  int Nx_, Ny_, Nz_, Nt_;
  int Nvol_;
  GaugeField1D U_dag;
  bfm_internal<Float>& bfm_obj_;

  void BasisConversion(FermionField& F_out, FermionField& F,int Conversion,int cb, int s = 0); 
  BFM_Storage(); //Hide default Constructor
public:
  BFM_Storage(bfm_internal<Float>&); //The allowed constructor
  void GaugeExport_to_BFM(GaugeField& U); // higher level
  void GaugeExport_to_BFM(const Field* U); // higher level
  void GaugeImport_from_BFM(Field* U, Matrix_t handle, int dir, int cb);

  // Fermion routines need also Dirac basis tranformation
  // Dirac representation (IroIro) <-> Chiral representation (Bagel)

  void FermionExport_to_BFM(FermionField& F, Fermion_t handle, int cb);
  void FermionImport_from_BFM(FermionField& F, Fermion_t handle, int cb);

  void FermionExport_to_BFM_5D(FermionField& F, Fermion_t handle, int cb, int s);
  void FermionImport_from_BFM_5D(FermionField& F, Fermion_t handle, int cb, int s, int Ls);


};

//#include "bfm_storage.cpp"

#endif
