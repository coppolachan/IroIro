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

#include "bfm.h"

class BFM_Storage {
  int Nx_, Ny_, Nz_, Nt_;
  GaugeField1D U_dag;
  bfm_dp& bfm_obj_;

  FermionField BasisConversion(FermionField& F,int Conversion); 
  BFM_Storage(); //Hide default Constructor
public:
  BFM_Storage(bfm_dp&); //The allowed constructor
  void GaugeExport_to_BFM(GaugeField& U);
  
  // Fermion routines need also Dirac basis tranformation
  // Dirac representation (IroIro) <-> Chiral representation (Bagel)

  void FermionExport_to_BFM(FermionField& F, Fermion_t handle, int cb);
  void FermionImport_from_BFM(FermionField& F, Fermion_t handle, int cb);


};


#endif
