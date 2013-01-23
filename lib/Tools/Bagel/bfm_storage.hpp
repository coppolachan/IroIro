//---------------------------------------------------------------------
/*! @file bfm_storage.hpp
  @brief Routines to convert IroIro fields to BFM (Bagel) type

  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef BFM_STORAGE_INCLUDED
#define BFM_STORAGE_INCLUDED

#include "include/common_fields.hpp"

class BFM_Storage {
  int Nx_, Ny_, Nz_, Nt_;
  GaugeField1D U_dag;

public:
  BFM_Storage(); //Constructor
  double* GaugeExport_to_BFM(GaugeField& U);
  
  // Fermion routines need also Dirac basis tranformation
  // Dirac representation (IroIro) <-> Chiral representation (Bagel)

  // Spinor basis change matrix
  // Matrix: 1  1
  //         1 -1
  // Need a 1/sqrt(2) on export
  double* FermionExport_to_BFM(double *F);
  double* FermionImport_from_BFM(double* F, int cb);


};


#endif
