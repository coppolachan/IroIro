//---------------------------------------------------------------------
/*! @file bfm_storage.cpp
  @brief Routines to convert IroIro fields to BFM (Bagel) type

  Class definitions
*/ 
//---------------------------------------------------------------------

#include "macros.hpp" //defines Color, Spin, Dimension constants
#include "commonPrms.h"

#include "bfm_storage.hpp"
#include "lib/Main/Geometry/shiftField.hpp"
#include "Tools/fieldUtils.hpp"

BFM_Storage::BFM_Storage():Nx_(CommonPrms::instance()->Nx()),
			   Ny_(CommonPrms::instance()->Ny()),
			   Nz_(CommonPrms::instance()->Nz()),
			   Nt_(CommonPrms::instance()->Nt()){};


void BFM_Storage::BasisConversion(FermionField& F, int ConvType){
  for(int site = 0; site<F.format.Nvol(); ++site){
    for(int col = 0; col<NC_; ++col){
      //write explicitly 
    }  
  }

}

double* BFM_Storage::GaugeExport_to_BFM(GaugeField& U){
  // BFM storage pattern N.B. : EVEN-ODD
  // double gauge [x%2][t][z][y][x/2] [row][column][realimag]  
  using namespace Mapping;
  using namespace FieldUtils;

  double* U_Ptr; 
  
  for (int mu = 0; mu < NDIM_; mu++) {
    GaugeField1D Umu = DirSlice(U, mu);
    int dir = 2*mu+1;
    U_Ptr = Umu.data.getaddr(0);
    //call some import function(U_ptr, dir)
    
    dir = 2*mu;
    U_dag = shiftField(Umu, mu, Backward()); // U(x-mu)
    for(int site = 0; site<U.format.Nvol(); ++site)
      U_dag.data[U_dag.format.cslice(0,site)] =
	mat_dag(U_dag,site).getva();
    
    U_Ptr = U_dag.data.getaddr(0);
    //call some import function(U_ptr, dir)
  }
};

double* BFM_Storage::FermionExport_to_BFM(double* U){
  // BFM storage pattern
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]
  double* F_ptr;

  //Convert Dirac basis: Dirac -> Chiral
  
  // Export




};

double* BFM_Storage::FermionImport_from_BFM(double* U, int check_board){
  

};

				
