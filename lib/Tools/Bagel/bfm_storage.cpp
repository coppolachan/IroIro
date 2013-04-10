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



BFM_Storage::BFM_Storage(bfm_dp& bfm):Nx_(CommonPrms::instance()->Nx()),
			   Ny_(CommonPrms::instance()->Ny()),
			   Nz_(CommonPrms::instance()->Nz()),
			   Nt_(CommonPrms::instance()->Nt()),
			   bfm_obj_(bfm){};


FermionField BFM_Storage::BasisConversion(FermionField& F, int ConvType){
  //sqrt(2)^-1 = 0.707106781 = a
  // Gamma matrices transformation (\Gamma)
  //  0   1   0  -1
  // -1   0   1   0     
  //  0   1   0   1
  // -1   0  -1   0
  // times "a" defined above to make the determinant 1
  // In compact form is
  //
  // A  -A
  // A   A
  // where A is a 2x2 matrix 
  //  0 1
  // -1 0
  //
  // The inverse (\Gamma^-1) is
  //
  // -A  -A
  //  A  -A
  // times "a" 

  double factor = 1.0;
  if (ConvType == CHIRAL_TO_DIRAC)
    factor = -1.0;

  FermionField F_out(F.Nvol());
  for(int site = 0; site<F.format.Nvol(); ++site){
    for(int c = 0; c<NC_; ++c){
      //write explicitly the spinor indexes multiplication
 
     // \psi_0 = a*(\psi_1 - \psi_3)
      F_out.data.set(F.format.index_r(c, 0, site) , 
		     0.707106781*(factor*F[F.format.index_r(c,1,site)]
				  - F[F.format.index_r(c,3,site)]));
      F_out.data.set(F.format.index_i(c, 0, site) , 
		     0.707106781*(factor*F[F.format.index_i(c,1,site)] 
				  - F[F.format.index_i(c,3,site)]));

      // \psi_1 = a*(\psi_2 - \psi_0)
      F_out.data.set(F.format.index_r(c, 1, site) ,
		     0.707106781*(F[F.format.index_r(c,2,site)] 
				  - factor*F[F.format.index_r(c,0,site)]));
      F_out.data.set(F.format.index_i(c, 1, site),  
		     0.707106781*(F[F.format.index_i(c,2,site)]
				  - factor*F[F.format.index_i(c,0,site)]));

      // \psi_2 = a*(\psi_1 + \psi_3)
      F_out.data.set(F.format.index_r(c, 2, site),  
		     0.707106781*(F[F.format.index_r(c,1,site)] 
				  + factor*F[F.format.index_r(c,3,site)]));
      F_out.data.set(F.format.index_i(c, 2, site),  
		     0.707106781*(F[F.format.index_i(c,1,site)]
				  + factor*F[F.format.index_i(c,3,site)]));

      // \psi_3 = -a*(\psi_0 + \psi_2)
      F_out.data.set(F.format.index_r(c, 3, site),  
		     -0.707106781*(F[F.format.index_r(c,0,site)] 
				   + factor*F[F.format.index_r(c,2,site)]));
      F_out.data.set(F.format.index_i(c, 3, site),  
		     -0.707106781*(F[F.format.index_i(c,0,site)] 
				   + factor*F[F.format.index_i(c,2,site)]));
    }  
  }

  return F_out;
}

void BFM_Storage::GaugeExport_to_BFM(GaugeField& U){
  // BFM storage pattern N.B. : EVEN-ODD
  // double gauge [x%2][t][z][y][x/2] [row][column][realimag]  
  using namespace Mapping;
  using namespace FieldUtils;

  double* U_Ptr; 

  GaugeField1D Umu;
  
  for (int mu = 0; mu < NDIM_; mu++) {

    Umu = DirSlice(U, mu);
    CCIO::cout<<"BFM_Storage: Before shiftfield call\n";
    U_dag = shiftField(Umu, mu, Backward()); // U(x-mu)

    int dir = 2*mu+1;
    CCIO::cout << "Direction : "<<dir<<"\n";
    U_Ptr = Umu.data.getaddr(0);
    bfm_obj_.importGauge(U_Ptr, dir); //BFM function

    dir = 2*mu;
    CCIO::cout << "Direction2 : "<<dir<<"\n";
    
    CCIO::cout<<"BFM_Storage: Before u_dag loop call\n";
    for(int site = 0; site<U.format.Nvol(); ++site)
      U_dag.data[U_dag.format.cslice(0,site)] =
	mat_dag(U_dag,site).getva();
    
    CCIO::cout<<"BFM_Storage: Before getting the pointer\n";
    U_Ptr = U_dag.data.getaddr(0);
    CCIO::cout<<"BFM_Storage: Before importGauge call\n";
    bfm_obj_.importGauge(U_Ptr, dir); //BFM function
  }
};

void BFM_Storage::FermionExport_to_BFM(FermionField& F, Fermion_t handle, int cb){
  // BFM storage pattern
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]
  double* F_ptr;

  //Convert Dirac basis: Dirac -> Chiral
  FermionField F_transformed = BasisConversion(F, DIRAC_TO_CHIRAL);

  F_ptr = F_transformed.data.getaddr(0);

  bfm_obj_.importFermion(F_ptr,handle, cb);
  
};


void BFM_Storage::FermionImport_from_BFM(FermionField&F, Fermion_t handle, int cb){
  // BFM storage pattern      
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]    
  FermionField F_temp;
  double* F_ptr = F_temp.data.getaddr(0);
  bfm_obj_.exportFermion(F_ptr,handle, cb);  

  //Convert Dirac basis: Chiral -> Dirac 
  F = BasisConversion(F_temp, CHIRAL_TO_DIRAC);

};

				
