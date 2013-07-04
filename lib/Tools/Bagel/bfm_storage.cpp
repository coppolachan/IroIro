
//---------------------------------------------------------------------
/*! @file bfm_storage.cpp
  @brief Routines to convert IroIro fields to BFM (Bagel) type

  Class definitions
*/ 
//---------------------------------------------------------------------

#include "macros.hpp" //defines Color, Spin, Dimension constants
#include "commonPrms.hpp"

#include "include/numerical_const.hpp"
#include "bfm_storage.hpp"
#include "lib/Geometry/shiftField.hpp"
#include "Tools/fieldUtils.hpp"



BFM_Storage::BFM_Storage(bfm_dp& bfm):Nx_(CommonPrms::instance()->Nx()),
				      Ny_(CommonPrms::instance()->Ny()),
				      Nz_(CommonPrms::instance()->Nz()),
				      Nt_(CommonPrms::instance()->Nt()),
				      Nvol_(CommonPrms::instance()->Nvol()),
				      bfm_obj_(bfm){};


void BFM_Storage::BasisConversion(FermionField& F_out, FermionField& F, int ConvType, int cb, int s = 0){
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

  for(int site = Nvol_*(2*s+cb)/2; site<Nvol_/2*(1+cb+2*s); ++site){
    for(int c = 0; c<NC_; ++c){
      //write explicitly the spinor indexes multiplication
 
     // \psi_0 = a*(\psi_1 - \psi_3)
      F_out.data.set(F.format.index_r(c, 0, site) , 
		     INVSQRT2*(factor*F[F.format.index_r(c,1,site)]
				  - F[F.format.index_r(c,3,site)]));
      F_out.data.set(F.format.index_i(c, 0, site) , 
		     INVSQRT2*(factor*F[F.format.index_i(c,1,site)] 
				  - F[F.format.index_i(c,3,site)]));

      // \psi_1 = a*(\psi_2 - \psi_0)
      F_out.data.set(F.format.index_r(c, 1, site) ,
		     INVSQRT2*(F[F.format.index_r(c,2,site)] 
				  - factor*F[F.format.index_r(c,0,site)]));
      F_out.data.set(F.format.index_i(c, 1, site),  
		     INVSQRT2*(F[F.format.index_i(c,2,site)]
				  - factor*F[F.format.index_i(c,0,site)]));

      // \psi_2 = a*(\psi_1 + \psi_3)
      F_out.data.set(F.format.index_r(c, 2, site),  
		     INVSQRT2*(F[F.format.index_r(c,1,site)] 
				  + factor*F[F.format.index_r(c,3,site)]));
      F_out.data.set(F.format.index_i(c, 2, site),  
		     INVSQRT2*(F[F.format.index_i(c,1,site)]
				  + factor*F[F.format.index_i(c,3,site)]));

      // \psi_3 = -a*(\psi_0 + \psi_2)
      F_out.data.set(F.format.index_r(c, 3, site),  
		     -INVSQRT2*(F[F.format.index_r(c,0,site)] 
				   + factor*F[F.format.index_r(c,2,site)]));
      F_out.data.set(F.format.index_i(c, 3, site),  
		     -INVSQRT2*(F[F.format.index_i(c,0,site)] 
				   + factor*F[F.format.index_i(c,2,site)]));
    }  
  }

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
    U_dag = shiftField(Umu, mu, Backward()); // U(x-mu)

    int dir = 2*mu+1;
    //    CCIO::cout << "Direction : "<<dir<<"\n";
    U_Ptr = Umu.data.getaddr(0);
    bfm_obj_.importGauge(U_Ptr, dir); //BFM function

    dir = 2*mu;
    //    CCIO::cout << "Direction2 : "<<dir<<"\n";
    
    for(int site = 0; site<U.format.Nvol(); ++site)
      U_dag.data[U_dag.format.cslice(0,site)] =
	mat_dag(U_dag,site).getva();

    U_Ptr = U_dag.data.getaddr(0);
    bfm_obj_.importGauge(U_Ptr, dir); //BFM function
  }
};

void BFM_Storage::GaugeExport_to_BFM(const Field* U){
  GaugeField GF = GaugeField(*U);
  GaugeExport_to_BFM(GF);
}

void BFM_Storage::FermionExport_to_BFM(FermionField& F, Fermion_t handle, int cb){
  // BFM storage pattern
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]
  double* F_ptr;

  //Convert Dirac basis: Dirac -> Chiral
  FermionField F_transformed(F.Nvol());
  BasisConversion(F_transformed, F, DIRAC_TO_CHIRAL,cb);
  
  F_ptr = F_transformed.data.getaddr(0);
  bfm_obj_.importFermion(F_ptr,handle, cb);
  
};

void BFM_Storage::FermionExport_to_BFM_5D(FermionField& F, Fermion_t handle, int cb, int s){
  // BFM storage pattern
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]
  double* F_ptr;
  int Ls = F.Nvol()/Nvol_;

  //Convert Dirac basis: Dirac -> Chiral
  FermionField F_transformed(F.Nvol());
  BasisConversion(F_transformed, F, DIRAC_TO_CHIRAL,cb,s);

  F_ptr = F_transformed.data.getaddr(0);
  
  bfm_obj_.impexFermion(F_ptr+Nvol_*F.Nin()*s,handle,1, cb,Ls-1-s);
  
};


void BFM_Storage::FermionImport_from_BFM(FermionField&F, Fermion_t handle, int cb){
  // BFM storage pattern      
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]    
  FermionField F_temp(F.Nvol());//to accomodate for 5d fermions
  double* F_ptr = F_temp.data.getaddr(0);
  bfm_obj_.exportFermion(F_ptr,handle, cb);  

  //Convert Dirac basis: Chiral -> Dirac 
  BasisConversion(F, F_temp, CHIRAL_TO_DIRAC,cb);

};

void BFM_Storage::FermionImport_from_BFM_5D(FermionField&F, Fermion_t handle, int cb, int s){
  // BFM storage pattern      
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]  
  FermionField F_temp(F.Nvol());//to accomodate for 5d fermions
  double* F_ptr = F_temp.data.getaddr(0);
  int Ls = F.Nvol()/Nvol_;
  
  bfm_obj_.impexFermion(F_ptr+Nvol_*F.Nin()*s,handle, 0, cb, Ls - 1-s);  

  //Convert Dirac basis: Chiral -> Dirac 
  BasisConversion(F, F_temp, CHIRAL_TO_DIRAC,cb,s);

};

				
