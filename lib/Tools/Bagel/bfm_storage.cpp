
//---------------------------------------------------------------------
/*! @file bfm_storage.cpp
  @brief Routines to convert IroIro fields to BFM (Bagel) type

  Class definitions
  
  Written by Guido Cossu
*/ 
//---------------------------------------------------------------------

#include "macros.hpp" //defines Color, Spin, Dimension constants
#include "commonPrms.hpp"

#include "include/numerical_const.hpp"
#include "bfm_storage.hpp"
#include "lib/Geometry/shiftField.hpp"
#include "Tools/fieldUtils.hpp"

#include "omp.h"

template <class Float>
BFM_Storage<Float>::BFM_Storage(bfm_internal<Float>& bfm):Nx_(CommonPrms::instance()->Nx()),
							  Ny_(CommonPrms::instance()->Ny()),
							  Nz_(CommonPrms::instance()->Nz()),
							  Nt_(CommonPrms::instance()->Nt()),
							  Nvol_(CommonPrms::instance()->Nvol()),
							  bfm_obj_(bfm){};

template <class Float>
void BFM_Storage<Float>::BasisConversion(FermionField& F_out, FermionField& F, int ConvType, int cb, int s){
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

template <class Float>
void BFM_Storage<Float>::GaugeExport_to_BFM(GaugeField& U){
  // BFM storage pattern N.B. : EVEN-ODD
  // double gauge [x%2][t][z][y][x/2] [row][column][realimag]  
  using namespace Mapping;
  using namespace FieldUtils;
  double* U_Ptr; 

  GaugeField1D Umu;

  CCIO::cout << "U (IroIro) norm :" << U.data.norm() << "\n";

  for (int mu = 0; mu < NDIM_; mu++) {
    //CCIO::cout << "Dirslice\n";
    Umu = DirSlice(U, mu);
    //CCIO::cout << "ShiftField \n";
    U_dag = shiftField(Umu, mu, Backward()); // U(x-mu)

    int dir = 2*mu+1;
    //CCIO::cout << "Direction : "<<dir<<"\n";
    U_Ptr = Umu.data.getaddr(0);
    //CCIO::cout << "ImportGauge BFM function...";
    bfm_obj_.importGauge(U_Ptr, dir); //BFM function
    //CCIO::cout << "complete\n";
    dir = 2*mu;
    //CCIO::cout << "Direction2 : "<<dir<<"\n";
    
    for(int site = 0; site<U.format.Nvol(); ++site)
      U_dag.data[U_dag.format.cslice(0,site)] =
	mat_dag(U_dag,site).getva();

    U_Ptr = U_dag.data.getaddr(0);

    //CCIO::cout << "ImportGauge BFM function 2...";
    bfm_obj_.importGauge(U_Ptr, dir); //BFM function
    //CCIO::cout << "complete\n";
    

  }

  double matnorm;

#pragma omp parallel for 
    for (int t=0; t < bfm_obj_.nthread ; t++)
      matnorm = bfm_obj_.normMat((Matrix_t)bfm_obj_.u);


  CCIO::cout << "DEBUG U (BFM) norm : " << matnorm << "\n"; 

};

template <class Float>
void BFM_Storage<Float>::GaugeExport_to_BFM(const Field* U){
  GaugeField GF = GaugeField(*U);
  GaugeExport_to_BFM(GF);
}

template <class Float>
void BFM_Storage<Float>::GaugeImport_from_BFM(Field* U, Matrix_t handle, int dir, int cb){
  int Ndircoco=36; /*4 directions stored*/
  int Ncoco=9;
  Format::Format_G GaugeFormat(Nvol_);
  double *U_ptr = U->getaddr(0);
  omp_set_num_threads(bfm_obj_.nthread);
#pragma omp parallel 
  {    
#pragma omp for 
    for (int site=0;site<Nx_*Ny_*Nz_*Nt_;site++ ) { 
      
      int x[4] ;
      int s=site;
      x[0]=s%Nx_;    s=s/Nx_;
      x[1]=s%Ny_;    s=s/Ny_;
      x[2]=s%Nz_;    s=s/Nz_;
      x[3]=s%Nt_;
      
      if ( ((x[0]+x[1]+x[2]+x[3])&0x1) == cb ) {      
	int bbase = dir*9;
	int siteIdxEO = SiteIndex::instance()->site(x[0], x[1], x[2], x[3]);
	for ( int coco=0;coco<9;coco++ ) { 
	  for ( int reim=0;reim<2;reim++ ) { 
	    double* bagel = (double *)handle;
	    int idx = GaugeFormat.index(2*coco+reim, siteIdxEO,dir);
	    int bidx = bfm_obj_.bagel_idx(x,reim,coco+bbase,Ndircoco,1);
	    // CCIO::cout << "idx = "<<idx<< " bidx = "<<bidx << "\n";
	    //CCIO::cout << "bagel: "<<bagel[bidx] <<"\n";
	    U_ptr[idx] = bagel[bidx];
	    
	  }
	}
      }
    }
  }
}


template <class Float>
void BFM_Storage<Float>::FermionExport_to_BFM(FermionField& F, Fermion_t handle, int cb){
  // BFM storage pattern
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]
  double* F_ptr;

  //Convert Dirac basis: Dirac -> Chiral
  FermionField F_transformed(F.Nvol());
  BasisConversion(F_transformed, F, DIRAC_TO_CHIRAL,cb);
  
  F_ptr = F_transformed.data.getaddr(0);
  bfm_obj_.importFermion(F_ptr,handle, cb);
  
};

template <class Float>
void BFM_Storage<Float>::FermionExport_to_BFM_5D(FermionField& F, Fermion_t handle, int cb, int s){
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

template <class Float>
void BFM_Storage<Float>::FermionImport_from_BFM(FermionField&F, Fermion_t handle, int cb){
  // BFM storage pattern      
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]    
  FermionField F_temp(F.Nvol());//to accomodate for 5d fermions
  double* F_ptr = F_temp.data.getaddr(0);
  bfm_obj_.exportFermion(F_ptr,handle, cb);  

  //Convert Dirac basis: Chiral -> Dirac 
  BasisConversion(F, F_temp, CHIRAL_TO_DIRAC,cb);

};

template <class Float>
void BFM_Storage<Float>::FermionImport_from_BFM_5D(FermionField&F, Fermion_t handle, int cb, int s, int Ls){
  // BFM storage pattern      
  // double psi   [x%2][t][z][y][x/2] [spin][color][realimag]  
  FermionField F_temp(Nvol_*Ls);//to accomodate for 5d fermions
  double* F_ptr = F_temp.data.getaddr(0);

  bfm_obj_.impexFermion(F_ptr+Nvol_*F.Nin()*s,handle, 0, cb, Ls - 1-s);  

  //Convert Dirac basis: Chiral -> Dirac 
  BasisConversion(F, F_temp, CHIRAL_TO_DIRAC,cb,s);

  

};

				
template class BFM_Storage<double>;
template class BFM_Storage<float>;
