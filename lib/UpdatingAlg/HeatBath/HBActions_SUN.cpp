//--------------------------------------------------------------------
/*! @file HBActions.cpp
 *
 * @brief Declarations of actions for Heat Bath updates
 *
 *
 *  Time-stamp: <2013-04-10 16:21:33 neo>
 */
//--------------------------------------------------------------------

#include "HBActions_SUN.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include "Tools/randNum.h"
#include "include/numerical_const.hpp"
#include "Main/Geometry/siteIndex_EvenOdd.hpp"

#include <valarray>
#include <complex>

inline double SU2_det_r(SUNmat& u) {
  return  u.r(0,0)*u.r(1,1)- u.i(0,0)*u.i(1,1)-
   u.r(0,1)*u.r(1,0)+ u.i(0,1)*u.i(1,0);

}

inline double SU2_det_i(SUNmat& u) {
  return  u.r(0,0)*u.i(1,1)+ u.i(0,0)*u.r(1,1)-
   u.r(0,1)*u.i(1,0)- u.i(0,1)*u.r(1,0);

}


// Projection for Cabibbo-Marinari
// Any general matrix V 2x2 complex
// can be written 
//   V = v0*1 + i \sigma_i * v_i
// v complex vector
// The real part of v is relevant for
// the heatbath update

inline double ReV0 (SUNmat& u) {
  return 0.5*(u.r(0,0)+u.r(1,1));
}

inline double ReV1 (SUNmat& u) {
  return 0.5*(u.r(0,1)-u.r(1,0));
}

inline double ReV2 (SUNmat& u) {
  return 0.5*(u.i(0,1)+u.i(1,0));
}

inline double ReV3 (SUNmat& u) {
  return 0.5*(u.i(0,0)-u.i(1,1));
}

inline double projected_det(SUNmat& u) {
  // Calculates the determinant of the matrix
  // proportional to SU(2) coming from the 
  // Cabibbo-Marinari algorithm
  return (ReV0(u)*ReV0(u)+ReV1(u)*ReV1(u)+
	  ReV2(u)*ReV2(u)+ReV3(u)*ReV3(u));
}

void HBAction_SUN::init(const RandNum& RNG) {
  //empty
  CCIO::cout << "INIT!!!\n";
};

double HBAction_SUN::calc_H() {
  //Number of plaquettes
  int Nplaq = CommonPrms::instance()->Lvol()*NDIM_*(NDIM_-1)/2.0;
  double plaq = stplfield_.plaquette(*u_);
  double Hgauge = Params_.beta*Nplaq*(1.0 -plaq);

  _Message(ACTION_VERB_LEVEL,"    [HeatBathSUNGaugeWilson] H = "<<Hgauge <<"\n");
  _Message(ACTION_VERB_LEVEL,"    [HeatBathSUNGaugeWilson] Plaquette = "<<plaq<<"\n");
  
  return Hgauge;

};

void HBAction_SUN::hb_update(const RandNum& RNG) {
  using namespace SUNmatUtils; 
  using namespace FieldUtils;

  if (NC_!=2) {
    CCIO::cout << "Error [HBActions_SUN.cpp]: hb_update defined only for SU(2)\n";
    abort();
  }

  GaugeField1D tmp, tmpup, tmpdn; 
  SUNmat staple, newlink;
  int sector = 0;

  double beta_adj_p, sigma;
  // Case of mixed actions (fundamental/adjoint)
  // Generate the auxiliary fields
  //Field zeta(Nvol_*2);
  Field fp(2*Nvol_*(NDIM_-2)*(NDIM_-1));
  Field zp(2*Nvol_*(NDIM_-2)*(NDIM_-1));
  double z_r, z_i;
  
  beta_adj_p = 0.0;
  if (Params_.beta_adj !=0.0) {
    beta_adj_p = Params_.beta_adj*NC_*NC_/(NC_*NC_-1);
    sigma = sqrt(1.0/(1.0*beta_adj_p));
  } else sigma = 1.0;//protects from overflow
  
  // Fill the field with complex numbers
  // distributed normally with variance = sigma
  /*
  for (int m = 0; m < NDIM_-1; m++) {
    for (int n = m+1; n < NDIM_; n++) {
      tmp = stplfield_.upper(*u_, m,n);
      for (int s = 0; s < Nvol_; s++) {
	RNG.get_complex_gauss(&z_r, &z_i, sigma);
	double plaqr = ReTr(mat(*u_,s,m)*mat_dag(tmp,s));
	//double plaqi = ImTr(mat(*u_,s,m)*mat_dag(tmp,s));
	zp.set(2*s+  2*Nvol_*sector  , z_r);
	zp.set(2*s+1+2*Nvol_*sector  , z_i);

	fp.set(2*s+  2*Nvol_*sector  , Params_.beta/NC_+2.0*beta_adj_p/NC_*(z_r+ plaqr/(NC_)));
	//fp.set(2*s+1+2*Nvol_*sector  ,             +2.0*beta_adj_p*(z_i+plaqi/NC_));
	fp.set(2*s+1+2*Nvol_*sector  ,             0.0);

	//CCIO::cout << "fp[s] = "<<fp[2*s+  2*Nvol_*sector] << "  fp[s+1] = "<< fp[2*s+1+  2*Nvol_*sector] <<"\n";
	//CCIO::cout << "z_r  "<<z_r  <<"\n";
	//CCIO::cout << "s= "<<s<< " plr = "<<plaqr << "  pli = "<< plaqi <<"\n";
      }
      sector++;
    }
  }
  */

  // SU(2) updating
 
  
  std::vector<int> Sites; // stores even/odd information
  std::valarray<double> rand_num;
  std::valarray<double> rand_vec;
  double lambda_sq, sum;
  double q0, q1, q2, q3;
  double u0, u1, u2, u3;
  double x0;
  int site;
  
  rand_num.resize(4);
  rand_vec.resize(3);
  
 

  CCIO::cout << "Heat-bath updating...\n";
  for (int eo = 0; eo < 2; eo++) { 
    
    // use even-odd checkerboard 
    if (eo) { 
      Sites = SiteIndex_EvenOdd::instance()->esec();
    } else {
      Sites = SiteIndex_EvenOdd::instance()->osec();
    }
    for(int m = 0; m < NDIM_; ++m){
      
      //Create Auxiliary field fp
      sector = 0;
      for (int mu = 0; mu < NDIM_-1; mu++) {
	for (int n = mu+1; n < NDIM_; n++) {
	  tmp = stplfield_.upper(*u_, mu,n);
	  for (int s = 0; s < Nvol_; s++) {
	    site = s;
	    RNG.get_complex_gauss(&z_r, &z_i, sigma);
	    double plaqr = ReTr(mat(*u_,site,mu)*mat_dag(tmp,site));
	    
	    fp.set(2*site+  2*Nvol_*sector  , Params_.beta/NC_+2.0*beta_adj_p/NC_*(z_r+plaqr/(NC_)));
	  
	    fp.set(2*site+1+2*Nvol_*sector  , 0.0);
	  }
	  sector++;
	  
	}
      }
      
      tmp = 0.0;
       // Fill the staple field, tmp
      for(int n=0; n< NDIM_; ++n) {
	if(n != m) 
	  {
	    tmp += stplfield_.upper_lower(*u_,m,n, fp);
	  }
      }
	  
      
      // Sweep among links 
      // first even (eo=0), then odd (eo=1)
      for (int s = 0; s < Nvol_/2; s++) {
	site = Sites[s];
	staple = mat(tmp, site);

	// Calculate xi = det (staple)^(1/2)
	double xi    = SU2_det_r(staple); // for SU(2) det is always real and >0
	//double xi = projected_det(staple);// more general, little overhead
	/*
	if (xi < 0 )
	  {
	    CCIO::cout << "Warning: bad xi = "<< xi <<"\n";
	    abort();
	  };
	*/	  
	
	// Generating the new link
	// Step 1
	// Generate 4 random numbers distributed uniformly in (0,1]
	// Kennedy Pendleton Phys. Lett. 156B (1985) 393
	xi = sqrt(xi);
	
	do {
	  RNG.get(rand_num); // generates in [0,1)
	  rand_num[0] = 1.0 - rand_num[0];
	  rand_num[1] = 1.0 - rand_num[1];
	  rand_num[2] = 1.0 - rand_num[2];
	  /*
	  lambda_sq = -1.0/(2.0*xi*Params_.beta)*
	    (log(rand_num[0])+
	     cos(2.0*PI*rand_num[1])*cos(2.0*PI*rand_num[1])*log(rand_num[2]));
	  */
	  lambda_sq = -1.0/(2.0*NC_*xi)*
	    (log(rand_num[0])+
	     cos(2.0*PI*rand_num[1])*cos(2.0*PI*rand_num[1])*log(rand_num[2]));
	}	while ((rand_num[3]*rand_num[3]) > (1.0 - lambda_sq));
	//CCIO::cout << "lambda = "<<lambda_sq  <<"\n";
	x0 = 1.0 - 2.0 * lambda_sq;
	
	// Step 2
	do {
	  RNG.get(rand_vec); // generates in [0,1)
	  rand_vec *= 2.0;
	  rand_vec -= 1.0; //now in [-1,1)
	  sum = rand_vec[0]*rand_vec[0]+rand_vec[1]*rand_vec[1]+rand_vec[2]*rand_vec[2];
	}	while (sum > 1);
	
	// Normalize to sqrt(1-x0^2)
	double final_norm = sqrt(1.0 - x0*x0);
	double norm = sqrt(sum);
	rand_vec *= final_norm/norm;
	
	/* Create SU(2) matrix:
	        x0 + i x3      x1 + i x2
	   u = 
	       -x1 + i x2      x0 - i x3
	   
	   x1, x2, x3 elements of rand_vec
	*/
	
	
	  // Original pure SU(2) update
	
	newlink.set(0,0, x0          ,  rand_vec[2]);
	newlink.set(0,1,  rand_vec[0],  rand_vec[1]);
	newlink.set(1,0, -rand_vec[0],  rand_vec[1]);
	newlink.set(1,1, x0          , -rand_vec[2]);
	newlink *= staple;
	newlink /= xi;
	

	/*
	q0 =  ReV0(staple)/xi;
	q1 =  ReV1(staple)/xi;
	q2 =  ReV2(staple)/xi;
	q3 =  ReV3(staple)/xi;
	
	u0 = x0*q0-rand_vec[2]*q3-rand_vec[0]*q1-rand_vec[1]*q2;
	u1 = x0*q1-rand_vec[2]*q2+rand_vec[0]*q0+rand_vec[1]*q3;
	u2 = x0*q2+rand_vec[2]*q1-rand_vec[0]*q3+rand_vec[1]*q0;
	u3 = x0*q3+rand_vec[2]*q0+rand_vec[0]*q2-rand_vec[1]*q1;
	
	newlink.set(0,0,  u0,  u3 );
        newlink.set(0,1,  u1,  u2 );
	newlink.set(1,0, -u1,  u2 );
	newlink.set(1,1,  u0, -u3 );
	*/
	
	SetMat(*u_, newlink, site, m);


      } // end loop on sites for updating U

      

    } // end loop on mu
  } // end loop on even/odd checkerboard (essential)
  

};


