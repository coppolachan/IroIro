/*
  @file StoutSmear.hpp
  @brief Defines Stout smearing class
*/
#include "Smearing/stoutSmear.hpp"
#include "Communicator/comm_io.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include <complex>

#include <omp.h>
#include "timings.hpp"

using namespace std;
typedef complex<double> dcomplex;

//====================================================================
void Smear_Stout::smear(GaugeField& u_smr, const GaugeField& u_in) const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  long double timing;
 

  register int Nvol = CommonPrms::instance()->Nvol();
  GaugeField u_tmp1, q_mu;

  _Message(DEBUG_VERB_LEVEL, "Stout smearing started\n");

  // Smear the configuration
  SmearBase->smear(u_tmp1,u_in);

  double* U_mu_ptr;
  double* utmp_ptr; 
  double* u_smr_ptr;
  double* q_mu_ptr;
  double* inmat;
  double* outmat;
  double trace;
  const int CC2 = 2*NC_*NC_;

  //FINE_TIMING_START(timing);  

 U_mu_ptr  = const_cast<GaugeField&>(u_in).data.getaddr(0);
 utmp_ptr  = u_tmp1.data.getaddr(0);
 u_smr_ptr = u_smr.data.getaddr(0);
 q_mu_ptr  = q_mu.data.getaddr(0);
 
#pragma omp parallel
 { 
   const int nid = omp_get_num_threads();
   const int tid = omp_get_thread_num();
   const int is = tid*Nvol/nid;
   const int ns = Nvol/nid;
   const int str2 = is*CC2*NDIM_;
   
   BGWilsonSU3_MatMult_ND(u_smr_ptr+str2, utmp_ptr+str2, U_mu_ptr+str2, ns*NDIM_);
 
   
   for(int site = is*NDIM_; site < (is+ns)*NDIM_; ++site){
     inmat  = u_smr_ptr + CC2*site;
     outmat = q_mu_ptr  + CC2*site;
       
       trace = (inmat[1]+inmat[9]+inmat[17])/NC_;// imtrace
       for(int a=0; a<NC_; ++a){
	 for(int b=a; b<NC_; ++b){
	   
	   int ab = 2*(NC_*a+b);
	   int ba = 2*(NC_*b+a);

	   *(outmat+ab)   = 0.5*(*(inmat+ab)  -*(inmat+ba));
	   *(outmat+ab+1) = 0.5*(*(inmat+ab+1)+*(inmat+ba+1));
	   
	   *(outmat+ba)   = -*(outmat+ab);
	   *(outmat+ba+1) = *(outmat+ab+1);

	 }
       }
       outmat[1] -= trace;
       outmat[9] -= trace;
       outmat[17] -= trace;
       
     }
   //}

 }

 //CCIO::cout << "q_mu norm: "<< q_mu.norm() << "\n";
 //  FINE_TIMING_END(timing);
 //  _Message(TIMING_VERB_LEVEL, "[Timing] - Smear_Stout::smear - q_mu timing = "
 //	   << timing << std::endl); 


 // FINE_TIMING_START(timing);  
  exponentiate_iQ(u_tmp1, q_mu);// already omp
  //  FINE_TIMING_END(timing);
  // _Message(TIMING_VERB_LEVEL, "[Timing] - Smear_Stout::smear - exp timing = "
  //	   << timing << std::endl); 

  //CCIO::cout << "exp norm: "<< u_tmp1.norm() << "\n";
  //FINE_TIMING_START(timing);  
#pragma omp parallel
  {
    int nid = omp_get_num_threads();
    int tid = omp_get_thread_num();
    
    int is = tid*Nvol / nid;
    int ns = Nvol / nid;
    register int jump = is*CC2*4;
    BGWilsonSU3_MatMult_NN(u_smr_ptr+jump, utmp_ptr+jump, U_mu_ptr+jump, ns*NDIM_);
  } 
  //FINE_TIMING_END(timing);
  // _Message(TIMING_VERB_LEVEL, "[Timing] - Smear_Stout::smear - final step timing = "
  //	   << timing << std::endl); 

  _Message(DEBUG_VERB_LEVEL, "Stout smearing completed \n");
}

void Smear_Stout::exponentiate_iQ(GaugeField& e_iQ,const GaugeField& iQ)const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;
  register int Nvol = e_iQ.format.Nvol();
  register int Nex  = e_iQ.format.Nex();
  const int CC2 = 2*NC_*NC_;
  GaugeField iQ_0;
  GaugeField iQ_1;
  GaugeField iQ_2;

  double* iQ_ptr = const_cast<GaugeField&>(iQ).data.getaddr(0);
  double* iQ0_ptr = iQ_0.data.getaddr(0);
  double* iQ1_ptr = iQ_1.data.getaddr(0);
  double* iQ2_ptr = iQ_2.data.getaddr(0);
  double* e_iQ_ptr = e_iQ.data.getaddr(0);

#pragma omp parallel 
  {
    int tid, nid;
    int ns,is;

    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    
    is = tid*Nvol / nid;
    ns = Nvol / nid;
    const int str2 = is*CC2*Nex;

    double c0, c1, c0max, u_val, w, theta, xi0, u2, w2, cosw, fden;
    dcomplex f0, f1, f2, h0, h1, h2, e2iu, emiu, ixi0, qt;
    
    BGWilsonSU3_MatUnity(iQ0_ptr+str2, ns*Nex);
    BGWilsonSU3_MatMult_NN(iQ1_ptr+str2, iQ_ptr+str2, iQ_ptr+str2, ns*Nex);
    BGWilsonSU3_MatMult_NN(iQ2_ptr+str2, iQ_ptr+str2, iQ1_ptr+str2, ns*Nex);// is iQ3

    for(int site = is*Nex; site < (is+ns)*Nex; ++site){

      register int shift = CC2*site;

	c0 = 0.0;
	c1 = 0.0;
	double* iQ3mat = iQ2_ptr  + shift;
	double* iQ1mat = iQ_ptr   + shift;
	double* iQ2mat = iQ1_ptr  + shift;
	double* iQ0mat = iQ0_ptr  + shift;
	double* eiQmat = e_iQ_ptr + shift;
	for(int cc = 0; cc < NC_; ++cc){
	  c0 += iQ3mat[2*cc*(NC_+1)+1];
	  c1 += iQ2mat[2*cc*(NC_+1) ];
	}
	c0 = -c0/3.0;
	c1 = -c1/2.0;
	c0max = 2.0 * pow(c1/3.0,1.5);
	
	theta = acos(c0/c0max);
	u_val = sqrt(c1/3.0) * cos(theta/3.0);
	w = sqrt(c1) * sin(theta/3.0);
	xi0 = func_xi0(w);
	u2 = u_val * u_val;
	w2 = w * w;
	cosw = cos(w);
	
	ixi0 = dcomplex(0.0,xi0);
	emiu = dcomplex(cos(u_val),-sin(u_val));
	e2iu = dcomplex(cos(2.0*u_val),sin(2.0*u_val));
	
	h0 = e2iu * dcomplex(u2-w2,0.0)
	  + emiu * (dcomplex(8.0*u2*cosw,0.0)
		    + dcomplex(2.0*u_val*(3.0*u2 + w2),0.0)*ixi0);
	h1 = dcomplex(2*u_val,0.0) * e2iu
	  - emiu*(dcomplex(2.0*u_val*cosw,0.0)
		  - dcomplex(3.0*u2-w2,0.0) * ixi0);
	h2 = e2iu - emiu * (dcomplex(cosw,0.0)
			    + dcomplex(3.0*u_val,0.0)*ixi0);
	
	fden = 1.0/(9.0*u2 - w2);
	f0 = h0 * dcomplex(fden,0.0);
	f1 = h1 * dcomplex(fden,0.0);
	f2 = h2 * dcomplex(fden,0.0);
	
	for(int cc = 0; cc < NC_*NC_; ++cc){
	  qt =  f0 * dcomplex(iQ0mat[2*cc], iQ0mat[2*cc+1])
	    + f1 * dcomplex(iQ1mat[2*cc+1],-iQ1mat[2*cc])
	    - f2 * dcomplex(iQ2mat[2*cc], iQ2mat[2*cc+1]);
	  eiQmat[2*cc] = qt.real();
	  eiQmat[2*cc+1] = qt.imag();

	}
      }
    }
 
}

