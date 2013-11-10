/*!
  @file smartConf.cpp
  @brief Defines the SmartConf class member functions
*/
#include "Smearing/smartConf.hpp"
#include "Geometry/autoMap.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"



#include <omp.h>

typedef std::complex<double> dcomplex;



void Trace(double* re, double* im, double* mat1, double* mat2){
  *re = 0;
  *im = 0;
  for (int c0 = 0; c0< NC_; ++c0){
    *re += mat1[6*c0  ]*mat2[2*c0   ] - mat1[6*c0+1]*mat2[2*c0+1 ];
    *re += mat1[6*c0+2]*mat2[2*c0+6 ] - mat1[6*c0+3]*mat2[2*c0+7 ];
    *re += mat1[6*c0+4]*mat2[2*c0+12] - mat1[6*c0+5]*mat2[2*c0+13];

    *im += mat1[6*c0  ]*mat2[2*c0+1 ] + mat1[6*c0+1]*mat2[2*c0   ];
    *im += mat1[6*c0+2]*mat2[2*c0+7 ] + mat1[6*c0+3]*mat2[2*c0+6 ];
    *im += mat1[6*c0+4]*mat2[2*c0+13] + mat1[6*c0+5]*mat2[2*c0+12];   

  }
}

//====================================================================
void SmartConf::smeared_force(GaugeField& SigmaTilde)const {
  GaugeField force = SigmaTilde;//actually = U*SigmaTilde

  int Nvol = CommonPrms::instance()->Nvol();
  const int CC2 = 2*NC_*NC_;

  double* ThinLinks_ptr;
  double* force_ptr;
  double* SigmaTilde_ptr;
  double* SmearedSet_ptr;
  
  SmearedSet_ptr = const_cast<Field&>(SmearedSet[smearingLevels-1].data).getaddr(0); 
  force_ptr      = force.data.getaddr(0);  
  ThinLinks_ptr  = ThinLinks->data.getaddr(0);  
  SigmaTilde_ptr = SigmaTilde.data.getaddr(0);

#pragma omp parallel
 { 
   const int nid = omp_get_num_threads();
   const int tid = omp_get_thread_num();
   const int is = tid*Nvol/nid;
   const int ns = Nvol/nid;
   const int str2 = is*CC2*NDIM_;
   BGWilsonSU3_MatMult_DN(force_ptr+str2, SmearedSet_ptr+str2, force_ptr+str2, ns*NDIM_);
 }

  for(int ismr = smearingLevels - 1; ismr > 0; --ismr)
    force = AnalyticSmearedForce(force,get_smeared_conf(ismr-1));
  
  force = AnalyticSmearedForce(force,*ThinLinks);

#pragma omp parallel
 { 
   const int nid = omp_get_num_threads();
   const int tid = omp_get_thread_num();
   const int is = tid*Nvol/nid;
   const int ns = Nvol/nid;
   const int str2 = is*CC2*NDIM_;
   BGWilsonSU3_MatMult_NN(SigmaTilde_ptr+str2, ThinLinks_ptr+str2, force_ptr+str2, ns*NDIM_);
 }
}
//====================================================================
GaugeField SmartConf::AnalyticSmearedForce(const GaugeField& SigmaKPrime,
					   const GaugeField& GaugeK) const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  int Nvol = CommonPrms::instance()->Nvol();
  GaugeField iQ, e_iQ;
  GaugeField C, iLambda;
  GaugeField SigmaK;
  


  StoutSmearing.BaseSmear(C,GaugeK);

  double* e_iQ_ptr = e_iQ.data.getaddr(0);
  double* iQ_ptr   = iQ.data.getaddr(0);
  double* SigmaK_ptr          = SigmaK.data.getaddr(0);
  double* SigmaKPrime_ptr     = const_cast<Field&>(SigmaKPrime.data).getaddr(0);
  double* C_ptr               = C.data.getaddr(0);
  double* iLambda_ptr         = iLambda.data.getaddr(0);

  double* inmat;
  double* outmat;
  double trace;
  const int CC2 = 2*NC_*NC_;

#pragma omp parallel
 { 
   const int nid = omp_get_num_threads();
   const int tid = omp_get_thread_num();
   const int is = tid*Nvol/nid;
   const int ns = Nvol/nid;
   const int str2 = is*CC2*NDIM_;
   BGWilsonSU3_MatMult_ND(e_iQ_ptr+ str2, C_ptr + str2, 
			   const_cast<Field&>(GaugeK.data).getaddr(0)+ str2, 
			   ns*NDIM_);

   //Traceless antihermitian part
   for(int site = is*NDIM_; site < (is+ns)*NDIM_; ++site){
     inmat  = e_iQ_ptr + CC2*site;
     outmat = iQ_ptr  + CC2*site;
       
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

 }
 
 set_iLambda(iLambda, e_iQ, iQ, SigmaKPrime, GaugeK);
 
#pragma omp parallel
 {
   const int nid = omp_get_num_threads();
   const int tid = omp_get_thread_num();
   const int is = tid*Nvol/nid;
   const int ns = Nvol/nid;
   const int str2 = is*CC2*NDIM_;
   
   BGWilsonSU3_MatMult_NN(SigmaK_ptr+str2, SigmaKPrime_ptr+str2, e_iQ_ptr+str2, ns*NDIM_);
   BGWilsonSU3_MatMultAdd_DN(SigmaK_ptr+str2, SigmaK_ptr+str2, C_ptr+str2, iLambda_ptr+str2, ns*NDIM_);
  }

  StoutSmearing.derivative(SigmaK, iLambda, GaugeK);
  
  return SigmaK;
}
//====================================================================
void SmartConf::set_iLambda(GaugeField& iLambda, 
			    GaugeField& e_iQ,
			    const GaugeField& iQ,
			    const GaugeField& Sigmap, 
			    const GaugeField& GaugeK)const{

 
  register int Nvol = CommonPrms::instance()->Nvol();
  register int Nex  = e_iQ.format.Nex();
  const int CC2 = 2*NC_*NC_; 

  GaugeField iQ_0;
  GaugeField iQ_1;
  GaugeField iQ_2;
  GaugeField USigmap;
  GaugeField iQUS;
  GaugeField iUSQ;
  double* iQ_ptr = const_cast<GaugeField&>(iQ).data.getaddr(0);
  double* Sigmap_ptr = const_cast<GaugeField&>(Sigmap).data.getaddr(0);
  double* GaugeK_ptr = const_cast<GaugeField&>(GaugeK).data.getaddr(0);
  double* iQ0_ptr = iQ_0.data.getaddr(0);
  double* iQ1_ptr = iQ_1.data.getaddr(0);
  double* iQ2_ptr = iQ_2.data.getaddr(0);
  double* e_iQ_ptr = e_iQ.data.getaddr(0);
  double* USigmap_ptr = USigmap.data.getaddr(0);
  double* iQUS_ptr = iQUS.data.getaddr(0);
  double* iUSQ_ptr = iUSQ.data.getaddr(0);
  double* iLambda_ptr = iLambda.data.getaddr(0);

#pragma omp parallel 
  {
    int tid, nid;
    int ns,is;
    
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    
    is = tid*Nvol / nid;
    ns = Nvol / nid;
    const int str2 = is*CC2*Nex;
    
    double c0, c1, c0max, theta;
    double u, w, u2, w2, cosw, xi0, xi1, fden;
    dcomplex f0, f1, f2, h0, h1, h2, e2iu, emiu, qt;
    dcomplex r01, r11, r21, r02, r12, r22, tr1, tr2;
    dcomplex b10, b11, b12, b20, b21, b22;
    
    double B1[2*NC_*NC_];
    double B2[2*NC_*NC_];
    double iGamma[2*NC_*NC_];
    
    
    BGWilsonSU3_MatUnity(iQ0_ptr+str2, ns*Nex);
    BGWilsonSU3_MatMult_NN(iQ1_ptr+str2, iQ_ptr+str2, iQ_ptr+str2, ns*Nex);
    BGWilsonSU3_MatMult_NN(iQ2_ptr+str2, iQ_ptr+str2, iQ1_ptr+str2, ns*Nex);// is iQ3
    BGWilsonSU3_MatMult_NN(USigmap_ptr+str2, GaugeK_ptr+str2, Sigmap_ptr+str2, ns*Nex);
    BGWilsonSU3_MatMult_NN(iQUS_ptr+str2, iQ_ptr+str2, USigmap_ptr+str2, ns*Nex);
    BGWilsonSU3_MatMult_NN(iUSQ_ptr+str2, USigmap_ptr+str2, iQ_ptr+str2, ns*Nex);
    
    
    for(int  site = is*Nex; site < (is+ns)*Nex; ++site){
      register int shift = CC2*site;
      c0 = 0.0;
      c1 = 0.0;
      double* iQ3mat = iQ2_ptr  + shift;
      double* iQ1mat = iQ_ptr   + shift;
      double* iQ2mat = iQ1_ptr  + shift;
      double* iQ0mat = iQ0_ptr  + shift;
      double* eiQmat = e_iQ_ptr + shift;
      double* USigmap_mat = USigmap_ptr + shift;
      double* iQUSmat = iQUS_ptr + shift;
      double* iUSQmat = iUSQ_ptr + shift;
      for(int cc = 0; cc < NC_; ++cc){
	c0 += iQ3mat[2*cc*(NC_+1)+1];
	c1 += iQ2mat[2*cc*(NC_+1) ];
      }
      c0 = -c0/3.0;
      c1 = -c1/2.0;
      c0max = 2.0 * pow(c1/3.0,1.5);
      
      theta = acos(c0/c0max);
      u = sqrt(c1/3.0) * cos(theta/3.0);
      w = sqrt(c1) * sin(theta/3.0);
      
      set_fj(f0,f1,f2,u,w);

      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  f0 * dcomplex(iQ0mat[2*cc], iQ0mat[2*cc+1])
	  + f1 * dcomplex(iQ1mat[2*cc+1],-iQ1mat[2*cc])
	  - f2 * dcomplex(iQ2mat[2*cc], iQ2mat[2*cc+1]);
	eiQmat[2*cc]   = qt.real();
	eiQmat[2*cc+1] = qt.imag();
      }
      
      xi0 = func_xi0(w);
      xi1 = func_xi1(w);
      u2 = u * u;
      w2 = w * w;
      cosw = cos(w);
      
      emiu = dcomplex(cos(u),-sin(u));
      e2iu = dcomplex(cos(2.0*u),sin(2.0*u));
      
      r01 = dcomplex(2.0*u,2.0*(u2-w2)) * e2iu
	+ emiu * dcomplex(16.0*u*cosw + 2.0*u*(3.0*u2+w2)*xi0,
			  -8.0*u2*cosw + 2.0*(9.0*u2+w2)*xi0);
      
      r11 = dcomplex(2.0,4.0*u) * e2iu
	+ emiu * dcomplex(-2.0*cosw + (3.0*u2-w2)*xi0,
			  2.0*u*cosw + 6.0*u*xi0);
      
      r21 = dcomplex(0.0,2.0) * e2iu
	+ emiu * dcomplex(-3.0*u*xi0, cosw - 3.0*xi0);
      
      r02 = dcomplex(-2.0,0.0) * e2iu
	+ emiu * dcomplex(-8.0*u2*xi0,
			  2.0*u*(cosw + xi0 + 3.0*u2*xi1));
      
      r12 = emiu * dcomplex(2.0*u*xi0,
			    -cosw - xi0 + 3.0*u2*xi1);
      
      r22 = emiu * dcomplex(xi0, -3.0*u*xi1);
      
      fden = 1.0/(2*(9.0*u2-w2)*(9.0*u2-w2));
      
      b10 =  dcomplex(2.0*u, 0.0)*r01 + dcomplex(3.0*u2-w2, 0.0)*r02
	- dcomplex(30.0*u2+2.0*w2, 0.0)*f0;
      
      b11 =  dcomplex(2.0*u, 0.0)*r11 + dcomplex(3.0*u2-w2, 0.0)*r12
	- dcomplex(30.0*u2+2.0*w2, 0.0)*f1;
      
      b12 =  dcomplex(2.0*u, 0.0)*r21 + dcomplex(3.0*u2-w2, 0.0)*r22
	- dcomplex(30.0*u2+2.0*w2, 0.0)*f2;
      
      b20 = r01 - dcomplex(3.0*u, 0.0)*r02 - dcomplex(24.0*u, 0.0)*f0;
      
      b21 = r11 - dcomplex(3.0*u, 0.0)*r12 - dcomplex(24.0*u, 0.0)*f1;
      
      b22 = r21 - dcomplex(3.0*u, 0.0)*r22 - dcomplex(24.0*u, 0.0)*f2;
      
      b10 *= dcomplex(fden,0.0);
      b11 *= dcomplex(fden,0.0);
      b12 *= dcomplex(fden,0.0);
      b20 *= dcomplex(fden,0.0);
      b21 *= dcomplex(fden,0.0);
      b22 *= dcomplex(fden,0.0);

      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  b10 * dcomplex(iQ0mat[2*cc], iQ0mat[2*cc+1])
	  + b11 * dcomplex(iQ1mat[2*cc+1],-iQ1mat[2*cc])
	  - b12 * dcomplex(iQ2mat[2*cc], iQ2mat[2*cc+1]);
	B1[2*cc]   = qt.real();
	B1[2*cc+1] = qt.imag();
	
	qt =  b20 * dcomplex(iQ0mat[2*cc], iQ0mat[2*cc+1])
	  + b21 * dcomplex(iQ1mat[2*cc+1],-iQ1mat[2*cc])
	  - b22 * dcomplex(iQ2mat[2*cc], iQ2mat[2*cc+1]);
	B2[2*cc]   = qt.real();
	B2[2*cc+1] = qt.imag();
      }
      
   
      double tr_r, tr_i;
      Trace(&tr_r, &tr_i, USigmap_mat, B1);
      tr1 = dcomplex(tr_r, tr_i);
      Trace(&tr_r, &tr_i, USigmap_mat, B2);
      tr2 = dcomplex(tr_r, tr_i);
      
      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  tr1 * dcomplex(iQ1mat[2*cc+1],-iQ1mat[2*cc])
	  - tr2 * dcomplex(iQ2mat[2*cc], iQ2mat[2*cc+1])
	  + f1  * dcomplex(USigmap_mat[2*cc], USigmap_mat[2*cc+1])
	  + f2  * dcomplex(iQUSmat[2*cc+1],-iQUSmat[2*cc])
	  + f2  * dcomplex(iUSQmat[2*cc+1],-iUSQmat[2*cc]);
	iGamma[2*cc]   = -qt.imag();
	iGamma[2*cc+1] =  qt.real();
      }

      double *outmat = iLambda_ptr  + CC2*site;
      double trace = (iGamma[1]+iGamma[9]+iGamma[17])/NC_;// imtrace
      for(int a=0; a<NC_; ++a){
	for(int b=a; b<NC_; ++b){
	  
	  int ab = 2*(NC_*a+b);
	  int ba = 2*(NC_*b+a);
	  
	  *(outmat+ab)   = 0.5*(*(iGamma+ab)  -*(iGamma+ba));
	  *(outmat+ab+1) = 0.5*(*(iGamma+ab+1)+*(iGamma+ba+1));
	  
	  *(outmat+ba)   = -*(outmat+ab);
	  *(outmat+ba+1) = *(outmat+ab+1);
	  
	}
      }
      outmat[1] -= trace;
      outmat[9] -= trace;
      outmat[17] -= trace;
      
      
    }
  }
  
}

//====================================================================
void SmartConf::set_uw(double& u, double& w,
		       const SUNmat& iQ2,
		       const SUNmat& iQ3)const{
  double c0 = 0.0;
  double c1 = 0.0;

  for(int cc = 0; cc < NC_; ++cc){
    c0 += iQ3.i(cc,cc);
    c1 += iQ2.r(cc,cc);
  }
  c0 = -c0/3.0;
  c1 = -c1/2.0;
  double c0max = 2.0*pow(c1/3.0,1.5);

  double theta = acos(c0/c0max);
  u = sqrt(c1/3.0) * cos(theta/3.0);
  w = sqrt(c1) * sin(theta/3.0);
}
