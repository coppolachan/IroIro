/*!
 * @file fieldUtils.cpp
 * @brief function definitions in FieldUtils
 */
#include "fieldUtils.hpp"
#include "sunMatUtils.hpp"
#include "include/macros.hpp"
#include "timings.hpp"
#include "include/messages_macros.hpp"
#include <omp.h>

#ifdef IBM_BGQ_WILSON
#include "bgqthread.h"
#endif

#ifdef SR16K_WILSON
#include "srmwilson.h"
#include "srmwilson_cmpl.hpp" 
#endif

class RandNum;

namespace FieldUtils{
  const GaugeField1D field_oprod(const FermionField& f1,
				 const FermionField& f2){
    using namespace SUNmatUtils;

    GaugeField1D f;
    SUNmat mat;
    int Nd = CommonPrms::instance()->Nd();
    int Nvol = CommonPrms::instance()->Nvol();

    for(int site=0; site<Nvol; ++site){
      mat.zero();
      for(int s=0; s<Nd; ++s)
	mat += outer_prod_t(vec(f1,s,site), vec(f2,s,site));
      SetMat(f,mat,site);
    }
    return f;
  }

  const GaugeField ReUnit(const GaugeField& U){
    using namespace SUNmatUtils;
    GaugeField Ur(U);
#ifdef IBM_BGQ_WILSON
    register int Nvol = CommonPrms::instance()->Nvol();
#pragma omp parallel 
    {
      int nid = omp_get_num_threads();
      int tid = omp_get_thread_num();
      int is = tid*Nvol/nid;
      int ns = Nvol/nid;
      
      for(int mu=0; mu<U.Nex(); ++mu)
	for(int site=is; site<is+ns; ++site)
	  SetMat(Ur,reunit(mat(U,site,mu)),site,mu);
    }
#else
    for(int mu=0; mu<U.Nex(); ++mu)
      for(int site=0; site<U.Nvol(); ++site)
	SetMat(Ur,reunit(mat(U,site,mu)),site,mu);
#endif
    return Ur;
  }
  
  const GaugeField1D ReUnit(const GaugeField1D& G){
    using namespace SUNmatUtils;
    GaugeField1D Gr(G);
    for(int site=0; site<G.Nvol(); ++site)
      SetMat(Gr,reunit(mat(G,site)),site);
    return Gr;
  }

#if defined IBM_BGQ_WILSON  || defined SR16K_WILSON
  void reunit_ALIGNED(double* out, double* in) {
    double nrm_i;
    for(int a=0; a<NC_; ++a){
      nrm_i =0.0;
      for(int cc=0; cc<2*NC_; ++cc){
	out[a*2*NC_ +cc]=in[a*2*NC_ +cc];
	nrm_i += in[a*2*NC_ +cc]*in[a*2*NC_ +cc];
      }
      
      nrm_i = 1.0/sqrt(nrm_i);
      
      for(int cc=0; cc<2*NC_; cc++)
	out[a*2*NC_ +cc]*=nrm_i;

      for(int b=a+1; b<NC_; ++b){
	double prr = 0.0;
	double pri = 0.0;
	for(int c=0; c<NC_; ++c){
	  int ac = a*NC_+c;      
	  int bc = b*NC_+c;
	  prr += in[2*ac]*in[2*bc  ] +in[2*ac+1]*in[2*bc+1];
	  pri += in[2*ac]*in[2*bc+1] -in[2*ac+1]*in[2*bc  ];
	}
	for(int c=0; c<NC_; ++c){
	  int ac = a*NC_+c;      
	  int bc = b*NC_+c;
	  out[2*bc  ] = in[2*bc  ]- prr*in[2*ac  ] -pri*in[2*ac+1];
	  out[2*bc+1] = in[2*bc+1]- prr*in[2*ac+1] +pri*in[2*ac  ];
	}
      }
    }
  }
#endif 

#ifdef IBM_BGQ_WILSON  
  void ExponentiateMult_BGQ(GaugeField& U,const GaugeField& G,double d,int N){
    using namespace SUNmatUtils;
    GaugeField temp; //allocate fields here instead of creating objects
    GaugeField unit, temp2; //same here
    register int Nvol = G.Nvol();

    BGQThread_Init(); //initializing BGQ fast threading routines

    double* temp_ptr  = temp.data.getaddr(0);
    double* G_ptr  = const_cast<GaugeField&>(G).data.getaddr(0);
    double* U_ptr  = U.data.getaddr(0);
    double* temp2_ptr = temp2.data.getaddr(0);
    double* unit_ptr = unit.data.getaddr(0);
    
#pragma omp parallel
    {
      int nid = omp_get_num_threads();
      int tid = omp_get_thread_num();
      int ns = Nvol/nid;
      int ns2 = ns*G.Nex();
      int is = tid*ns;
      int Ncc = G.Nin();
      register int jump2 = tid*ns2*Ncc;

      BGWilsonSU3_MatUnity(unit_ptr+jump2,ns2);
      BGWilsonSU3_MatUnity(temp_ptr+jump2,ns2);
      
      //#pragma omp barrier   
      BGQThread_Barrier(0,nid);
      
      for (int i = N; i>=1;--i){
	BGWilsonSU3_MatMultScalar(temp_ptr+jump2,d/i,ns2);
	BGWilsonSU3_MatMultAdd_NN(temp2_ptr+jump2,unit_ptr+jump2, 
				  temp_ptr+jump2,G_ptr+jump2,ns2);
	BGWilsonSU3_MatEquate(temp_ptr+jump2,temp2_ptr+jump2,ns2);
      }

      //#pragma omp barrier   
      BGQThread_Barrier(0,nid);
      for(int mu=0; mu<G.Nex(); ++mu)
	for(int site=is; site<is+ns; ++site)
	  reunit_ALIGNED(temp_ptr+(site+Nvol*mu)*Ncc,
		       temp2_ptr+(site+Nvol*mu)*Ncc);
      //#pragma omp barrier  
      BGQThread_Barrier(0,nid);
      BGWilsonSU3_MatMult_NN(temp2_ptr+jump2,temp_ptr+jump2,U_ptr+jump2,ns2);

      //#pragma omp barrier          
      BGQThread_Barrier(0,nid);
      for(int mu=0; mu<G.Nex(); ++mu)
	for(int site=is; site<is+ns; ++site)
	  reunit_ALIGNED(U_ptr+(site+Nvol*mu)*Ncc,
		       temp2_ptr+(site+Nvol*mu)*Ncc);
    }
  } 
#elif defined SR16K_WILSON  
  void ExponentiateMult_SR16K(GaugeField& U,const GaugeField& G,double d,int N){

    using namespace SUNmatUtils;
    GaugeField temp; //allocate fields here instead of creating objects
    GaugeField unit, temp2; //same here
    register int Nvol = G.Nvol();

    double* temp_ptr  = temp.data.getaddr(0);
    double* G_ptr  = const_cast<GaugeField&>(G).data.getaddr(0);
    double* U_ptr  = U.data.getaddr(0);
    double* temp2_ptr = temp2.data.getaddr(0);
    double* unit_ptr = unit.data.getaddr(0);

    int ns2 = Nvol*G.Nex();
    int Ncc = G.Nin();

    SRCMPL::SRWilsonSU3_MatUnity(unit_ptr,ns2);
    SRCMPL::SRWilsonSU3_MatUnity(temp_ptr,ns2);

    for(int i=N; i>=1;--i){
      SRCMPL::SRWilsonSU3_MatMultScalar(temp_ptr,d/i,ns2);
      SRWilsonSU3_MatMultAdd_NN(temp2_ptr,unit_ptr,temp_ptr,G_ptr,ns2);
      SRCMPL::SRWilsonSU3_MatEquate(temp_ptr,temp2_ptr,ns2);
    }

    for(int nx=0; nx<ns2; ++nx)
      reunit_ALIGNED(temp_ptr+nx*Ncc,temp2_ptr+nx*Ncc);

    //temp = ReUnit(temp2);
    SRWilsonSU3_MatMult_NN(temp2_ptr,temp_ptr,U_ptr,ns2);

    for(int nx=0; nx<ns2; ++nx)
      reunit_ALIGNED(U_ptr+nx*Ncc,temp2_ptr+nx*Ncc);

    //U = ReUnit(temp2);
  } 
#endif

  const GaugeField Exponentiate(const GaugeField& G,double d,int N){
    using namespace SUNmatUtils;
    GaugeField temp;
#ifdef IBM_BGQ_WILSON

    GaugeField unit, temp2;
    register int Nvol = G.Nvol();

    for(int mu=0; mu<G.Nex(); ++mu)
      for(int site=0; site<Nvol; ++site)
        SetMat(unit, unity(),site,mu);
    
    temp = unit;
    double* temp_ptr  = temp.data.getaddr(0);
    double* G_ptr  = const_cast<GaugeField&>(G).data.getaddr(0);
    double* temp2_ptr = temp2.data.getaddr(0);
    double* unit_ptr = unit.data.getaddr(0);

    int ns = Nvol*G.Nex();
    for(int i=N; i>=1;--i){
      BGWilsonSU3_MatMultScalar(temp_ptr,d/i,ns);
      BGWilsonSU3_MatMultAdd_NN(temp2_ptr,unit_ptr,temp_ptr,G_ptr,ns);
      BGWilsonSU3_MatEquate(temp_ptr,temp2_ptr,ns);
    }
    return ReUnit(temp2);

#elif defined SR16K_WILSON
    CCIO::cout<<"FieldUtils::Exponenciate called\n";
    GaugeField unit, temp2;
    register int Nvol = G.Nvol();

    for(int mu=0; mu<G.Nex(); ++mu)
      for(int site=0; site<Nvol; ++site)
        SetMat(unit, unity(),site,mu);
    
    temp = unit;
    double* temp_ptr  = temp.data.getaddr(0);
    double* G_ptr  = const_cast<GaugeField&>(G).data.getaddr(0);
    double* temp2_ptr = temp2.data.getaddr(0);
    double* unit_ptr = unit.data.getaddr(0);

    int ns = Nvol*G.Nex();
    for(int i=N; i>=1;--i){
      SRCMPL::SRWilsonSU3_MatMultScalar(temp_ptr,d/i,ns);
      SRWilsonSU3_MatMultAdd_NN(temp2_ptr,unit_ptr,temp_ptr,G_ptr,ns);
      SRCMPL::SRWilsonSU3_MatEquate(temp_ptr,temp2_ptr,ns);
    }
    return ReUnit(temp2);

#else 
    for(int mu=0; mu<G.Nex(); ++mu)
      for(int site=0; site<G.Nvol(); ++site)
	SetMat(temp, exponential(mat(G,site,mu)*d,N),site,mu);

    return ReUnit(temp);

#endif
  }

  const GaugeField TracelessAntihermite(const GaugeField& G){
    using namespace SUNmatUtils;
    GaugeField TAField(G.Nvol());
    for(int mu=0; mu<G.Nex(); ++mu){
#pragma omp parallel for
      for(int site=0; site<G.Nvol(); ++site){
	SetMat(TAField,anti_hermite_traceless(mat(G,site,mu)),site,mu);
      }
    }
    return TAField;
  }

  const GaugeField1D TracelessAntihermite(const GaugeField1D& G){
    using namespace SUNmatUtils;
    GaugeField1D TAField(G.Nvol());
    for(int site=0; site<G.Nvol(); ++site)
      SetMat(TAField,anti_hermite_traceless(mat(G,site)),site);
    return TAField;
  }

#if defined IBM_BGQ_WILSON || defined SR16K_WILSON
  //assumes some ordering of the matrices
  void DirSliceALIGNED(GaugeField1D &G, const GaugeField& F, int dir){
    register int Nvol = F.Nvol();
    register int Ncc = F.Nin()/2;

    register double _Complex* pV0 = (double _Complex*)G.data.getaddr(0);
    register double _Complex* pW0 = (double _Complex*)const_cast<Field&>(F.data).getaddr(0)+ Ncc*Nvol*dir;

    for(int i=0;i<Nvol*Ncc; ++i) *(pV0 +i) = *(pW0 +i);
  }
#endif
  GaugeField1D DirSlice(const GaugeField& F, int dir){
    return GaugeField1D(Field(F.data[F.format.ex_slice(dir)]));
  }

  void SetSlice(GaugeField& G, const GaugeField1D& Gslice, int dir){
    G.data.set(G.format.ex_slice(dir), Gslice.data.getva());
  }

  void AddSlice(GaugeField& G, const GaugeField1D& Gslice, int dir){
    G.data.add(G.format.ex_slice(dir), Gslice.data.getva());
  }

  void SetMat(GaugeField& F, const SUNmat& mat, int site, int dir){
    F.data.set(F.format.islice(site,dir), mat.getva());
  }
  void SetMat(GaugeField1D& F, const SUNmat& mat, int site){
    F.data.set(F.format.islice(site), mat.getva());
  }

  void AddMat(GaugeField& F, const SUNmat& mat, int site, int dir){
    F.data.add(F.format.islice(site,dir), mat.getva());
  }
  void AddMat(GaugeField1D& F, const SUNmat& mat, int site){
    F.data.add(F.format.islice(site), mat.getva());
  }

  void SetVec(FermionField& F, const SUNvec& vec, int spin, int site){
    F.data.set(F.format.cslice(spin,site), vec.getva());
  }
  void SetVec(FermionField1sp& F, const SUNvec& vec, int site){
    F.data.set(F.format.islice(site), vec.getva());
  }

  void AddVec(FermionField& F, const SUNvec& vec, int spin, int site){
    F.data.add(F.format.cslice(spin,site), vec.getva());
  }
  void AddVec(FermionField1sp& F, const SUNvec& vec, int site){
    F.data.add(F.format.islice(site), vec.getva());
  }
}
