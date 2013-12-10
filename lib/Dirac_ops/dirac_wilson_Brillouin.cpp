/* !@filename dirac_wilson_Brillouin.cpp
 * @brief implementation of Dirac_Wilson_Brillouin class
 *  Time-stamp: <2013-12-10 17:32:06 noaki>
 */
#include "dirac_wilson_Brillouin.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include "Communicator/comm_io.hpp"
#include "Fields/field_expressions.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#include "Architecture_Optimized/dirac_wilson_Brillouin_BGQ.code"
#else
#include "dirac_wilson_Brillouin_regular.code"
#endif

using namespace SUNvecUtils;
using namespace std;

const Field Dirac_Wilson_Brillouin::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_Wilson_Brillouin::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

int Dirac_Wilson_Brillouin::sgm(int mu,int nu,int rho)const{
  for(int sigma=0; sigma<NDIM_; ++sigma){
    if(sigma !=mu && sigma !=nu && sigma !=rho) return sigma;
  }
}

const Field Dirac_Wilson_Brillouin::gamma5(const Field& f)const{ 
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site)
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
		   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}

void  Dirac_Wilson_Brillouin::Wsetup(const Field* u){ 
  //W_ = xxxx;
}
void  Dirac_Wilson_Brillouin::Fsetup(const Field* u){ 
  //F_ = yyyy;
}
