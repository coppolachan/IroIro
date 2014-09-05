/* !@filename dirac_wilson_Brillouin.cpp
 * @brief implementation of Dirac_Wilson_Brillouin class
 *  Time-stamp: <2013-09-26 15:57:51 noaki>
 */
/*
  In order to compile this file, one need the python script.
  The python script in a directory python.
  You have to execute
  ./python/bri_spinor_gammas.py > lib/Dirac_ops/dirac_wilson_Brillouin_spinor.code
  Should this be included into Makefile?

  Y-G.Cho
*/

#include "dirac_wilson_Brillouin_OSS.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include "Tools/fieldUtils.hpp"
#include "Communicator/comm_io.hpp"
#include "Fields/field_expressions.hpp"
#include "Geometry/shiftField.hpp"
#include <typeinfo>
#include <stdlib.h>
#include <omp.h>
#include "dirac_wilson_Brillouin_OSS_spinor.code"
#include <spi/include/kernel/process.h> //for GetTimeBase

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#endif

using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace SUNvecUtils;
using namespace std;
/*
const Field Dirac_Wilson_Brillouin_OSS::gammaX(const Field& f)const{
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site){
    gammaXcore(w.getaddr(ff_.index(0,site)),
	       const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  }
  return w;
}

const Field Dirac_Wilson_Brillouin_OSS::gammaY(const Field& f)const{ 
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site){
    gammaYcore(w.getaddr(ff_.index(0,site)),
	       const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  }
  return w;
}

const Field Dirac_Wilson_Brillouin_OSS::gammaZ(const Field& f)const{ 
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site){
    gammaZcore(w.getaddr(ff_.index(0,site)),
	       const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  }
  return w;
}

const Field Dirac_Wilson_Brillouin_OSS::gammaT(const Field& f)const{ 
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site){
    gammaTcore(w.getaddr(ff_.index(0,site)),
	       const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  }
  return w;
}

const Field Dirac_Wilson_Brillouin_OSS::gamma5(const Field& f)const{ 
  Field w(fsize_);
  for(int site=0;site<Nvol_;++site)
    gamma5core(w.getaddr(ff_.index(0,site)),
	       const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}
*/

const Field Dirac_Wilson_Brillouin_OSS::gamma5(const Field& f)const{ 
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site)
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
                   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}

void  Dirac_Wilson_Brillouin_OSS::Wsetup(const GaugeField& gu){ 
  CCIO::cout<<"start generating off-axis link variable."<<endl;
  for(int dir=0;dir<80;++dir){
    W_.set(Wfmt_.ex_slice(dir),((this->*hopterms[dir])(gu)).getva());
  }
  CCIO::cout<<"finished generating off-axis link variable."<<endl;
}

const Field Dirac_Wilson_Brillouin_OSS::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

void Dirac_Wilson_Brillouin_OSS::mult_std(Field& w,const Field& f)const{
  using namespace FieldExpression;
  w = mult_bri(f);
  w *= kbr_;
  w += f;
}

//mult for the improved Brillouin
void Dirac_Wilson_Brillouin_OSS::mult_imp(Field& w,const Field& f)const{
  using namespace FieldExpression;
  /*
  Field fn = lap_bri();

  w = lap_bri(fn);
  w -= 7.50*fn;
  w *= 0.125;

  fn -= 15.0/4.0*f; 
  fn *= -1.0/12.0;
  fn += f;

  Field dr = der_iso(fn); 

  fn = lap_bri(dr);
  fn /= -12.0;
  fn += (21.0/16.0)*dr;
  
  w += fn;
  w *= kbr_;
  w += f;
  */
}

const Field Dirac_Wilson_Brillouin_OSS::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

const Field Dirac_Wilson_Brillouin_OSS::mult_col(const Field& f,int dir)const{
  Field w(fsize_);
#ifdef IBM_BGQ_WILSON
  for(int s=0;s<ND_;++s){
    double* pU = const_cast<Field&>(W_).getaddr(Wfmt_.index(0,0,dir));
    double* pf = const_cast<Field&>(f).getaddr(ff_.index_r(0,s,0));
    double* pF = w.getaddr(ff_.index_r(0,s,0));
    BGWilsonSU3_MultU_1D_S(pF,ND_,pU,1,pf,ND_,Nvol_);
  }
#else
for(int site=0;site<Nvol;++site){
  const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,site));	  
  const double* U = const_cast<Field&>(W_).getaddr(Wfmt_.index(0,site,dir));
  for(int s=0; s<ND_; ++s){
    for(int c=0; c<NC_; ++c){
      double vr = 0.0, vi = 0.0;
      for(int c1=0; c1<NC_; ++c1){
	vr += U[re(c,c1)]*v[sr(s,c1)] -U[im(c,c1)]*v[si(s,c1)];
	vi += U[im(c,c1)]*v[sr(s,c1)] +U[re(c,c1)]*v[si(s,c1)];
      }
      w.add(ff_.index_r(c,s,site),vr);
      w.add(ff_.index_i(c,s,site),vi);
    }
  }
 }
#endif
 return w;
}

void Dirac_Wilson_Brillouin_OSS::spinor_col(Field& w,const Field& f,int dir)const{
  using namespace FieldExpression;
  Field tmp  = (this->*spinor[dir])(mult_col(f,dir));
  double* pf = tmp.getaddr(0);
  double* pw = w.getaddr(0);
  BGWilsonLA_Add(pw,pf,Nvol_);
  //double* pf = const_cast<Field>((this->*spinor[dir])(mult_col(f,dir))).getaddr(0);
  //BGWilsonLA_Add(pw,const_cast<Field>((this->*spinor[dir])(mult_col(f,dir))).getaddr(0),Nvol_);

}

/*
const Field Dirac_Wilson_Brillouin_OSS::spinor_col(const Field& f,int dir)const{
  using namespace FieldExpression;
  Field w(fsize_);
  w = (this->*spinor[dir])(mult_col(f,dir));
  return w;
}
*/

// shift field
void Dirac_Wilson_Brillouin_OSS::shift_field(Field &w, const Field& f,int dir,int ifb) const{
  double* pf = const_cast<Field&>(f).getaddr(0);
  double* pw = w.getaddr(0);
  BGWilson_Shift_Dir(pw,pf,Nin_*sizeof(double),dir,ifb);
}

void Dirac_Wilson_Brillouin_OSS::shift_fwd(Field& w,const Field& f,int dir) const{
  double* pf = const_cast<Field&>(f).getaddr(0);
  double* pw = w.getaddr(0);
  BGWilson_Shift_Dir(pw,pf,Nin_*sizeof(double),dir,BGWILSON_FORWARD);
}

void Dirac_Wilson_Brillouin_OSS::shift_bwd(Field& w,const Field& f,int dir) const{
  double* pf = const_cast<Field&>(f).getaddr(0);
  double* pw = w.getaddr(0);
  BGWilson_Shift_Dir(pw,pf,Nin_*sizeof(double),dir,BGWILSON_BACKWARD);
}

// Sum(mu)(Gamma(mu)*nab^{iso}_{mu}) - 0.5*lap^{bri} + 0.5 * a 0hop term
const Field Dirac_Wilson_Brillouin_OSS::mult_bri(const Field& f)const{ 
  Field w(fsize_);
  //int time1 = GetTimeBase();
  {// (+,0,0,0)
    Field fpX(fsize_);
    shift_fwd(fpX,f,XDIR);  spinor_col(w,fpX,0); 
    {// (+,+,0,0)
      Field fpXpY(fsize_);
      shift_fwd(fpXpY,fpX,YDIR);  spinor_col(w,fpXpY,8); 
      {// (+,+,+,0)
	Field fpXpYpZ(fsize_);
	shift_fwd(fpXpYpZ,fpXpY,ZDIR);  spinor_col(w,fpXpYpZ,32);
	{// (+,+,+,+)
	  Field fpXpYpZpT(fsize_);
	  shift_fwd(fpXpYpZpT,fpXpYpZ,TDIR);  spinor_col(w,fpXpYpZpT,64);
	}
	{// (+,+,+,-)
	  Field fpXpYpZmT(fsize_);
	  shift_bwd(fpXpYpZmT,fpXpYpZ,TDIR);  spinor_col(w,fpXpYpZmT,66);  
	}
      }
      {// (+,+,0,+)
	Field fpXpYpT(fsize_);
	shift_fwd(fpXpYpT,fpXpY,TDIR);	spinor_col(w,fpXpYpT,34);
      }
      {// (+,+,-,0)
	Field fpXpYmZ(fsize_);
	shift_bwd(fpXpYmZ,fpXpY,ZDIR);	spinor_col(w,fpXpYmZ,40);
	{// (+,+,-,+)
	  Field fpXpYmZpT(fsize_);
	  shift_fwd(fpXpYmZpT,fpXpYmZ,TDIR);  spinor_col(w,fpXpYmZpT,68); 
	}
	{// (+,+,-,-)
	  Field fpXpYmZmT(fsize_);
	  shift_bwd(fpXpYmZmT,fpXpYmZ,TDIR);  spinor_col(w,fpXpYmZmT,72);
	}
      }
      {// (+,+,0,-)
	Field fpXpYmT(fsize_);
	shift_bwd(fpXpYmT,fpXpY,TDIR);	spinor_col(w,fpXpYmT,42);
      }
    }
    {// (+,-,0,0)
      Field fpXmY(fsize_);
      shift_bwd(fpXmY,fpX,YDIR);  spinor_col(w,fpXmY,10);
      {// (+,-,+,0)
	Field fpXmYpZ(fsize_);
	shift_fwd(fpXmYpZ,fpXmY,ZDIR);	spinor_col(w,fpXmYpZ,48); 
	{// (+,-,+,+)
	  Field fpXmYpZpT(fsize_);
	  shift_fwd(fpXmYpZpT,fpXmYpZ,TDIR);  spinor_col(w,fpXmYpZpT,70); 
	}
	{// (+,-,+,-)
	  Field fpXmYpZmT(fsize_);
	  shift_bwd(fpXmYpZmT,fpXmYpZ,TDIR);  spinor_col(w,fpXmYpZmT,74);
	}
      }
      {// (+,-,0,+)
	Field fpXmYpT(fsize_);
	shift_fwd(fpXmYpT,fpXmY,TDIR); 	spinor_col(w,fpXmYpT,50);
      }
      {// (+,-,-,0)
	Field fpXmYmZ(fsize_);
	shift_bwd(fpXmYmZ,fpXmY,ZDIR); 	spinor_col(w,fpXmYmZ,56); 
	{
	  Field fpXmYmZpT(fsize_);
	  shift_fwd(fpXmYmZpT,fpXmYmZ,TDIR);  spinor_col(w,fpXmYmZpT,76); 
	}
	{
	  Field fpXmYmZmT(fsize_);
	  shift_bwd(fpXmYmZmT,fpXmYmZ,TDIR);  spinor_col(w,fpXmYmZmT,78); 
	}
      }
      {// (+,-,0,-)
	Field fpXmYmT(fsize_);
	shift_bwd(fpXmYmT,fpXmY,TDIR); 	spinor_col(w,fpXmYmT,58); 
      }
    }
    {// (+,0,+,0) 
      Field fpXpZ(fsize_);
      shift_fwd(fpXpZ,fpX,ZDIR);  spinor_col(w,fpXpZ,12);
      {// (+,0,+,+) 
	Field fpXpZpT(fsize_);
	shift_fwd(fpXpZpT,fpXpZ,TDIR); 	spinor_col(w,fpXpZpT,36);
      }
      {// (+,0,+,-) 
	Field fpXpZmT(fsize_);
	shift_bwd(fpXpZmT,fpXpZ,TDIR); 	spinor_col(w,fpXpZmT,44); 
      }
    }
    { // (+,0,-,0)
      Field fpXmZ(fsize_);
      shift_bwd(fpXmZ,fpX,ZDIR);  spinor_col(w,fpXmZ,14);
      {// (+,0,-,+)
	Field fpXmZpT(fsize_);
	shift_fwd(fpXmZpT,fpXmZ,TDIR); 	spinor_col(w,fpXmZpT,52);
      }
      {// (+,0,-,-)
	Field fpXmZmT(fsize_);
	shift_bwd(fpXmZmT,fpXmZ,TDIR); 	spinor_col(w,fpXmZmT,60);
      }
    }
    { // (+,0,0,+)
      Field fpXpT(fsize_);
      shift_fwd(fpXpT,fpX,TDIR);  spinor_col(w,fpXpT,16);
    }
    { // (+,0,0,-)
      Field fpXmT(fsize_);
      shift_bwd(fpXmT,fpX,TDIR);  spinor_col(w,fpXmT,18);
    }
  }
  {// (-,0,0,0) 
    Field fmX(fsize_);    
    shift_bwd(fmX,f,XDIR);  spinor_col(w,fmX,1);
    
    {// (-,+,0,0) 
      Field fmXpY(fsize_);
      shift_fwd(fmXpY,fmX,YDIR);  spinor_col(w,fmXpY,9);
      {// (-,+,0,+) 
	Field fmXpYpT(fsize_);
	shift_fwd(fmXpYpT,fmXpY,TDIR);	spinor_col(w,fmXpYpT,35);
      }
      {// (-,+,+,0) 
	Field fmXpYpZ(fsize_);
	shift_fwd(fmXpYpZ,fmXpY,ZDIR);	spinor_col(w,fmXpYpZ,33);
	{// (-,+,+,+)
	  Field fmXpYpZpT(fsize_);
	  shift_fwd(fmXpYpZpT,fmXpYpZ,TDIR);  spinor_col(w,fmXpYpZpT,65);
	}
	{// (-,+,+,-)
	  Field fmXpYpZmT(fsize_);
	  shift_bwd(fmXpYpZmT,fmXpYpZ,TDIR);  spinor_col(w,fmXpYpZmT,67);
	}
      }
      {// (-,+,-,0) 
	Field fmXpYmZ(fsize_);	
	shift_bwd(fmXpYmZ,fmXpY,ZDIR);	spinor_col(w,fmXpYmZ,41);
	{// (-,+,-,+)
	  Field fmXpYmZpT(fsize_);	
	  shift_fwd(fmXpYmZpT,fmXpYmZ,TDIR);  spinor_col(w,fmXpYmZpT,69);
	}
	{// (-,+,-,-)
	  Field fmXpYmZmT(fsize_);	
	  shift_bwd(fmXpYmZmT,fmXpYmZ,TDIR);  spinor_col(w,fmXpYmZmT,73);
	}
      }
      {// (-,+,0,-)
	Field fmXpYmT(fsize_);	
	shift_bwd(fmXpYmT,fmXpY,TDIR);	spinor_col(w,fmXpYmT,43);
      }
    }
    {// (-,-,0,0)
      Field fmXmY(fsize_);
      shift_bwd(fmXmY,fmX,YDIR);  spinor_col(w,fmXmY,11);
      {// (-,-,+,0)
	Field fmXmYpZ(fsize_);
	shift_fwd(fmXmYpZ,fmXmY,ZDIR);	spinor_col(w,fmXmYpZ,49);
	{// (-,-,+,+)
	  Field fmXmYpZpT(fsize_);
	  shift_fwd(fmXmYpZpT,fmXmYpZ,TDIR);  spinor_col(w,fmXmYpZpT,71);
	}
	{// (-,-,+,-)
	  Field fmXmYpZmT(fsize_);
	  shift_bwd(fmXmYpZmT,fmXmYpZ,TDIR);  spinor_col(w,fmXmYpZmT,75);
	}
      }
      {// (-,-,-,0)
	Field fmXmYmZ(fsize_);
	shift_bwd(fmXmYmZ,fmXmY,ZDIR);	spinor_col(w,fmXmYmZ,57);
	{// (-,-,-,+)
	  Field fmXmYmZpT(fsize_);
	  shift_fwd(fmXmYmZpT,fmXmYmZ,TDIR);  spinor_col(w,fmXmYmZpT,77);
	}
	{// (-,-,-,-)
	  Field fmXmYmZmT(fsize_);
	  shift_bwd(fmXmYmZmT,fmXmYmZ,TDIR);  spinor_col(w,fmXmYmZmT,79);
	}
      }
      {// (-,-,0,+)
	Field fmXmYpT(fsize_);
	shift_fwd(fmXmYpT,fmXmY,TDIR);	spinor_col(w,fmXmYpT,51);
      }
      {// (-,-,0,-)
	Field fmXmYmT(fsize_);
	shift_bwd(fmXmYmT,fmXmY,TDIR);	spinor_col(w,fmXmYmT,59);
      }
    }
    {// (-,0,+,0)
      Field fmXpZ(fsize_);
      shift_fwd(fmXpZ,fmX,ZDIR);  spinor_col(w,fmXpZ,13);

      {// (-,0,+,+)
	Field fmXpZpT(fsize_);
	shift_fwd(fmXpZpT,fmXpZ,TDIR); spinor_col(w,fmXpZpT,37);
      }
      {// (-,0,+,-)
	Field fmXpZmT(fsize_);
	shift_bwd(fmXpZmT,fmXpZ,TDIR);	spinor_col(w,fmXpZmT,45);
      }
    }
    {// (-,0,-,0)
      Field fmXmZ(fsize_);
      shift_bwd(fmXmZ,fmX,ZDIR);  spinor_col(w,fmXmZ,15);
      {// (-,0,-,+)
	Field fmXmZpT(fsize_);
	shift_fwd(fmXmZpT,fmXmZ,TDIR);	spinor_col(w,fmXmZpT,53);
      }
      {// (-,0,-,-)
	Field fmXmZmT(fsize_);
	shift_bwd(fmXmZmT,fmXmZ,TDIR);	spinor_col(w,fmXmZmT,61);
      }
    }
    {// (-,0,0,+)
      Field fmXpT(fsize_);
      shift_fwd(fmXpT,fmX,TDIR);  spinor_col(w,fmXpT,17);
    }
    {// (-,0,0,-)
      Field fmXmT(fsize_);
      shift_bwd(fmXmT,fmX,TDIR);  spinor_col(w,fmXmT,19);
    }
  }
  {// (0,+,0,0)
    Field fpY(fsize_);
    shift_fwd(fpY,f,YDIR);  spinor_col(w,fpY,2); 
    {// (0,+,+,0)
      Field fpYpZ(fsize_);
      shift_fwd(fpYpZ,fpY,ZDIR); spinor_col(w,fpYpZ,20);
      {// (0,+,+,+)
	Field fpYpZpT(fsize_);
	shift_fwd(fpYpZpT,fpYpZ,TDIR); 	spinor_col(w,fpYpZpT,38);
      }
      {// (0,+,+,-)
	Field fpYpZmT(fsize_);
	shift_bwd(fpYpZmT,fpYpZ,TDIR); 	spinor_col(w,fpYpZmT,46);
      }
    }
    {// (0,+,-,0)
      Field fpYmZ(fsize_);
      shift_bwd(fpYmZ,fpY,ZDIR); spinor_col(w,fpYmZ,22);
      {// (0,+,-,+)
	Field fpYmZpT(fsize_);
	shift_fwd(fpYmZpT,fpYmZ,TDIR); spinor_col(w,fpYmZpT,54);
      }	
      {// (0,+,-,-)
	Field fpYmZmT(fsize_);
	shift_bwd(fpYmZmT,fpYmZ,TDIR);	spinor_col(w,fpYmZmT,62);
      }
    }
    {// (0,+,0,+)
      Field fpYpT(fsize_);
      shift_fwd(fpYpT,fpY,TDIR);  spinor_col(w,fpYpT,24);
    }
    {// (0,+,0,-)
      Field fpYmT(fsize_);
      shift_bwd(fpYmT,fpY,TDIR);  spinor_col(w,fpYmT,26);
    }
  }
  {// (0,-,0,0)
    Field fmY(fsize_);
    shift_bwd(fmY,f,YDIR);  spinor_col(w,fmY,3);
    {// (0,-,+,0)
      Field fmYpZ(fsize_);
      shift_fwd(fmYpZ,fmY,ZDIR);  spinor_col(w,fmYpZ,21);
      {// (0,-,+,+)
	Field fmYpZpT(fsize_);
	shift_fwd(fmYpZpT,fmYpZ,TDIR);	spinor_col(w,fmYpZpT,39);
      }
      {// (0,-,+,-)
	Field fmYpZmT(fsize_);
	shift_bwd(fmYpZmT,fmYpZ,TDIR);	spinor_col(w,fmYpZmT,47);
      }
    }
    {// (0,-,-,0)
      Field fmYmZ(fsize_);
      shift_bwd(fmYmZ,fmY,ZDIR);  spinor_col(w,fmYmZ,23);
      {// (0,-,-,+)
	Field fmYmZpT(fsize_);
	shift_fwd(fmYmZpT,fmYmZ,TDIR);	spinor_col(w,fmYmZpT,55);
      }
      {// (0,-,-,-)
	Field fmYmZmT(fsize_);
	shift_bwd(fmYmZmT,fmYmZ,TDIR);	spinor_col(w,fmYmZmT,63);
      }
    }
    {// (0,-,0,+)
      Field fmYpT(fsize_);
      shift_fwd(fmYpT,fmY,TDIR);  spinor_col(w,fmYpT,25);
    }
    {// (0,-,0,-)
      Field fmYmT(fsize_);
      shift_bwd(fmYmT,fmY,TDIR);  spinor_col(w,fmYmT,27); 
    }
  }
  {// (0,0,+,0)
    Field fpZ(fsize_);
    shift_fwd(fpZ,f,ZDIR);  spinor_col(w,fpZ,4);    

    {// (0,0,+,+)
      Field fpZpT(fsize_);
      shift_fwd(fpZpT,fpZ,TDIR);  spinor_col(w,fpZpT,28);
    }
    {// (0,0,+,-)
      Field fpZmT(fsize_);
      shift_bwd(fpZmT,fpZ,TDIR);  spinor_col(w,fpZmT,30);
    }
  }
  {// (0,0,-,0)
    Field fmZ(fsize_);
    shift_bwd(fmZ,f,ZDIR);  spinor_col(w,fmZ,5);

    {// (0,0,-,+)
      Field fmZpT(fsize_);
      shift_fwd(fmZpT,fmZ,TDIR);  spinor_col(w,fmZpT,29);
    }
    {// (0,0,-,-)
      Field fmZmT(fsize_);      
      shift_bwd(fmZmT,fmZ,TDIR);  spinor_col(w,fmZmT,31);
    }    
  }
  {// (0,0,0,+)
    Field fpT(fsize_);
    shift_fwd(fpT,f,TDIR); spinor_col(w,fpT,6);
  }
  {// (0,0,0,-)
    Field fmT(fsize_);
    shift_bwd(fmT,f,TDIR); spinor_col(w,fmT,7);
  }
  //int time2 = GetTimeBase();
  //CCIO::cout<<"time of mutl_bri : "<<time2-time1<<endl;
  return w;
}


/*
const Field Dirac_Wilson_Brillouin_OSS::mult_bri(const Field& f)const{ 
  Field w(fsize_);

  int time1 = GetTimeBase();
  {// (+,0,0,0)
    Field fpX(fsize_);
    shift_fwd(fpX,f,XDIR);  w += spinor_col(fpX,0); 
    {// (+,+,0,0)
      Field fpXpY(fsize_);
      shift_fwd(fpXpY,fpX,YDIR);  w += spinor_col(fpXpY,8); 
      {// (+,+,+,0)
	Field fpXpYpZ(fsize_);
	shift_fwd(fpXpYpZ,fpXpY,ZDIR);  w += spinor_col(fpXpYpZ,32);
	{// (+,+,+,+)
	  Field fpXpYpZpT(fsize_);
	  shift_fwd(fpXpYpZpT,fpXpYpZ,TDIR);  w += spinor_col(fpXpYpZpT,64);
	}
	{// (+,+,+,-)
	  Field fpXpYpZmT(fsize_);
	  shift_bwd(fpXpYpZmT,fpXpYpZ,TDIR);  w += spinor_col(fpXpYpZmT,66);  
	}
      }
      {// (+,+,0,+)
	Field fpXpYpT(fsize_);
	shift_fwd(fpXpYpT,fpXpY,TDIR);	w += spinor_col(fpXpYpT,34);
      }
      {// (+,+,-,0)
	Field fpXpYmZ(fsize_);
	shift_bwd(fpXpYmZ,fpXpY,ZDIR);	w += spinor_col(fpXpYmZ,40);
	{// (+,+,-,+)
	  Field fpXpYmZpT(fsize_);
	  shift_fwd(fpXpYmZpT,fpXpYmZ,TDIR);  w += spinor_col(fpXpYmZpT,68); 
	}
	{// (+,+,-,-)
	  Field fpXpYmZmT(fsize_);
	  shift_bwd(fpXpYmZmT,fpXpYmZ,TDIR);  w += spinor_col(fpXpYmZmT,72);
	}
      }
      {// (+,+,0,-)
	Field fpXpYmT(fsize_);
	shift_bwd(fpXpYmT,fpXpY,TDIR);	w += spinor_col(fpXpYmT,42);
      }
    }
    {// (+,-,0,0)
      Field fpXmY(fsize_);
      shift_bwd(fpXmY,fpX,YDIR);  w += spinor_col(fpXmY,10);
      {// (+,-,+,0)
	Field fpXmYpZ(fsize_);
	shift_fwd(fpXmYpZ,fpXmY,ZDIR);	w += spinor_col(fpXmYpZ,48); 
	{// (+,-,+,+)
	  Field fpXmYpZpT(fsize_);
	  shift_fwd(fpXmYpZpT,fpXmYpZ,TDIR);  w += spinor_col(fpXmYpZpT,70); 
	}
	{// (+,-,+,-)
	  Field fpXmYpZmT(fsize_);
	  shift_bwd(fpXmYpZmT,fpXmYpZ,TDIR);  w += spinor_col(fpXmYpZmT,74);
	}
      }
      {// (+,-,0,+)
	Field fpXmYpT(fsize_);
	shift_fwd(fpXmYpT,fpXmY,TDIR); 	w += spinor_col(fpXmYpT,50);
      }
      {// (+,-,-,0)
	Field fpXmYmZ(fsize_);
	shift_bwd(fpXmYmZ,fpXmY,ZDIR); 	w += spinor_col(fpXmYmZ,56); 
	{
	  Field fpXmYmZpT(fsize_);
	  shift_fwd(fpXmYmZpT,fpXmYmZ,TDIR);  w += spinor_col(fpXmYmZpT,76); 
	}
	{
	  Field fpXmYmZmT(fsize_);
	  shift_bwd(fpXmYmZmT,fpXmYmZ,TDIR);  w += spinor_col(fpXmYmZmT,78); 
	}
      }
      {// (+,-,0,-)
	Field fpXmYmT(fsize_);
	shift_bwd(fpXmYmT,fpXmY,TDIR); 	w += spinor_col(fpXmYmT,58); 
      }
    }
    {// (+,0,+,0) 
      Field fpXpZ(fsize_);
      shift_fwd(fpXpZ,fpX,ZDIR);  w += spinor_col(fpXpZ,12);
      {// (+,0,+,+) 
	Field fpXpZpT(fsize_);
	shift_fwd(fpXpZpT,fpXpZ,TDIR); 	w += spinor_col(fpXpZpT,36);
      }
      {// (+,0,+,-) 
	Field fpXpZmT(fsize_);
	shift_bwd(fpXpZmT,fpXpZ,TDIR); 	w += spinor_col(fpXpZmT,44); 
      }
    }
    { // (+,0,-,0)
      Field fpXmZ(fsize_);
      shift_bwd(fpXmZ,fpX,ZDIR);  w += spinor_col(fpXmZ,14);
      {// (+,0,-,+)
	Field fpXmZpT(fsize_);
	shift_fwd(fpXmZpT,fpXmZ,TDIR); 	w += spinor_col(fpXmZpT,52);
      }
      {// (+,0,-,-)
	Field fpXmZmT(fsize_);
	shift_bwd(fpXmZmT,fpXmZ,TDIR); 	w += spinor_col(fpXmZmT,60);
      }
    }
    { // (+,0,0,+)
      Field fpXpT(fsize_);
      shift_fwd(fpXpT,fpX,TDIR);  w += spinor_col(fpXpT,16);
    }
    { // (+,0,0,-)
      Field fpXmT(fsize_);
      shift_bwd(fpXmT,fpX,TDIR);  w += spinor_col(fpXmT,18);
    }
  }
  {// (-,0,0,0) 
    Field fmX(fsize_);    
    shift_bwd(fmX,f,XDIR);  w += spinor_col(fmX,1);
    
    {// (-,+,0,0) 
      Field fmXpY(fsize_);
      shift_fwd(fmXpY,fmX,YDIR);  w += spinor_col(fmXpY,9);
      {// (-,+,0,+) 
	Field fmXpYpT(fsize_);
	shift_fwd(fmXpYpT,fmXpY,TDIR);	w += spinor_col(fmXpYpT,35);
      }
      {// (-,+,+,0) 
	Field fmXpYpZ(fsize_);
	shift_fwd(fmXpYpZ,fmXpY,ZDIR);	w += spinor_col(fmXpYpZ,33);
	{// (-,+,+,+)
	  Field fmXpYpZpT(fsize_);
	  shift_fwd(fmXpYpZpT,fmXpYpZ,TDIR);  w += spinor_col(fmXpYpZpT,65);
	}
	{// (-,+,+,-)
	  Field fmXpYpZmT(fsize_);
	  shift_bwd(fmXpYpZmT,fmXpYpZ,TDIR);  w += spinor_col(fmXpYpZmT,67);
	}
      }
      {// (-,+,-,0) 
	Field fmXpYmZ(fsize_);	
	shift_bwd(fmXpYmZ,fmXpY,ZDIR);	w += spinor_col(fmXpYmZ,41);
	{// (-,+,-,+)
	  Field fmXpYmZpT(fsize_);	
	  shift_fwd(fmXpYmZpT,fmXpYmZ,TDIR);  w += spinor_col(fmXpYmZpT,69);
	}
	{// (-,+,-,-)
	  Field fmXpYmZmT(fsize_);	
	  shift_bwd(fmXpYmZmT,fmXpYmZ,TDIR);  w += spinor_col(fmXpYmZmT,73);
	}
      }
      {// (-,+,0,-)
	Field fmXpYmT(fsize_);	
	shift_bwd(fmXpYmT,fmXpY,TDIR);	w += spinor_col(fmXpYmT,43);
      }
    }
    {// (-,-,0,0)
      Field fmXmY(fsize_);
      shift_bwd(fmXmY,fmX,YDIR);  w += spinor_col(fmXmY,11);
      {// (-,-,+,0)
	Field fmXmYpZ(fsize_);
	shift_fwd(fmXmYpZ,fmXmY,ZDIR);	w += spinor_col(fmXmYpZ,49);
	{// (-,-,+,+)
	  Field fmXmYpZpT(fsize_);
	  shift_fwd(fmXmYpZpT,fmXmYpZ,TDIR);  w += spinor_col(fmXmYpZpT,71);
	}
	{// (-,-,+,-)
	  Field fmXmYpZmT(fsize_);
	  shift_bwd(fmXmYpZmT,fmXmYpZ,TDIR);  w += spinor_col(fmXmYpZmT,75);
	}
      }
      {// (-,-,-,0)
	Field fmXmYmZ(fsize_);
	shift_bwd(fmXmYmZ,fmXmY,ZDIR);	w += spinor_col(fmXmYmZ,57);
	{// (-,-,-,+)
	  Field fmXmYmZpT(fsize_);
	  shift_fwd(fmXmYmZpT,fmXmYmZ,TDIR);  w += spinor_col(fmXmYmZpT,77);
	}
	{// (-,-,-,-)
	  Field fmXmYmZmT(fsize_);
	  shift_bwd(fmXmYmZmT,fmXmYmZ,TDIR);  w += spinor_col(fmXmYmZmT,79);
	}
      }
      {// (-,-,0,+)
	Field fmXmYpT(fsize_);
	shift_fwd(fmXmYpT,fmXmY,TDIR);	w += spinor_col(fmXmYpT,51);
      }
      {// (-,-,0,-)
	Field fmXmYmT(fsize_);
	shift_bwd(fmXmYmT,fmXmY,TDIR);	w += spinor_col(fmXmYmT,59);
      }
    }
    {// (-,0,+,0)
      Field fmXpZ(fsize_);
      shift_fwd(fmXpZ,fmX,ZDIR);  w += spinor_col(fmXpZ,13);

      {// (-,0,+,+)
	Field fmXpZpT(fsize_);
	shift_fwd(fmXpZpT,fmXpZ,TDIR); w += spinor_col(fmXpZpT,37);
      }
      {// (-,0,+,-)
	Field fmXpZmT(fsize_);
	shift_bwd(fmXpZmT,fmXpZ,TDIR);	w += spinor_col(fmXpZmT,45);
      }
    }
    {// (-,0,-,0)
      Field fmXmZ(fsize_);
      shift_bwd(fmXmZ,fmX,ZDIR);  w += spinor_col(fmXmZ,15);
      {// (-,0,-,+)
	Field fmXmZpT(fsize_);
	shift_fwd(fmXmZpT,fmXmZ,TDIR);	w += spinor_col(fmXmZpT,53);
      }
      {// (-,0,-,-)
	Field fmXmZmT(fsize_);
	shift_bwd(fmXmZmT,fmXmZ,TDIR);	w += spinor_col(fmXmZmT,61);
      }
    }
    {// (-,0,0,+)
      Field fmXpT(fsize_);
      shift_fwd(fmXpT,fmX,TDIR);  w += spinor_col(fmXpT,17);
    }
    {// (-,0,0,-)
      Field fmXmT(fsize_);
      shift_bwd(fmXmT,fmX,TDIR);  w += spinor_col(fmXmT,19);
    }
  }
  {// (0,+,0,0)
    Field fpY(fsize_);
    shift_fwd(fpY,f,YDIR);  w += spinor_col(fpY,2); 
    {// (0,+,+,0)
      Field fpYpZ(fsize_);
      shift_fwd(fpYpZ,fpY,ZDIR); w += spinor_col(fpYpZ,20);
      {// (0,+,+,+)
	Field fpYpZpT(fsize_);
	shift_fwd(fpYpZpT,fpYpZ,TDIR); 	w += spinor_col(fpYpZpT,38);
      }
      {// (0,+,+,-)
	Field fpYpZmT(fsize_);
	shift_bwd(fpYpZmT,fpYpZ,TDIR); 	w += spinor_col(fpYpZmT,46);
      }
    }
    {// (0,+,-,0)
      Field fpYmZ(fsize_);
      shift_bwd(fpYmZ,fpY,ZDIR); w += spinor_col(fpYmZ,22);
      {// (0,+,-,+)
	Field fpYmZpT(fsize_);
	shift_fwd(fpYmZpT,fpYmZ,TDIR); w += spinor_col(fpYmZpT,54);
      }	
      {// (0,+,-,-)
	Field fpYmZmT(fsize_);
	shift_bwd(fpYmZmT,fpYmZ,TDIR);	w += spinor_col(fpYmZmT,62);
      }
    }
    {// (0,+,0,+)
      Field fpYpT(fsize_);
      shift_fwd(fpYpT,fpY,TDIR);  w += spinor_col(fpYpT,24);
    }
    {// (0,+,0,-)
      Field fpYmT(fsize_);
      shift_bwd(fpYmT,fpY,TDIR);  w += spinor_col(fpYmT,26);
    }
  }
  {// (0,-,0,0)
    Field fmY(fsize_);
    shift_bwd(fmY,f,YDIR);  w += spinor_col(fmY,3);
    {// (0,-,+,0)
      Field fmYpZ(fsize_);
      shift_fwd(fmYpZ,fmY,ZDIR);  w += spinor_col(fmYpZ,21);
      {// (0,-,+,+)
	Field fmYpZpT(fsize_);
	shift_fwd(fmYpZpT,fmYpZ,TDIR);	w += spinor_col(fmYpZpT,39);
      }
      {// (0,-,+,-)
	Field fmYpZmT(fsize_);
	shift_bwd(fmYpZmT,fmYpZ,TDIR);	w += spinor_col(fmYpZmT,47);
      }
    }
    {// (0,-,-,0)
      Field fmYmZ(fsize_);
      shift_bwd(fmYmZ,fmY,ZDIR);  w += spinor_col(fmYmZ,23);
      {// (0,-,-,+)
	Field fmYmZpT(fsize_);
	shift_fwd(fmYmZpT,fmYmZ,TDIR);	w += spinor_col(fmYmZpT,55);
      }
      {// (0,-,-,-)
	Field fmYmZmT(fsize_);
	shift_bwd(fmYmZmT,fmYmZ,TDIR);	w += spinor_col(fmYmZmT,63);
      }
    }
    {// (0,-,0,+)
      Field fmYpT(fsize_);
      shift_fwd(fmYpT,fmY,TDIR);  w += spinor_col(fmYpT,25);
    }
    {// (0,-,0,-)
      Field fmYmT(fsize_);
      shift_bwd(fmYmT,fmY,TDIR);  w += spinor_col(fmYmT,27); 
    }
  }
  {// (0,0,+,0)
    Field fpZ(fsize_);
    shift_fwd(fpZ,f,ZDIR);  w += spinor_col(fpZ,4);    

    {// (0,0,+,+)
      Field fpZpT(fsize_);
      shift_fwd(fpZpT,fpZ,TDIR);  w += spinor_col(fpZpT,28);
    }
    {// (0,0,+,-)
      Field fpZmT(fsize_);
      shift_bwd(fpZmT,fpZ,TDIR);  w += spinor_col(fpZmT,30);
    }
  }
  {// (0,0,-,0)
    Field fmZ(fsize_);
    shift_bwd(fmZ,f,ZDIR);  w += spinor_col(fmZ,5);

    {// (0,0,-,+)
      Field fmZpT(fsize_);
      shift_fwd(fmZpT,fmZ,TDIR);  w += spinor_col(fmZpT,29);
    }
    {// (0,0,-,-)
      Field fmZmT(fsize_);      
      shift_bwd(fmZmT,fmZ,TDIR);  w += spinor_col(fmZmT,31);
    }    
  }
  {// (0,0,0,+)
    Field fpT(fsize_);
    shift_fwd(fpT,f,TDIR); w += spinor_col(fpT,6);
  }
  {// (0,0,0,-)
    Field fmT(fsize_);
    shift_bwd(fmT,f,TDIR); w += spinor_col(fmT,7);
  }
  int time2 = GetTimeBase();
  CCIO::cout<<"time of mutl_bri : "<<time2-time1<<endl;
  return w;
}
*/

const Field (Dirac_Wilson_Brillouin_OSS::*Dirac_Wilson_Brillouin_OSS::spinor[])(const Field&)const
={
  &Dirac_Wilson_Brillouin_OSS::spinor0, &Dirac_Wilson_Brillouin_OSS::spinor1,
  &Dirac_Wilson_Brillouin_OSS::spinor2, &Dirac_Wilson_Brillouin_OSS::spinor3,
  &Dirac_Wilson_Brillouin_OSS::spinor4, &Dirac_Wilson_Brillouin_OSS::spinor5,
  &Dirac_Wilson_Brillouin_OSS::spinor6, &Dirac_Wilson_Brillouin_OSS::spinor7,
  &Dirac_Wilson_Brillouin_OSS::spinor8, &Dirac_Wilson_Brillouin_OSS::spinor9,
  &Dirac_Wilson_Brillouin_OSS::spinor10, &Dirac_Wilson_Brillouin_OSS::spinor11,
  &Dirac_Wilson_Brillouin_OSS::spinor12, &Dirac_Wilson_Brillouin_OSS::spinor13,
  &Dirac_Wilson_Brillouin_OSS::spinor14, &Dirac_Wilson_Brillouin_OSS::spinor15,
  &Dirac_Wilson_Brillouin_OSS::spinor16, &Dirac_Wilson_Brillouin_OSS::spinor17,
  &Dirac_Wilson_Brillouin_OSS::spinor18, &Dirac_Wilson_Brillouin_OSS::spinor19,
  &Dirac_Wilson_Brillouin_OSS::spinor20, &Dirac_Wilson_Brillouin_OSS::spinor21,
  &Dirac_Wilson_Brillouin_OSS::spinor22, &Dirac_Wilson_Brillouin_OSS::spinor23,
  &Dirac_Wilson_Brillouin_OSS::spinor24, &Dirac_Wilson_Brillouin_OSS::spinor25,
  &Dirac_Wilson_Brillouin_OSS::spinor26, &Dirac_Wilson_Brillouin_OSS::spinor27,
  &Dirac_Wilson_Brillouin_OSS::spinor28, &Dirac_Wilson_Brillouin_OSS::spinor29,
  &Dirac_Wilson_Brillouin_OSS::spinor30, &Dirac_Wilson_Brillouin_OSS::spinor31,
  &Dirac_Wilson_Brillouin_OSS::spinor32, &Dirac_Wilson_Brillouin_OSS::spinor33,
  &Dirac_Wilson_Brillouin_OSS::spinor34, &Dirac_Wilson_Brillouin_OSS::spinor35,
  &Dirac_Wilson_Brillouin_OSS::spinor36, &Dirac_Wilson_Brillouin_OSS::spinor37,
  &Dirac_Wilson_Brillouin_OSS::spinor38, &Dirac_Wilson_Brillouin_OSS::spinor39,
  &Dirac_Wilson_Brillouin_OSS::spinor40, &Dirac_Wilson_Brillouin_OSS::spinor41,
  &Dirac_Wilson_Brillouin_OSS::spinor42, &Dirac_Wilson_Brillouin_OSS::spinor43,
  &Dirac_Wilson_Brillouin_OSS::spinor44, &Dirac_Wilson_Brillouin_OSS::spinor45,
  &Dirac_Wilson_Brillouin_OSS::spinor46, &Dirac_Wilson_Brillouin_OSS::spinor47,
  &Dirac_Wilson_Brillouin_OSS::spinor48, &Dirac_Wilson_Brillouin_OSS::spinor49,
  &Dirac_Wilson_Brillouin_OSS::spinor50, &Dirac_Wilson_Brillouin_OSS::spinor51,
  &Dirac_Wilson_Brillouin_OSS::spinor52, &Dirac_Wilson_Brillouin_OSS::spinor53,
  &Dirac_Wilson_Brillouin_OSS::spinor54, &Dirac_Wilson_Brillouin_OSS::spinor55,
  &Dirac_Wilson_Brillouin_OSS::spinor56, &Dirac_Wilson_Brillouin_OSS::spinor57,
  &Dirac_Wilson_Brillouin_OSS::spinor58, &Dirac_Wilson_Brillouin_OSS::spinor59,
  &Dirac_Wilson_Brillouin_OSS::spinor60, &Dirac_Wilson_Brillouin_OSS::spinor61,
  &Dirac_Wilson_Brillouin_OSS::spinor62, &Dirac_Wilson_Brillouin_OSS::spinor63,
  &Dirac_Wilson_Brillouin_OSS::spinor64, &Dirac_Wilson_Brillouin_OSS::spinor65,
  &Dirac_Wilson_Brillouin_OSS::spinor66, &Dirac_Wilson_Brillouin_OSS::spinor67,
  &Dirac_Wilson_Brillouin_OSS::spinor68, &Dirac_Wilson_Brillouin_OSS::spinor69,
  &Dirac_Wilson_Brillouin_OSS::spinor70, &Dirac_Wilson_Brillouin_OSS::spinor71,
  &Dirac_Wilson_Brillouin_OSS::spinor72, &Dirac_Wilson_Brillouin_OSS::spinor73,
  &Dirac_Wilson_Brillouin_OSS::spinor74, &Dirac_Wilson_Brillouin_OSS::spinor75,
  &Dirac_Wilson_Brillouin_OSS::spinor76, &Dirac_Wilson_Brillouin_OSS::spinor77,
  &Dirac_Wilson_Brillouin_OSS::spinor78, &Dirac_Wilson_Brillouin_OSS::spinor79,
};

const GaugeField1D (Dirac_Wilson_Brillouin_OSS::*Dirac_Wilson_Brillouin_OSS::hopterms[])(const GaugeField&)const
={
  &Dirac_Wilson_Brillouin_OSS::hop_pX,&Dirac_Wilson_Brillouin_OSS::hop_mX, // 0, 1
  &Dirac_Wilson_Brillouin_OSS::hop_pY,&Dirac_Wilson_Brillouin_OSS::hop_mY, // 2, 3
  &Dirac_Wilson_Brillouin_OSS::hop_pZ,&Dirac_Wilson_Brillouin_OSS::hop_mZ, // 4, 5
  &Dirac_Wilson_Brillouin_OSS::hop_pT,&Dirac_Wilson_Brillouin_OSS::hop_mT, // 6, 7
  &Dirac_Wilson_Brillouin_OSS::hop_pXpY,&Dirac_Wilson_Brillouin_OSS::hop_mXpY, // 8, 9 
  &Dirac_Wilson_Brillouin_OSS::hop_pXmY,&Dirac_Wilson_Brillouin_OSS::hop_mXmY, // 10, 11
  &Dirac_Wilson_Brillouin_OSS::hop_pXpZ,&Dirac_Wilson_Brillouin_OSS::hop_mXpZ, // 12, 13
  &Dirac_Wilson_Brillouin_OSS::hop_pXmZ,&Dirac_Wilson_Brillouin_OSS::hop_mXmZ, // 14, 15
  &Dirac_Wilson_Brillouin_OSS::hop_pXpT,&Dirac_Wilson_Brillouin_OSS::hop_mXpT, // 16, 17
  &Dirac_Wilson_Brillouin_OSS::hop_pXmT,&Dirac_Wilson_Brillouin_OSS::hop_mXmT, // 18, 19 
  &Dirac_Wilson_Brillouin_OSS::hop_pYpZ,&Dirac_Wilson_Brillouin_OSS::hop_mYpZ, // 20, 21
  &Dirac_Wilson_Brillouin_OSS::hop_pYmZ,&Dirac_Wilson_Brillouin_OSS::hop_mYmZ, // 22, 23
  &Dirac_Wilson_Brillouin_OSS::hop_pYpT,&Dirac_Wilson_Brillouin_OSS::hop_mYpT, // 24, 25
  &Dirac_Wilson_Brillouin_OSS::hop_pYmT,&Dirac_Wilson_Brillouin_OSS::hop_mYmT, // 26, 27
  &Dirac_Wilson_Brillouin_OSS::hop_pZpT,&Dirac_Wilson_Brillouin_OSS::hop_mZpT, // 28, 29 
  &Dirac_Wilson_Brillouin_OSS::hop_pZmT,&Dirac_Wilson_Brillouin_OSS::hop_mZmT, // 30, 31
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYpZ,&Dirac_Wilson_Brillouin_OSS::hop_mXpYpZ,// 32, 33
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYpT,&Dirac_Wilson_Brillouin_OSS::hop_mXpYpT,// 34, 35
  &Dirac_Wilson_Brillouin_OSS::hop_pXpZpT,&Dirac_Wilson_Brillouin_OSS::hop_mXpZpT,// 36, 37
  &Dirac_Wilson_Brillouin_OSS::hop_pYpZpT,&Dirac_Wilson_Brillouin_OSS::hop_mYpZpT,// 38, 39
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYmZ,&Dirac_Wilson_Brillouin_OSS::hop_mXpYmZ,// 40, 41
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYmT,&Dirac_Wilson_Brillouin_OSS::hop_mXpYmT, // 42, 43
  &Dirac_Wilson_Brillouin_OSS::hop_pXpZmT,&Dirac_Wilson_Brillouin_OSS::hop_mXpZmT, // 44, 45
  &Dirac_Wilson_Brillouin_OSS::hop_pYpZmT,&Dirac_Wilson_Brillouin_OSS::hop_mYpZmT, // 46, 47
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYpZ,&Dirac_Wilson_Brillouin_OSS::hop_mXmYpZ, // 48, 49
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYpT,&Dirac_Wilson_Brillouin_OSS::hop_mXmYpT, // 50, 51
  &Dirac_Wilson_Brillouin_OSS::hop_pXmZpT,&Dirac_Wilson_Brillouin_OSS::hop_mXmZpT, // 52, 53
  &Dirac_Wilson_Brillouin_OSS::hop_pYmZpT,&Dirac_Wilson_Brillouin_OSS::hop_mYmZpT, // 54, 55
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYmZ,&Dirac_Wilson_Brillouin_OSS::hop_mXmYmZ, // 56, 57
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYmT,&Dirac_Wilson_Brillouin_OSS::hop_mXmYmT, // 58, 59
  &Dirac_Wilson_Brillouin_OSS::hop_pXmZmT,&Dirac_Wilson_Brillouin_OSS::hop_mXmZmT, // 60, 61
  &Dirac_Wilson_Brillouin_OSS::hop_pYmZmT,&Dirac_Wilson_Brillouin_OSS::hop_mYmZmT, // 62, 63
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYpZpT,&Dirac_Wilson_Brillouin_OSS::hop_mXpYpZpT, // 64, 65
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYpZmT,&Dirac_Wilson_Brillouin_OSS::hop_mXpYpZmT, // 66, 67
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYmZpT,&Dirac_Wilson_Brillouin_OSS::hop_mXpYmZpT, // 68, 69
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYpZpT,&Dirac_Wilson_Brillouin_OSS::hop_mXmYpZpT, // 70, 71
  &Dirac_Wilson_Brillouin_OSS::hop_pXpYmZmT,&Dirac_Wilson_Brillouin_OSS::hop_mXpYmZmT, // 72, 73
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYpZmT,&Dirac_Wilson_Brillouin_OSS::hop_mXmYpZmT, // 74, 75
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYmZpT,&Dirac_Wilson_Brillouin_OSS::hop_mXmYmZpT, // 76, 77
  &Dirac_Wilson_Brillouin_OSS::hop_pXmYmZmT,&Dirac_Wilson_Brillouin_OSS::hop_mXmYmZmT, // 78, 79
  };

//1-hopping terms * 8
const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMu(const GaugeField& gu,int mu) const{
  return DirSlice(gu,mu);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mMu(const GaugeField& gu,int mu) const{
  using namespace Mapping;
  GaugeField1D Umu_mN =shiftField(DirSlice(gu,mu),mu,Backward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.set(tmp.format.islice(site),(mat_dag(Umu_mN,site)).getva());
  }
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pX(const GaugeField& gu) const{
  return hop_pMu(gu,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mX(const GaugeField& gu) const{
  return hop_mMu(gu,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pY(const GaugeField& gu) const{
  return hop_pMu(gu,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mY(const GaugeField& gu) const{
  return hop_mMu(gu,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pZ(const GaugeField& gu) const{
  return hop_pMu(gu,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mZ(const GaugeField& gu) const{
  return hop_mMu(gu,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pT(const GaugeField& gu) const{
  return hop_pMu(gu,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mT(const GaugeField& gu) const{
  return hop_mMu(gu,TDIR);
}

//2-hopping terms * 24
const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMupNu(const GaugeField& gu,int mu,int nu) const{
  using namespace Mapping;
  GaugeField1D Umu = DirSlice(gu,mu);
  GaugeField1D Unu = DirSlice(gu,nu);
  GaugeField1D Umu_pN = shiftField(Umu,nu,Forward());
  GaugeField1D Unu_pM = shiftField(Unu,mu,Forward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(Umu,site)*mat(Unu_pM,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(Unu,site)*mat(Umu_pN,site)).getva());
  }
  tmp.data*=0.5;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpY(const GaugeField& gu) const{
  return hop_pMupNu(gu,XDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpZ(const GaugeField& gu) const{
  return hop_pMupNu(gu,XDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpT(const GaugeField& gu) const{
  return hop_pMupNu(gu,XDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYpZ(const GaugeField& gu) const{
  return hop_pMupNu(gu,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYpT(const GaugeField& gu) const{
  return hop_pMupNu(gu,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pZpT(const GaugeField& gu) const{
  return hop_pMupNu(gu,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mMumNu(const GaugeField& gu,int mu,int nu) const{
  using namespace Mapping;
  GaugeField1D Umu_mM = shiftField(DirSlice(gu,mu),mu,Backward());
  GaugeField1D Unu_mN = shiftField(DirSlice(gu,nu),nu,Backward());
  GaugeField1D Umu_mMmN = shiftField(Umu_mM,nu,Backward());
  GaugeField1D Unu_mMmN = shiftField(Unu_mN,mu,Backward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat_dag(Umu_mM,site)*mat_dag(Unu_mMmN,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat_dag(Unu_mN,site)*mat_dag(Umu_mMmN,site)).getva());
  }

  tmp.data*=0.5;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmY(const GaugeField& gu) const{
  return hop_mMumNu(gu,XDIR,YDIR);
}


const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmZ(const GaugeField& gu) const{
  return hop_mMumNu(gu,XDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmT(const GaugeField& gu) const{
  return hop_mMumNu(gu,XDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYmZ(const GaugeField& gu) const{
  return hop_mMumNu(gu,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYmT(const GaugeField& gu) const{
  return hop_mMumNu(gu,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mZmT(const GaugeField& gu) const{
  return hop_mMumNu(gu,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMumNu(const GaugeField& gu,int mu,int nu) const{
  using namespace Mapping; 
  GaugeField1D Umu = DirSlice(gu,mu);
  GaugeField1D Unu = DirSlice(gu,nu);
  GaugeField1D Umu_mNu = shiftField(Umu,nu,Backward());
  GaugeField1D Unu_mNu = shiftField(Unu,nu,Backward());
  GaugeField1D Unu_pMumNu = shiftField(Unu_mNu,mu,Forward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(Umu,site)*mat_dag(Unu_pMumNu,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat_dag(Unu_mNu,site)*mat(Umu_mNu,site)).getva());
  }
  tmp.data*=0.5;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmY(const GaugeField& gu) const{
  return hop_pMumNu(gu,XDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmZ(const GaugeField& gu) const{
  return hop_pMumNu(gu,XDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmT(const GaugeField& gu) const{
  return hop_pMumNu(gu,XDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYmZ(const GaugeField& gu) const{
  return hop_pMumNu(gu,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYmT(const GaugeField& gu) const{
  return hop_pMumNu(gu,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pZmT(const GaugeField& gu) const{
  return hop_pMumNu(gu,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpY(const GaugeField& gu) const{
  return hop_pMumNu(gu,YDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpZ(const GaugeField& gu) const{
  return hop_pMumNu(gu,ZDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpT(const GaugeField& gu) const{
  return hop_pMumNu(gu,TDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYpZ(const GaugeField& gu) const{
return hop_pMumNu(gu,ZDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYpT(const GaugeField& gu) const{
return hop_pMumNu(gu,TDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mZpT(const GaugeField& gu) const{
return hop_pMumNu(gu,TDIR,ZDIR);
}

//3-hopping terms * 32 
const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMupNupRho(const GaugeField& gu,int mu,int nu,int rho) const{
  using namespace Mapping;
  GaugeField1D V_pmu_pnu  = hop_pMupNu(gu,mu, nu);
  GaugeField1D V_pmu_prho = hop_pMupNu(gu,mu,rho);
  GaugeField1D V_pnu_prho = hop_pMupNu(gu,nu,rho);
  GaugeField1D Umu_pnu_prho = shiftField(shiftField(DirSlice(gu, mu), nu,Forward()),rho,Forward());
  GaugeField1D Unu_pmu_prho = shiftField(shiftField(DirSlice(gu, nu), mu,Forward()),rho,Forward());
  GaugeField1D Urho_pmu_pnu = shiftField(shiftField(DirSlice(gu,rho), mu,Forward()), nu,Forward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu,site) *mat(Urho_pmu_pnu,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_prho,site)*mat(Unu_pmu_prho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pnu_prho,site)*mat(Umu_pnu_prho,site)).getva());
  }
  tmp.data/=3.0;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYpZ(const GaugeField& gu) const{
  return hop_pMupNupRho(gu,XDIR,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYpT(const GaugeField& gu) const{
  return hop_pMupNupRho(gu,XDIR,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpZpT(const GaugeField& gu) const{
  return hop_pMupNupRho(gu,XDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYpZpT(const GaugeField& gu) const{
  return hop_pMupNupRho(gu,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMupNumRho(const GaugeField& gu,int mu,int nu,int rho) const{
  using namespace Mapping;
  GaugeField1D V_pmu_pnu  = hop_pMupNu(gu,mu, nu);
  GaugeField1D V_pmu_mrho = hop_pMumNu(gu,mu,rho);
  GaugeField1D V_pnu_mrho = hop_pMumNu(gu,nu,rho);
  GaugeField1D Umu_pnu_mrho = shiftField(shiftField(DirSlice(gu, mu),nu, Forward()),rho,Backward());
  GaugeField1D Unu_pmu_mrho = shiftField(shiftField(DirSlice(gu, nu),mu, Forward()),rho,Backward());
  GaugeField1D Urho_pmu_pnu_mrho = shiftField(shiftField(shiftField(DirSlice(gu,rho),mu, Forward()), nu, Forward()),rho,Backward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat( V_pmu_pnu,site)*mat_dag(Urho_pmu_pnu_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_mrho,site)*    mat(Unu_pmu_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pnu_mrho,site)*    mat(Umu_pnu_mrho,site)).getva());
  }
  tmp.data/=3.0;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYmZ(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,XDIR,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYmT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,XDIR,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpZmT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,XDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYpZmT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYpZ(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,XDIR,ZDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYpT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,XDIR,TDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmZpT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,XDIR,TDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYmZpT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,YDIR,TDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYpZ(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,ZDIR,YDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYpT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,TDIR,YDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpZpT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,TDIR,ZDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYpZpT(const GaugeField& gu) const{
  return hop_pMupNumRho(gu,TDIR,ZDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMumNumRho(const GaugeField& gu,int mu,int nu,int rho) const{
  using namespace Mapping;
  GaugeField1D V_pmu_mnu  = hop_pMumNu(gu,mu, nu);
  GaugeField1D V_pmu_mrho = hop_pMumNu(gu,mu,rho);
  GaugeField1D V_mnu_mrho = hop_mMumNu(gu,nu,rho);
  GaugeField1D Umu_mnu_mrho = shiftField(shiftField(DirSlice(gu,mu), nu,Backward()),rho,Backward());
  GaugeField1D Unu_pmu_mrho_mnu  = shiftField(shiftField(shiftField(DirSlice(gu,nu), mu,Forward()), rho,Backward()),nu,Backward());
  GaugeField1D Urho_pmu_mnu_mrho = shiftField(shiftField(shiftField(DirSlice(gu,rho),mu,Forward()),  nu,Backward()),rho,Backward());
  GaugeField1D tmp;

  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_mnu,site) *mat_dag(Urho_pmu_mnu_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_mrho,site)*mat_dag(Unu_pmu_mrho_mnu,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mnu_mrho,site)*    mat(Umu_mnu_mrho,site)).getva());
  }
  tmp.data/=3.0;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYmZ(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,XDIR,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYmT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,XDIR,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmZmT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,XDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pYmZmT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYmZ(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,YDIR,XDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYmT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,YDIR,XDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpZmT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,ZDIR,XDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYpZmT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,ZDIR,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYpZ(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,ZDIR,YDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYpT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,TDIR,YDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmZpT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,TDIR,ZDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYmZpT(const GaugeField& gu) const{
  return hop_pMumNumRho(gu,TDIR,ZDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mMumNumRho(const GaugeField& gu,int mu,int nu,int rho) const{
  using namespace Mapping;
  GaugeField1D V_mmu_mnu  = hop_mMumNu(gu,mu, nu);
  GaugeField1D V_mmu_mrho = hop_mMumNu(gu,mu,rho);
  GaugeField1D V_mnu_mrho = hop_mMumNu(gu,nu,rho);
  GaugeField1D Umu_mnu_mrho_mmu = shiftField(shiftField(shiftField(DirSlice(gu, mu),nu,Backward()),rho,Backward()),mu,Backward());
  GaugeField1D Unu_mmu_mrho_mnu = shiftField(shiftField(shiftField(DirSlice(gu, nu),mu,Backward()),rho,Backward()),nu,Backward());
  GaugeField1D Urho_mmu_mnu_mrho = shiftField(shiftField(shiftField(DirSlice(gu,rho),mu,Backward()),nu,Backward()),rho,Backward());
  GaugeField1D tmp;
  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site), (mat(V_mmu_mnu,site)*mat_dag(Urho_mmu_mnu_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mmu_mrho,site)*mat_dag(Unu_mmu_mrho_mnu,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mnu_mrho,site)*mat_dag(Umu_mnu_mrho_mmu,site)).getva());
  }
  tmp.data/=3.0;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYmZ(const GaugeField& gu) const{
  return hop_mMumNumRho(gu,XDIR,YDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYmT(const GaugeField& gu) const{
  return hop_mMumNumRho(gu,XDIR,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmZmT(const GaugeField& gu) const{
  return hop_mMumNumRho(gu,XDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mYmZmT(const GaugeField& gu) const{
  return hop_mMumNumRho(gu,YDIR,ZDIR,TDIR);
}


//4-hopping terms * 16
const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMupNupRhopSgm(const GaugeField& gu,int mu,int nu,int rho,int sgm) const{
  using namespace Mapping;
  GaugeField1D V_pmu_pnu_prho   = hop_pMupNupRho(gu,mu, nu,rho);
  GaugeField1D V_pmu_pnu_psgm   = hop_pMupNupRho(gu,mu, nu,sgm);
  GaugeField1D V_pmu_prho_psgm  = hop_pMupNupRho(gu,mu,rho,sgm);
  GaugeField1D V_pnu_prho_psgm  = hop_pMupNupRho(gu,nu,rho,sgm);
  GaugeField1D Umu_pnu_prho_psgm = shiftField(shiftField(shiftField(DirSlice(gu, mu),nu,Forward()),rho,Forward()),sgm,Forward());
  GaugeField1D Unu_pmu_prho_psgm = shiftField(shiftField(shiftField(DirSlice(gu, nu),rho,Forward()),mu,Forward()),sgm,Forward());
  GaugeField1D Urho_pmu_pnu_psgm = shiftField(shiftField(shiftField(DirSlice(gu,rho),mu,Forward()), nu,Forward()),sgm,Forward());
  GaugeField1D Usgm_pmu_pnu_prho = shiftField(shiftField(shiftField(DirSlice(gu,sgm),mu,Forward()), nu,Forward()),rho,Forward());
  GaugeField1D tmp;

  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu_prho, site)*mat(Usgm_pmu_pnu_prho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu_psgm, site)*mat(Urho_pmu_pnu_psgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_prho_psgm,site)*mat(Unu_pmu_prho_psgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pnu_prho_psgm,site)*mat(Umu_pnu_prho_psgm,site)).getva());
  }
  tmp.data*=0.25;

  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYpZpT(const GaugeField& gu) const{
  return hop_pMupNupRhopSgm(gu,XDIR,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMupNupRhomSgm(const GaugeField& gu,int mu,int nu,int rho,int sgm) const{
  using namespace Mapping;
  GaugeField1D V_pmu_pnu_prho   = hop_pMupNupRho(gu,mu,nu,rho);
  GaugeField1D V_pmu_pnu_msgm   = hop_pMupNumRho(gu,mu,nu,sgm);
  GaugeField1D V_pmu_prho_msgm  = hop_pMupNumRho(gu,mu,rho,sgm);
  GaugeField1D V_pnu_prho_msgm  = hop_pMupNumRho(gu,nu,rho,sgm);
  GaugeField1D Umu_pnu_prho_msgm = shiftField(shiftField(shiftField(DirSlice(gu, mu),nu,Forward()),rho,Forward()),sgm,Backward());
  GaugeField1D Unu_pmu_prho_msgm = shiftField(shiftField(shiftField(DirSlice(gu, nu),mu,Forward()),rho,Forward()),sgm,Backward());
  GaugeField1D Urho_pmu_pnu_msgm = shiftField(shiftField(shiftField(DirSlice(gu,rho),mu,Forward()), nu,Forward()),sgm,Backward());
  GaugeField1D Usgm_pmu_pnu_prho_msgm
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,sgm),mu,Forward()), nu,Forward()),rho, Forward()),sgm,Backward());
  GaugeField1D tmp;

  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu_prho,site) *mat_dag(Usgm_pmu_pnu_prho_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu_msgm,site) *mat(Urho_pmu_pnu_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_prho_msgm,site)*mat(Unu_pmu_prho_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pnu_prho_msgm,site)*mat(Umu_pnu_prho_msgm,site)).getva());
  }
  tmp.data*=0.25;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYpZmT(const GaugeField& gu) const{
  return hop_pMupNupRhomSgm(gu,XDIR,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYmZpT(const GaugeField& gu) const{
  return hop_pMupNupRhomSgm(gu,XDIR,YDIR,TDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYpZpT(const GaugeField& gu) const{
  return hop_pMupNupRhomSgm(gu,XDIR,TDIR,ZDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYpZpT(const GaugeField& gu) const{
  return hop_pMupNupRhomSgm(gu,TDIR,YDIR,ZDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMupNumRhomSgm(const GaugeField& gu,int mu,int nu,int rho,int sgm) const{
  using namespace Mapping;
  GaugeField1D V_pmu_pnu_mrho   = hop_pMupNumRho(gu,mu,nu,rho);
  GaugeField1D V_pmu_pnu_msgm   = hop_pMupNumRho(gu,mu,nu,sgm);
  GaugeField1D V_pmu_mrho_msgm  = hop_pMumNumRho(gu,mu,rho,sgm);
  GaugeField1D V_pnu_mrho_msgm  = hop_pMumNumRho(gu,nu,rho,sgm);
  GaugeField1D Umu_pnu_mrho_msgm = shiftField(shiftField(shiftField(DirSlice(gu, mu),nu,Forward()),rho,Backward()),sgm,Backward());
  GaugeField1D Unu_pmu_mrho_msgm = shiftField(shiftField(shiftField(DirSlice(gu, nu),mu,Forward()),rho,Backward()),sgm,Backward());
  GaugeField1D Urho_pmu_pnu_msgm_mrho
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,rho),mu,Forward()),nu,Forward()),sgm,Backward()),rho,Backward());
  GaugeField1D Usgm_pmu_pnu_mrho_msgm
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,sgm),mu,Forward()),nu,Forward()),rho,Backward()),sgm,Backward());
  GaugeField1D tmp;

  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu_mrho,site) *mat_dag(Usgm_pmu_pnu_mrho_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_pnu_msgm,site) *mat_dag(Urho_pmu_pnu_msgm_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_mrho_msgm,site)*mat(Unu_pmu_mrho_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pnu_mrho_msgm,site)*mat(Umu_pnu_mrho_msgm,site)).getva());
  }
  tmp.data*=0.25;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXpYmZmT(const GaugeField& gu) const{
  return hop_pMupNumRhomSgm(gu,XDIR,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYpZmT(const GaugeField& gu) const{
  return hop_pMupNumRhomSgm(gu,XDIR,ZDIR,YDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYmZpT(const GaugeField& gu) const{
  return hop_pMupNumRhomSgm(gu,XDIR,TDIR,ZDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYpZmT(const GaugeField& gu) const{
  return hop_pMupNumRhomSgm(gu,ZDIR,YDIR,XDIR,TDIR);
}
const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYmZpT(const GaugeField& gu) const{
  return hop_pMupNumRhomSgm(gu,YDIR,TDIR,XDIR,ZDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYpZpT(const GaugeField& gu) const{
  return hop_pMupNumRhomSgm(gu,ZDIR,TDIR,XDIR,YDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pMumNumRhomSgm(const GaugeField& gu,int mu,int nu,int rho,int sgm) const{
  using namespace Mapping;
  GaugeField1D V_pmu_mnu_mrho   = hop_pMumNumRho(gu,mu,nu,rho);
  GaugeField1D V_pmu_mnu_msgm   = hop_pMumNumRho(gu,mu,nu,sgm);
  GaugeField1D V_pmu_mrho_msgm  = hop_pMumNumRho(gu,mu,rho,sgm);
  GaugeField1D V_mnu_mrho_msgm  = hop_mMumNumRho(gu,nu,rho,sgm);
  GaugeField1D Umu_mnu_mrho_msgm = shiftField(shiftField(shiftField(DirSlice(gu, mu),nu,Backward()),rho,Backward()),sgm,Backward());
  GaugeField1D Unu_pmu_mrho_msgm_mnu
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu, nu),mu, Forward()),rho,Backward()),sgm,Backward()),nu, Backward());
  GaugeField1D Urho_pmu_mnu_msgm_mrho
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,rho),mu, Forward()), nu,Backward()),sgm,Backward()),rho,Backward());
  GaugeField1D Usgm_pmu_mnu_mrho_msgm
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,sgm),mu, Forward()), nu,Backward()),rho,Backward()),sgm,Backward());
  GaugeField1D tmp;

  for(int site=0;site<Nvol_;++site){
    tmp.data.set(tmp.format.islice(site),(mat(V_pmu_mnu_mrho,site) *mat_dag(Usgm_pmu_mnu_mrho_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_mnu_msgm,site) *mat_dag(Urho_pmu_mnu_msgm_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_pmu_mrho_msgm,site)*mat_dag(Unu_pmu_mrho_msgm_mnu,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mnu_mrho_msgm,site)*mat(Umu_mnu_mrho_msgm,site)).getva());
  }
  tmp.data*=0.25;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_pXmYmZmT(const GaugeField& gu) const{
  return hop_pMumNumRhomSgm(gu,XDIR,YDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXpYmZmT(const GaugeField& gu) const{
  return hop_pMumNumRhomSgm(gu,YDIR,XDIR,ZDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYpZmT(const GaugeField& gu) const{
  return hop_pMumNumRhomSgm(gu,ZDIR,YDIR,XDIR,TDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYmZpT(const GaugeField& gu) const{
  return hop_pMumNumRhomSgm(gu,TDIR,YDIR,ZDIR,XDIR);
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mMumNumRhomSgm(const GaugeField& gu,int mu,int nu,int rho,int sgm) const{
  using namespace Mapping;
  GaugeField1D V_mmu_mnu_mrho   = hop_mMumNumRho(gu,mu,nu,rho);
  GaugeField1D V_mmu_mnu_msgm   = hop_mMumNumRho(gu,mu,nu,sgm);
  GaugeField1D V_mmu_mrho_msgm  = hop_mMumNumRho(gu,mu,rho,sgm);
  GaugeField1D V_mnu_mrho_msgm  = hop_mMumNumRho(gu,nu,rho,sgm);
  GaugeField1D Umu_mnu_mrho_msgm_mmu
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu, mu),nu,Backward()),rho,Backward()),sgm,Backward()),mu,Backward());
  GaugeField1D Unu_pmu_mrho_msgm_mnu
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu, nu),mu,Backward()),rho,Backward()),sgm,Backward()),nu,Backward());
  GaugeField1D Urho_pmu_mnu_msgm_mrho
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,rho),mu,Backward()), nu,Backward()),sgm,Backward()),rho,Backward());
  GaugeField1D Usgm_pmu_mnu_mrho_msgm
    = shiftField(shiftField(shiftField(shiftField(DirSlice(gu,sgm),mu,Backward()), nu,Backward()),rho,Backward()),sgm,Backward());
  GaugeField1D tmp;

  for(int site=0;site<Nvol_;++site){
    tmp.data.add(tmp.format.islice(site),(mat(V_mmu_mnu_mrho,site) *mat_dag(Usgm_pmu_mnu_mrho_msgm,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mmu_mnu_msgm,site) *mat_dag(Urho_pmu_mnu_msgm_mrho,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mmu_mrho_msgm,site)*mat_dag(Unu_pmu_mrho_msgm_mnu,site)).getva());
    tmp.data.add(tmp.format.islice(site),(mat(V_mnu_mrho_msgm,site)*mat_dag(Umu_mnu_mrho_msgm_mmu,site)).getva());
  }
  tmp.data*=0.25;
  return tmp;
}

const GaugeField1D Dirac_Wilson_Brillouin_OSS::hop_mXmYmZmT(const GaugeField& gu) const{
  return hop_mMumNumRhomSgm(gu,XDIR,YDIR,ZDIR,TDIR);
}
