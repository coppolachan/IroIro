/* !@filename dirac_wilson_Brillouin_OSS.cpp
 * @brief implementation of Dirac_Wilson_Brillouin_OSS class
 *  Time-stamp: <2014-01-20 12:51:11 noaki>
 */
#include "dirac_wilson_Brillouin_OSS.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include "Tools/fieldUtils.hpp"
#include "Communicator/comm_io.hpp"
#include "Fields/field_expressions.hpp"
#include "Geometry/shiftField.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#endif

using namespace Mapping;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace SUNvecUtils;
using namespace std;

const Field (Dirac_Wilson_Brillouin_OSS::*Dirac_Wilson_Brillouin_OSS::der[])
()const = {&Dirac_Wilson_Brillouin_OSS::der_iso_X,
	   &Dirac_Wilson_Brillouin_OSS::der_iso_Y,
	   &Dirac_Wilson_Brillouin_OSS::der_iso_Z,
	   &Dirac_Wilson_Brillouin_OSS::der_iso_T,};

const Field Dirac_Wilson_Brillouin_OSS::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

void Dirac_Wilson_Brillouin_OSS::mult_std(Field& w,const Field& f)const{
  using namespace FieldExpression;
  CCIO::cout<<"mult(non-improved Brillouin_OSS) was called!"<<std::endl;
  Fsetup(FermionField(f));
  w = der_iso();
  w -= 0.5*lap_bri();
  w *= kbr_;
  w += f;
}

//mult for the improved Brillouin_OSS
void Dirac_Wilson_Brillouin_OSS::mult_imp(Field& w,const Field& f)const{
  using namespace FieldExpression;
  /*
  Fsetup(f);
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

int Dirac_Wilson_Brillouin_OSS::sgm(int mu,int nu,int rho)const{
  for(int sigma=0; sigma<NDIM_; ++sigma){
    if(sigma !=mu && sigma !=nu && sigma !=rho) return sigma;
  }
}

const Field Dirac_Wilson_Brillouin_OSS::gamma(const Field& f,int mu)const{
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site)
    (dm_.*(dm_.gamma[mu]))(w.getaddr(ff_.index(0,site)),
                           const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}

const Field Dirac_Wilson_Brillouin_OSS::gamma5(const Field& f)const{ 
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site)
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
		   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}

void Dirac_Wilson_Brillouin_OSS::mult_col(FermionField& f,int dir)const{
  for(int site=0;site<Nvol_;++site){
    const double* v = f.data.getaddr(ff_.index(0,site));	  
    const double* U = const_cast<Field&>(W_).getaddr(Wfmt_.index(0,site,dir));
    for(int s=0; s<ND_; ++s){
      for(int c=0; c<NC_; ++c){
	double vr = 0.0, vi = 0.0;
	for(int c1=0; c1<NC_; ++c1){
	  vr += U[re(c,c1)]*v[sr(s,c1)] -U[im(c,c1)]*v[si(s,c1)];
	  vi += U[im(c,c1)]*v[sr(s,c1)] +U[re(c,c1)]*v[si(s,c1)];
	}
	F_[dir].add(ff_.index_r(c,s,site),vr);
	F_[dir].add(ff_.index_i(c,s,site),vi);
      }
    }
  }
}

void Dirac_Wilson_Brillouin_OSS::Wsetup(const GaugeField& gu){ 
  GaugeField1D g;
  for(int dir=0; dir<80; ++dir){
    (this->*ghops[dir])(g,gu);
    W_.set(Wfmt_.ex_slice(dir),g.getva());
  }
}

const Field Dirac_Wilson_Brillouin_OSS::der_iso()const{
  Field h(fsize_);
  for(int mu=0; mu<NDIM_; ++mu) h += gamma((this->*der[mu])(),mu);
  return h;
}

//const Field Dirac_Wilson_Brillouin_OSS::lap_bri(const Field& f)const{
const Field Dirac_Wilson_Brillouin_OSS::lap_bri()const{
  using namespace FieldExpression;
  Field h(fsize_),tmp(fsize_);
  double c[NDIM_+1] = {-240.0/64.0, 8.0/64.0, 4.0/64.0, 2.0/64.0 ,1.0/64.0};
  // 1-hopping terms
  for(int i=0; i<8; ++i)  h += c[1]*F_[i];
  // 2-hopping terms  
  for(int i=8; i<32; ++i) h += c[2]*F_[i];
  // 3-hopping terms  
  for(int i=32; i<64; ++i) h += c[3]*F_[i];
  // 4-hopping terms  
  for(int i=64; i<80; ++i) h += c[4]*F_[i];
  return h;
}

const Field Dirac_Wilson_Brillouin_OSS::der_iso_X()const{
  using namespace FieldExpression;
  Field h(fsize_);
  double c[NDIM_+1] = {0.0, 64.0/432.0, 16.0/432.0, 4.0/432.0, 1.0/432.0};
  
  // 1-hopping terms
  h += c[1]*F_[0];  h -= c[1]*F_[1];

  // 2-hopping terms	 	
  for(int i=8;i<=30;i+=2){
    if(i!=20 && i!=22 && i!=24 && i!=26 && i!=28 && i!=30){
      h += c[2]*F_[i]; h -= c[2]*F_[i+1];
    }
  }
  // 3-hopping terms
  for(int i=32;i<=62;i+=2){
    if(i!=38 && i!=46 && i!=54 && i!=62){
      h += c[3]*F_[i];  h-= c[3]*F_[i+1];
    }	
  }
  // 4-hopping terms
  for(int i=64;i<=78;i+=2){
    h += c[4]*F_[i];  h-= c[4]*F_[i+1];
  }
  return h;
}

const Field Dirac_Wilson_Brillouin_OSS::der_iso_Y()const{
  using namespace FieldExpression;
  using namespace Mapping;
  Field h(fsize_);
  double c[NDIM_+1] = {0.0, 64.0/432.0, 16.0/432.0, 4.0/432.0, 1.0/432.0};
  
  // 1-hopping terms
  h += c[1]*F_[2];  h -= c[1]*F_[3];

  // 2-hopping terms	 	
  h += c[2]*F_[  8]; h -= c[2]*F_[ 10];
  h += c[2]*F_[  9]; h -= c[2]*F_[ 11];
  h += c[2]*F_[ 20]; h -= c[2]*F_[ 21];
  h += c[2]*F_[ 22]; h -= c[2]*F_[ 23];
  h += c[2]*F_[ 24]; h -= c[2]*F_[ 25];
  h += c[2]*F_[ 26]; h -= c[2]*F_[ 27];

  // 3-hopping terms
  h += c[3]*F_[32];  h-= c[3]*F_[48];
  h += c[3]*F_[33];  h-= c[3]*F_[49];
  h += c[3]*F_[34];  h-= c[3]*F_[50];
  h += c[3]*F_[35];  h-= c[3]*F_[51];
  h += c[3]*F_[38];  h-= c[3]*F_[39];
  h += c[3]*F_[40];  h-= c[3]*F_[56];
  h += c[3]*F_[41];  h-= c[3]*F_[57];
  h += c[3]*F_[42];  h-= c[3]*F_[58];
  h += c[3]*F_[43];  h-= c[3]*F_[59];
  h += c[3]*F_[46];  h-= c[3]*F_[47];
  h += c[3]*F_[54];  h-= c[3]*F_[55];
  h += c[3]*F_[62];  h-= c[3]*F_[63];
  
  // 4-hopping terms
  h += c[4]*F_[64];  h-= c[4]*F_[70];
  h += c[4]*F_[65];  h-= c[4]*F_[71];
  h += c[4]*F_[66];  h-= c[4]*F_[74];
  h += c[4]*F_[67];  h-= c[4]*F_[75];
  h += c[4]*F_[68];  h-= c[4]*F_[76];
  h += c[4]*F_[69];  h-= c[4]*F_[77];
  h += c[4]*F_[72];  h-= c[4]*F_[78];
  h += c[4]*F_[73];  h-= c[4]*F_[79];

  return h;
}

const Field Dirac_Wilson_Brillouin_OSS::der_iso_Z()const{
  using namespace FieldExpression;
  using namespace Mapping;
  Field h(fsize_);
  double c[NDIM_+1] = {0.0, 64.0/432.0, 16.0/432.0, 4.0/432.0, 1.0/432.0};
  
  // 1-hopping terms
  h += c[1]*F_[4];  h -= c[1]*F_[5];

  // 2-hopping terms	 	
  h += c[2]*F_[12]; h -= c[2]*F_[14];
  h += c[2]*F_[13]; h -= c[2]*F_[15];
  h += c[2]*F_[20]; h -= c[2]*F_[22];
  h += c[2]*F_[21]; h -= c[2]*F_[23];
  h += c[2]*F_[28]; h -= c[2]*F_[29];
  h += c[2]*F_[30]; h -= c[2]*F_[31];

  // 3-hopping terms
  h += c[3]*F_[32];  h-= c[3]*F_[40];
  h += c[3]*F_[33];  h-= c[3]*F_[41];
  h += c[3]*F_[36];  h-= c[3]*F_[52];
  h += c[3]*F_[37];  h-= c[3]*F_[53];
  h += c[3]*F_[38];  h-= c[3]*F_[54];
  h += c[3]*F_[39];  h-= c[3]*F_[55];
  h += c[3]*F_[44];  h-= c[3]*F_[60];
  h += c[3]*F_[45];  h-= c[3]*F_[61];
  h += c[3]*F_[46];  h-= c[3]*F_[62];
  h += c[3]*F_[47];  h-= c[3]*F_[63];
  h += c[3]*F_[48];  h-= c[3]*F_[56];
  h += c[3]*F_[49];  h-= c[3]*F_[57];
  
  // 4-hopping terms
  h += c[4]*F_[64];  h-= c[4]*F_[68];
  h += c[4]*F_[65];  h-= c[4]*F_[69];
  h += c[4]*F_[66];  h-= c[4]*F_[72];
  h += c[4]*F_[67];  h-= c[4]*F_[73];
  h += c[4]*F_[70];  h-= c[4]*F_[76];
  h += c[4]*F_[71];  h-= c[4]*F_[77];
  h += c[4]*F_[74];  h-= c[4]*F_[78];
  h += c[4]*F_[75];  h-= c[4]*F_[79];

  return h;
}

const Field Dirac_Wilson_Brillouin_OSS::der_iso_T()const{
  using namespace FieldExpression;
  using namespace Mapping;
  Field h(fsize_);
  double c[NDIM_+1] = {0.0, 64.0/432.0, 16.0/432.0, 4.0/432.0, 1.0/432.0};
  
  // 1-hopping terms
  h += c[1]*F_[6];  h -= c[1]*F_[7];

  // 2-hopping terms	 	
  h += c[2]*F_[16]; h -= c[2]*F_[18];
  h += c[2]*F_[17]; h -= c[2]*F_[19];
  h += c[2]*F_[24]; h -= c[2]*F_[26];
  h += c[2]*F_[25]; h -= c[2]*F_[27];
  h += c[2]*F_[28]; h -= c[2]*F_[30];
  h += c[2]*F_[29]; h -= c[2]*F_[31];

  // 3-hopping terms
  h += c[3]*F_[34];  h-= c[3]*F_[42];
  h += c[3]*F_[35];  h-= c[3]*F_[43];
  h += c[3]*F_[36];  h-= c[3]*F_[44];
  h += c[3]*F_[37];  h-= c[3]*F_[45];
  h += c[3]*F_[38];  h-= c[3]*F_[46];
  h += c[3]*F_[39];  h-= c[3]*F_[47];
  h += c[3]*F_[50];  h-= c[3]*F_[58];
  h += c[3]*F_[51];  h-= c[3]*F_[59];
  h += c[3]*F_[52];  h-= c[3]*F_[60];
  h += c[3]*F_[53];  h-= c[3]*F_[61];
  h += c[3]*F_[54];  h-= c[3]*F_[62];
  h += c[3]*F_[55];  h-= c[3]*F_[63];
  
  // 4-hopping terms
  h += c[4]*F_[64];  h-= c[4]*F_[66];
  h += c[4]*F_[65];  h-= c[4]*F_[67];
  h += c[4]*F_[68];  h-= c[4]*F_[72];
  h += c[4]*F_[69];  h-= c[4]*F_[73];
  h += c[4]*F_[70];  h-= c[4]*F_[74];
  h += c[4]*F_[71];  h-= c[4]*F_[75];
  h += c[4]*F_[76];  h-= c[4]*F_[78];
  h += c[4]*F_[77];  h-= c[4]*F_[79];

  return h;
}

void (Dirac_Wilson_Brillouin_OSS::*Dirac_Wilson_Brillouin_OSS::ghops[])
(GaugeField1D&,const GaugeField&)const ={
  &Dirac_Wilson_Brillouin_OSS::gpX,&Dirac_Wilson_Brillouin_OSS::gmX, // 0, 1
  &Dirac_Wilson_Brillouin_OSS::gpY,&Dirac_Wilson_Brillouin_OSS::gmY, // 2, 3
  &Dirac_Wilson_Brillouin_OSS::gpZ,&Dirac_Wilson_Brillouin_OSS::gmZ, // 4, 5
  &Dirac_Wilson_Brillouin_OSS::gpT,&Dirac_Wilson_Brillouin_OSS::gmT, // 6, 7
  &Dirac_Wilson_Brillouin_OSS::gpXpY,&Dirac_Wilson_Brillouin_OSS::gmXpY, // 8, 9 
  &Dirac_Wilson_Brillouin_OSS::gpXmY,&Dirac_Wilson_Brillouin_OSS::gmXmY, // 10, 11
  &Dirac_Wilson_Brillouin_OSS::gpXpZ,&Dirac_Wilson_Brillouin_OSS::gmXpZ, // 12, 13
  &Dirac_Wilson_Brillouin_OSS::gpXmZ,&Dirac_Wilson_Brillouin_OSS::gmXmZ, // 14, 15
  &Dirac_Wilson_Brillouin_OSS::gpXpT,&Dirac_Wilson_Brillouin_OSS::gmXpT, // 16, 17
  &Dirac_Wilson_Brillouin_OSS::gpXmT,&Dirac_Wilson_Brillouin_OSS::gmXmT, // 18, 19 
  &Dirac_Wilson_Brillouin_OSS::gpYpZ,&Dirac_Wilson_Brillouin_OSS::gmYpZ, // 20, 21
  &Dirac_Wilson_Brillouin_OSS::gpYmZ,&Dirac_Wilson_Brillouin_OSS::gmYmZ, // 22, 23
  &Dirac_Wilson_Brillouin_OSS::gpYpT,&Dirac_Wilson_Brillouin_OSS::gmYpT, // 24, 25
  &Dirac_Wilson_Brillouin_OSS::gpYmT,&Dirac_Wilson_Brillouin_OSS::gmYmT, // 26, 27
  &Dirac_Wilson_Brillouin_OSS::gpZpT,&Dirac_Wilson_Brillouin_OSS::gmZpT, // 28, 29 
  &Dirac_Wilson_Brillouin_OSS::gpZmT,&Dirac_Wilson_Brillouin_OSS::gmZmT, // 30, 31
  &Dirac_Wilson_Brillouin_OSS::gpXpYpZ,&Dirac_Wilson_Brillouin_OSS::gmXpYpZ,// 32, 33
  &Dirac_Wilson_Brillouin_OSS::gpXpYpT,&Dirac_Wilson_Brillouin_OSS::gmXpYpT,// 34, 35
  &Dirac_Wilson_Brillouin_OSS::gpXpZpT,&Dirac_Wilson_Brillouin_OSS::gmXpZpT,// 36, 37
  &Dirac_Wilson_Brillouin_OSS::gpYpZpT,&Dirac_Wilson_Brillouin_OSS::gmYpZpT,// 38, 39
  &Dirac_Wilson_Brillouin_OSS::gpXpYmZ,&Dirac_Wilson_Brillouin_OSS::gmXpYmZ,// 40, 41
  &Dirac_Wilson_Brillouin_OSS::gpXpYmT,&Dirac_Wilson_Brillouin_OSS::gmXpYmT, // 42, 43
  &Dirac_Wilson_Brillouin_OSS::gpXpZmT,&Dirac_Wilson_Brillouin_OSS::gmXpZmT, // 44, 45
  &Dirac_Wilson_Brillouin_OSS::gpYpZmT,&Dirac_Wilson_Brillouin_OSS::gmYpZmT, // 46, 47
  &Dirac_Wilson_Brillouin_OSS::gpXmYpZ,&Dirac_Wilson_Brillouin_OSS::gmXmYpZ, // 48, 49
  &Dirac_Wilson_Brillouin_OSS::gpXmYpT,&Dirac_Wilson_Brillouin_OSS::gmXmYpT, // 50, 51
  &Dirac_Wilson_Brillouin_OSS::gpXmZpT,&Dirac_Wilson_Brillouin_OSS::gmXmZpT, // 52, 53
  &Dirac_Wilson_Brillouin_OSS::gpYmZpT,&Dirac_Wilson_Brillouin_OSS::gmYmZpT, // 54, 55
  &Dirac_Wilson_Brillouin_OSS::gpXmYmZ,&Dirac_Wilson_Brillouin_OSS::gmXmYmZ, // 56, 57
  &Dirac_Wilson_Brillouin_OSS::gpXmYmT,&Dirac_Wilson_Brillouin_OSS::gmXmYmT, // 58, 59
  &Dirac_Wilson_Brillouin_OSS::gpXmZmT,&Dirac_Wilson_Brillouin_OSS::gmXmZmT, // 60, 61
  &Dirac_Wilson_Brillouin_OSS::gpYmZmT,&Dirac_Wilson_Brillouin_OSS::gmYmZmT, // 62, 63
  &Dirac_Wilson_Brillouin_OSS::gpXpYpZpT,&Dirac_Wilson_Brillouin_OSS::gmXpYpZpT, // 64, 65
  &Dirac_Wilson_Brillouin_OSS::gpXpYpZmT,&Dirac_Wilson_Brillouin_OSS::gmXpYpZmT, // 66, 67
  &Dirac_Wilson_Brillouin_OSS::gpXpYmZpT,&Dirac_Wilson_Brillouin_OSS::gmXpYmZpT, // 68, 69
  &Dirac_Wilson_Brillouin_OSS::gpXmYpZpT,&Dirac_Wilson_Brillouin_OSS::gmXmYpZpT, // 70, 71
  &Dirac_Wilson_Brillouin_OSS::gpXpYmZmT,&Dirac_Wilson_Brillouin_OSS::gmXpYmZmT, // 72, 73
  &Dirac_Wilson_Brillouin_OSS::gpXmYpZmT,&Dirac_Wilson_Brillouin_OSS::gmXmYpZmT, // 74, 75
  &Dirac_Wilson_Brillouin_OSS::gpXmYmZpT,&Dirac_Wilson_Brillouin_OSS::gmXmYmZpT, // 76, 77
  &Dirac_Wilson_Brillouin_OSS::gpXmYmZmT,&Dirac_Wilson_Brillouin_OSS::gmXmYmZmT, // 78, 79
};

//1-hopping terms * 8
void Dirac_Wilson_Brillouin_OSS::
gpX(GaugeField1D& g,const GaugeField& gu) const{ DirSlice(g,gu,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpY(GaugeField1D& g,const GaugeField& gu) const{ DirSlice(g,gu,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpZ(GaugeField1D& g,const GaugeField& gu) const{ DirSlice(g,gu,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpT(GaugeField1D& g,const GaugeField& gu) const{ DirSlice(g,gu,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmX(GaugeField1D& g,const GaugeField& gu) const{ AntiDirSlice(g,gu,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmY(GaugeField1D& g,const GaugeField& gu) const{ AntiDirSlice(g,gu,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmZ(GaugeField1D& g,const GaugeField& gu) const{ AntiDirSlice(g,gu,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmT(GaugeField1D& g,const GaugeField& gu) const{ AntiDirSlice(g,gu,TDIR);}

//2-hopping terms * 24
/// (+,+)
void Dirac_Wilson_Brillouin_OSS::
gpXpY(GaugeField1D& g,const GaugeField& gu)const{diagPPslice(g,gu,XDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXpZ(GaugeField1D& g,const GaugeField& gu)const{diagPPslice(g,gu,XDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXpT(GaugeField1D& g,const GaugeField& gu)const{diagPPslice(g,gu,XDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYpZ(GaugeField1D& g,const GaugeField& gu)const{diagPPslice(g,gu,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYpT(GaugeField1D& g,const GaugeField& gu)const{diagPPslice(g,gu,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPslice(g,gu,ZDIR,TDIR);}
/// (-,-)
void Dirac_Wilson_Brillouin_OSS::
gmXmY(GaugeField1D& g,const GaugeField& gu)const{diagMMslice(g,gu,XDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXmZ(GaugeField1D& g,const GaugeField& gu)const{diagMMslice(g,gu,XDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXmT(GaugeField1D& g,const GaugeField& gu)const{diagMMslice(g,gu,XDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYmZ(GaugeField1D& g,const GaugeField& gu)const{diagMMslice(g,gu,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYmT(GaugeField1D& g,const GaugeField& gu)const{diagMMslice(g,gu,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmZmT(GaugeField1D& g,const GaugeField& gu)const{diagMMslice(g,gu,ZDIR,TDIR);}
/// (+,-)
void Dirac_Wilson_Brillouin_OSS::
gpXmY(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,XDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXmZ(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,XDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXmT(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,XDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYmZ(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYmT(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,ZDIR,TDIR);}
/// (-,+)
void Dirac_Wilson_Brillouin_OSS::
gmXpY(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,YDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXpZ(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,ZDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXpT(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,TDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYpZ(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,ZDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYpT(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,TDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPMslice(g,gu,TDIR,ZDIR);}

//=====================
//3-hopping terms * 32 
/// (+,+,+)
void Dirac_Wilson_Brillouin_OSS::
gpXpYpZ(GaugeField1D& g,const GaugeField& gu)const{diagPPPslice(g,gu,XDIR,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXpYpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPslice(g,gu,XDIR,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPslice(g,gu,XDIR,ZDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPslice(g,gu,YDIR,ZDIR,TDIR);}
/// (+,+,-)
void Dirac_Wilson_Brillouin_OSS::
gpXpYmZ(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,XDIR,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXpYmT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,XDIR,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,XDIR,ZDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,YDIR,ZDIR,TDIR);}
/// (+,-,+)
void Dirac_Wilson_Brillouin_OSS::
gpXmYpZ(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,XDIR,ZDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXmYpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,XDIR,TDIR,YDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,XDIR,TDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,YDIR,TDIR,ZDIR);}
/// (-,+,+)
void Dirac_Wilson_Brillouin_OSS::
gmXpYpZ(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,ZDIR,YDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXpYpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,TDIR,YDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,TDIR,ZDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMslice(g,gu,TDIR,ZDIR,YDIR);}
/// (+,-,-)
void Dirac_Wilson_Brillouin_OSS::
gpXmYmZ(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,XDIR,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXmYmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,XDIR,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpXmZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,XDIR,ZDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gpYmZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,YDIR,ZDIR,TDIR);}
/// (-,+,-)
void Dirac_Wilson_Brillouin_OSS::
gmXpYmZ(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,YDIR,XDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXpYmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,YDIR,XDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,ZDIR,XDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,ZDIR,YDIR,TDIR);}
/// (-,-,+)
void Dirac_Wilson_Brillouin_OSS::
gmXmYpZ(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,ZDIR,YDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXmYpT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,TDIR,YDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,TDIR,ZDIR,XDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPMMslice(g,gu,TDIR,ZDIR,YDIR);}
/// (-,-,-)
void Dirac_Wilson_Brillouin_OSS::
gmXmYmZ(GaugeField1D& g,const GaugeField& gu)const{diagMMMslice(g,gu,XDIR,YDIR,ZDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXmYmT(GaugeField1D& g,const GaugeField& gu)const{diagMMMslice(g,gu,XDIR,YDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmXmZmT(GaugeField1D& g,const GaugeField& gu)const{diagMMMslice(g,gu,XDIR,ZDIR,TDIR);}

void Dirac_Wilson_Brillouin_OSS::
gmYmZmT(GaugeField1D& g,const GaugeField& gu)const{diagMMMslice(g,gu,YDIR,ZDIR,TDIR);}

//4-hopping terms * 16
/// (+,+,+,+)
void Dirac_Wilson_Brillouin_OSS::
gpXpYpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPPslice(g,gu,XDIR,YDIR,ZDIR,TDIR);}
/// (+,+,+,-)
void Dirac_Wilson_Brillouin_OSS::
gpXpYpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPPPMslice(g,gu,XDIR,YDIR,ZDIR,TDIR);}
/// (+,+,-,+)
void Dirac_Wilson_Brillouin_OSS::
gpXpYmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPMslice(g,gu,XDIR,YDIR,TDIR,ZDIR);}
/// (+,-,+,+)
void Dirac_Wilson_Brillouin_OSS::
gpXmYpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPMslice(g,gu,XDIR,ZDIR,TDIR,YDIR);}
/// (-,+,+,+)
void Dirac_Wilson_Brillouin_OSS::
gmXpYpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPPMslice(g,gu,TDIR,YDIR,ZDIR,XDIR);}
/// (+,+,-,-)
void Dirac_Wilson_Brillouin_OSS::
gpXpYmZmT(GaugeField1D& g,const GaugeField& gu)const{diagPPMMslice(g,gu,XDIR,YDIR,ZDIR,TDIR);}
/// (+,-,+,-)
void Dirac_Wilson_Brillouin_OSS::
gpXmYpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPPMMslice(g,gu,XDIR,ZDIR,YDIR,TDIR);}
/// (+,-,-,+)
void Dirac_Wilson_Brillouin_OSS::
gpXmYmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMMslice(g,gu,XDIR,TDIR,YDIR,ZDIR);}
/// (-,+,+,-)
void Dirac_Wilson_Brillouin_OSS::
gmXpYpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPPMMslice(g,gu,YDIR,ZDIR,XDIR,TDIR);}
/// (-,+,-,+)
void Dirac_Wilson_Brillouin_OSS::
gmXpYmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMMslice(g,gu,YDIR,TDIR,XDIR,ZDIR);}
/// (-,-,+,+)
void Dirac_Wilson_Brillouin_OSS::
gmXmYpZpT(GaugeField1D& g,const GaugeField& gu)const{diagPPMMslice(g,gu,ZDIR,TDIR,XDIR,YDIR);}
/// (+,-,-,-)
void Dirac_Wilson_Brillouin_OSS::
gpXmYmZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMMslice(g,gu,XDIR,YDIR,ZDIR,TDIR);}
/// (-,+,-,-)
void Dirac_Wilson_Brillouin_OSS::
gmXpYmZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMMslice(g,gu,YDIR,XDIR,ZDIR,TDIR);}
/// (-,-,+,-)
void Dirac_Wilson_Brillouin_OSS::
gmXmYpZmT(GaugeField1D& g,const GaugeField& gu)const{diagPMMMslice(g,gu,ZDIR,YDIR,XDIR,TDIR);}
/// (-,-,-,+)
void Dirac_Wilson_Brillouin_OSS::
gmXmYmZpT(GaugeField1D& g,const GaugeField& gu)const{diagPMMMslice(g,gu,TDIR,YDIR,ZDIR,XDIR);}
/// (-,-,-,-)
void Dirac_Wilson_Brillouin_OSS::
gmXmYmZmT(GaugeField1D& g,const GaugeField& gu)const{diagMMMMslice(g,gu,XDIR,YDIR,ZDIR,TDIR);}


#ifdef IBM_BGQ_WILSON
void Dirac_Wilson_Brillouin_OSS::shift_fwd(FermionField& w,const FermionField& f,int dir) const{
  double* pf = const_cast<Field&>(f.data).getaddr(0);
  double* pw = const_cast<Field&>(w.data).getaddr(0);
  BGWilson_Mult_Shift_Dir(pw,gptr_,pf,dir,BGWILSON_FORWARD);
}

void Dirac_Wilson_Brillouin_OSS::shift_bwd(FermionField& w,const FermionField& f,int dir) const{
  double* pf = const_cast<Field&>(f.data).getaddr(0);
  double* pw = const_cast<Field&>(w.data).getaddr(0);
  BGWilson_Mult_Shift_Dir(pw,gptr_,pf,dir,BGWILSON_BACKWARD);
}
#else
void Dirac_Wilson_Brillouin_OSS::shift_fwd(FermionField& w,const FermionField& f,int dir) const{
  Mapping::shiftField(w,f,dir,Mapping::Forward());
}

void Dirac_Wilson_Brillouin_OSS::shift_bwd(FermionField& w,const FermionField& f,int dir) const{
  Mapping::shiftField(w,f,dir,Mapping::Backward());
}
#endif

void Dirac_Wilson_Brillouin_OSS::Fsetup(const FermionField& f)const{ 
  {// (+,0,0,0)
    FermionField fpX(fsize_);
    shift_fwd(fpX,f,XDIR);  mult_col(fpX,0); 
    {// (+,+,0,0)
      FermionField fpXpY(fsize_);
      shift_fwd(fpXpY,fpX,YDIR);  mult_col(fpXpY,8); 
      {// (+,+,+,0)
	FermionField fpXpYpZ(fsize_);
	shift_fwd(fpXpYpZ,fpXpY,ZDIR);  mult_col(fpXpYpZ,32);
	{// (+,+,+,+)
	  FermionField fpXpYpZpT(fsize_);
	  shift_fwd(fpXpYpZpT,fpXpYpZ,TDIR);  mult_col(fpXpYpZpT,64);
	}
	{// (+,+,+,-)
	  FermionField fpXpYpZmT(fsize_);
	  shift_bwd(fpXpYpZmT,fpXpYpZ,TDIR);  mult_col(fpXpYpZmT,66);  
	}
      }
      {// (+,+,0,+)
	FermionField fpXpYpT(fsize_);
	shift_fwd(fpXpYpT,fpXpY,TDIR);	mult_col(fpXpYpT,34);
      }
      {// (+,+,-,0)
	FermionField fpXpYmZ(fsize_);
	shift_bwd(fpXpYmZ,fpXpY,ZDIR);	mult_col(fpXpYmZ,40);
	{// (+,+,-,+)
	  FermionField fpXpYmZpT(fsize_);
	  shift_fwd(fpXpYmZpT,fpXpYmZ,TDIR);  mult_col(fpXpYmZpT,68); 
	}
	{// (+,+,-,-)
	  FermionField fpXpYmZmT(fsize_);
	  shift_bwd(fpXpYmZmT,fpXpYmZ,TDIR);  mult_col(fpXpYmZmT,72);
	}
      }
      {// (+,+,0,-)
	FermionField fpXpYmT(fsize_);
	shift_bwd(fpXpYmT,fpXpY,TDIR);	mult_col(fpXpYmT,42);
      }
    }
    {// (+,-,0,0)
      FermionField fpXmY(fsize_);
      shift_bwd(fpXmY,fpX,YDIR);  mult_col(fpXmY,10);
      {// (+,-,+,0)
	FermionField fpXmYpZ(fsize_);
	shift_fwd(fpXmYpZ,fpXmY,ZDIR);	mult_col(fpXmYpZ,48); 
	{// (+,-,+,+)
	  FermionField fpXmYpZpT(fsize_);
	  shift_fwd(fpXmYpZpT,fpXmYpZ,TDIR);  mult_col(fpXmYpZpT,70); 
	}
	{// (+,-,+,-)
	  FermionField fpXmYpZmT(fsize_);
	  shift_bwd(fpXmYpZmT,fpXmYpZ,TDIR);  mult_col(fpXmYpZmT,73);
	}
      }
      {// (+,-,0,+)
	FermionField fpXmYpT(fsize_);
	shift_fwd(fpXmYpT,fpXmY,TDIR); 	mult_col(fpXmYpT,50);
      }
      {// (+,-,-,0)
	FermionField fpXmYmZ(fsize_);
	shift_bwd(fpXmYmZ,fpXmY,ZDIR); 	mult_col(fpXmYmZ,56); 
      }
      {// (+,-,0,-)
	FermionField fpXmYmT(fsize_);
	shift_bwd(fpXmYmT,fpXmY,TDIR); 	mult_col(fpXmYmT,58); 
      }
    }
    {// (+,0,+,0) 
      FermionField fpXpZ(fsize_);
      shift_fwd(fpXpZ,fpX,ZDIR);  mult_col(fpXpZ,12);
      {// (+,0,+,+) 
	FermionField fpXpZpT(fsize_);
	shift_fwd(fpXpZpT,fpXpZ,TDIR); 	mult_col(fpXpZpT,36);
      }
      {// (+,0,+,-) 
	FermionField fpXpZmT(fsize_);
	shift_bwd(fpXpZmT,fpXpZ,TDIR); 	mult_col(fpXpZmT,44); 
      }
    }
    { // (+,0,-,0)
      FermionField fpXmZ(fsize_);
      shift_bwd(fpXmZ,fpX,ZDIR);  mult_col(fpXmZ,14);
      {// (+,0,-,+)
	FermionField fpXmZpT(fsize_);
	shift_fwd(fpXmZpT,fpXmZ,TDIR); 	mult_col(fpXmZpT,52);
      }
      {// (+,0,-,-)
	FermionField fpXmZmT(fsize_);
	shift_bwd(fpXmZmT,fpXmZ,TDIR); 	mult_col(fpXmZmT,60);
      }
    }
    { // (+,0,0,+)
      FermionField fpXpT(fsize_);
      shift_fwd(fpXpT,fpX,TDIR);  mult_col(fpXpT,16);
    }
    { // (+,0,0,-)
      FermionField fpXmT(fsize_);
      shift_bwd(fpXmT,fpX,TDIR);  mult_col(fpXmT,18);
    }
  }
  {// (-,0,0,0) 
    FermionField fmX(fsize_);    
    shift_bwd(fmX,f,XDIR);  mult_col(fmX,1);

    {// (-,+,0,0) 
      FermionField fmXpY(fsize_);
      shift_fwd(fmXpY,fmX,YDIR);  mult_col(fmXpY,9);
      {// (-,+,0,+) 
	FermionField fmXpYpT(fsize_);
	shift_fwd(fmXpYpT,fmXpY,TDIR);	mult_col(fmXpYpT,35);
      }
      {// (-,+,+,0) 
	FermionField fmXpYpZ(fsize_);
	shift_fwd(fmXpYpZ,fmXpY,ZDIR);	mult_col(fmXpYpZ,33);
	{// (-,+,+,+)
	  FermionField fmXpYpZpT(fsize_);
	  shift_fwd(fmXpYpZpT,fmXpYpZ,TDIR);  mult_col(fmXpYpZpT,65);
	}
	{// (-,+,+,-)
	  FermionField fmXpYpZmT(fsize_);
	  shift_bwd(fmXpYpZmT,fmXpYpZ,TDIR);  mult_col(fmXpYpZmT,67);
	}
      }
      {// (-,+,-,0) 
	FermionField fmXpYmZ(fsize_);	
	shift_bwd(fmXpYmZ,fmXpY,ZDIR);	mult_col(fmXpYmZ,41);
	{// (-,+,-,+)
	  FermionField fmXpYmZpT(fsize_);	
	  shift_fwd(fmXpYmZpT,fmXpYmZ,TDIR);  mult_col(fmXpYmZpT,69);
	}
	{// (-,+,-,-)
	  FermionField fmXpYmZmT(fsize_);	
	  shift_bwd(fmXpYmZmT,fmXpYmZ,TDIR);  mult_col(fmXpYmZmT,73);
	}
      }
      {// (-,+,0,-)
	FermionField fmXpYmT(fsize_);	
	shift_bwd(fmXpYmT,fmXpY,TDIR);	mult_col(fmXpYmT,43);
      }
    }
    {// (-,-,0,0)
      FermionField fmXmY(fsize_);
      shift_bwd(fmXmY,fmX,YDIR);  mult_col(fmXmY,11);
      {// (-,-,+,0)
	FermionField fmXmYpZ(fsize_);
	shift_fwd(fmXmYpZ,fmXmY,ZDIR);	mult_col(fmXmYpZ,49);
	{// (-,-,+,+)
	  FermionField fmXmYpZpT(fsize_);
	  shift_fwd(fmXmYpZpT,fmXmYpZ,TDIR);  mult_col(fmXmYpZpT,71);
	}
	{// (-,-,+,-)
	  FermionField fmXmYpZmT(fsize_);
	  shift_bwd(fmXmYpZmT,fmXmYpZ,TDIR);  mult_col(fmXmYpZmT,75);
	}
      }
      {// (-,-,-,0)
	FermionField fmXmYmZ(fsize_);
	shift_bwd(fmXmYmZ,fmXmY,ZDIR);	mult_col(fmXmYmZ,57);
	{// (-,-,-,+)
	  FermionField fmXmYmZpT(fsize_);
	  shift_fwd(fmXmYmZpT,fmXmYmZ,TDIR);  mult_col(fmXmYmZpT,77);
	}
	{// (-,-,-,-)
	  FermionField fmXmYmZmT(fsize_);
	  shift_bwd(fmXmYmZmT,fmXmYmZ,TDIR);  mult_col(fmXmYmZmT,79);
	}
      }
      {// (-,-,0,+)
	FermionField fmXmYpT(fsize_);
	shift_fwd(fmXmYpT,fmXmY,TDIR);	mult_col(fmXmYpT,51);
      }
      {// (-,-,0,-)
	FermionField fmXmYmT(fsize_);
	shift_bwd(fmXmYmT,fmXmY,TDIR);	mult_col(fmXmYmT,59);
      }
    }
    {// (-,0,+,0)
      FermionField fmXpZ(fsize_);
      shift_fwd(fmXpZ,fmX,ZDIR);  mult_col(fmXpZ,13);

      {// (-,0,+,+)
	FermionField fmXpZpT(fsize_);
	shift_fwd(fmXpZpT,fmXpZ,TDIR); mult_col(fmXpZpT,37);
      }
      {// (-,0,+,-)
	FermionField fmXpZmT(fsize_);
	shift_bwd(fmXpZmT,fmXpZ,TDIR);	mult_col(fmXpZmT,45);
      }
    }
    {// (-,0,-,0)
      FermionField fmXmZ(fsize_);
      shift_bwd(fmXmZ,fmX,ZDIR);  mult_col(fmXmZ,15);
      {// (-,0,-,+)
	FermionField fmXmZpT(fsize_);
	shift_fwd(fmXmZpT,fmXmZ,TDIR);	mult_col(fmXmZpT,53);
      }
      {// (-,0,-,-)
	FermionField fmXmZmT(fsize_);
	shift_bwd(fmXmZmT,fmXmZ,TDIR);	mult_col(fmXmZmT,61);
      }
    }
    {// (-,0,0,+)
      FermionField fmXpT(fsize_);
      shift_fwd(fmXpT,fmX,TDIR);  mult_col(fmXpT,17);
    }
    {// (-,0,0,-)
      FermionField fmXmT(fsize_);
      shift_bwd(fmXmT,fmX,TDIR);  mult_col(fmXmT,19);
    }
  }
  {// (0,+,0,0)
    FermionField fpY(fsize_);
    shift_fwd(fpY,f,YDIR);  mult_col(fpY,2); 
    {// (0,+,+,0)
      FermionField fpYpZ(fsize_);
      shift_fwd(fpYpZ,fpY,ZDIR); mult_col(fpYpZ,20);
      {// (0,+,+,+)
	FermionField fpYpZpT(fsize_);
	shift_fwd(fpYpZpT,fpYpZ,TDIR); 	mult_col(fpYpZpT,38);
      }
      {// (0,+,+,-)
	FermionField fpYpZmT(fsize_);
	shift_bwd(fpYpZmT,fpYpZ,TDIR); 	mult_col(fpYpZmT,46);
      }
    }
    {// (0,+,-,0)
      FermionField fpYmZ(fsize_);
      shift_bwd(fpYmZ,fpY,ZDIR); mult_col(fpYmZ,22);
      {// (0,+,-,+)
	FermionField fpYmZpT(fsize_);
	shift_fwd(fpYmZpT,fpYmZ,TDIR); mult_col(fpYmZpT,54);
      }	
      {// (0,+,-,-)
	FermionField fpYmZmT(fsize_);
	shift_bwd(fpYmZmT,fpYmZ,TDIR);	mult_col(fpYmZmT,62);
      }
    }
    {// (0,+,0,+)
      FermionField fpYpT(fsize_);
      shift_fwd(fpYpT,fpY,TDIR);  mult_col(fpYpT,24);
    }
    {// (0,+,0,-)
      FermionField fpYmT(fsize_);
      shift_bwd(fpYmT,fpY,TDIR);  mult_col(fpYmT,26);
    }
  }
  {// (0,-,0,0)
    FermionField fmY(fsize_);
    shift_bwd(fmY,f,YDIR);  mult_col(fmY,3);
    {// (0,-,+,0)
      FermionField fmYpZ(fsize_);
      shift_fwd(fmYpZ,fmY,ZDIR);  mult_col(fmYpZ,21);
      {// (0,-,+,+)
	FermionField fmYpZpT(fsize_);
	shift_fwd(fmYpZpT,fmYpZ,TDIR);	mult_col(fmYpZpT,39);
      }
      {// (0,-,+,-)
	FermionField fmYpZmT(fsize_);
	shift_bwd(fmYpZmT,fmYpZ,TDIR);	mult_col(fmYpZmT,47);
      }
    }
    {// (0,-,-,0)
      FermionField fmYmZ(fsize_);
      shift_bwd(fmYmZ,fmY,ZDIR);  mult_col(fmYmZ,23);
      {// (0,-,-,+)
	FermionField fmYmZpT(fsize_);
	shift_fwd(fmYmZpT,fmYmZ,TDIR);	mult_col(fmYmZpT,55);
      }
      {// (0,-,-,-)
	FermionField fmYmZmT(fsize_);
	shift_bwd(fmYmZmT,fmYmZ,TDIR);	mult_col(fmYmZmT,63);
      }
    }
    {// (0,-,0,+)
      FermionField fmYpT(fsize_);
      shift_fwd(fmYpT,fmY,TDIR);  mult_col(fmYpT,25);
    }
    {// (0,-,0,-)
      FermionField fmYmT(fsize_);
      shift_bwd(fmYmT,fmY,TDIR);  mult_col(fmYmT,27); 
    }
  }
  {// (0,0,+,0)
    FermionField fpZ(fsize_);
    shift_fwd(fpZ,f,ZDIR);  mult_col(fpZ,4);    

    {// (0,0,+,+)
      FermionField fpZpT(fsize_);
      shift_fwd(fpZpT,fpZ,TDIR);  mult_col(fpZpT,28);
    }
    {// (0,0,+,-)
      FermionField fpZmT(fsize_);
      shift_bwd(fpZmT,fpZ,TDIR);  mult_col(fpZmT,30);
    }
  }
  {// (0,0,-,0)
    FermionField fmZ(fsize_);
    shift_bwd(fmZ,f,ZDIR);  mult_col(fmZ,5);

    {// (0,0,-,+)
      FermionField fmZpT(fsize_);
      shift_fwd(fmZpT,fmZ,TDIR);  mult_col(fmZpT,29);
    }
    {// (0,0,-,-)
      FermionField fmZmT(fsize_);      
      shift_bwd(fmZmT,fmZ,TDIR);  mult_col(fmZmT,31);
    }    
  }
  {// (0,0,0,+)
    FermionField fpT(fsize_);
    shift_fwd(fpT,f,TDIR); mult_col(fpT,6);
  }
  {// (0,0,0,-)
    FermionField fmT(fsize_);
    shift_bwd(fmT,f,TDIR); mult_col(fmT,7);
  }
}


