//---------------------------------------------------------------------
// eigenProc_Zolotarev.cpp
//---------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

// Forward declaration of Jacobi_elliptic
void Jacobi_elliptic(double uu, double emmc, 
		     double &sn, double &cn, double &dn);

double sign_Zolotarev(double x, const vector<double>& bl,
		      const vector<double>& cl){
//  cl[2*Np], bl[Np]: coefficients of rational approx.

  int Np = bl.size();
  double x2R = 0.0;
  for(int l=0; l<Np; l++) x2R += bl[l]/(x*x +cl[2*l]);
  return x*x2R*(x*x+cl[2*Np-1]);
}

void poly_Zolotarev(vector<double>& bl,vector<double>& cl,
		    double bmax,double& UK){
//    Return the coefficients of Zolotarev's approximation
//    of sign function to rational function.
//
//    Np: number of poles (2*Np is number of cl)
//    bmax: range of argument [1,bmax]
//    UK: complete elliptic integral of the 1st kind
//    cl[2*Np],bl[Np]: coefficients of Zolotarev rational approx.

  int Np = bl.size();
  int Nprec = 14;
  int imax = 20;

  double rk = sqrt(1.0 -1.0/(bmax*bmax));
  double emmc =1.0 -rk*rk;
  printf("emmc = %22.14e\n", emmc);

  double sn, cn, dn;

  // Determination of K
  double Dsr = 10.0;
  double u = 0.0;
  for(int iprec=0; iprec<Nprec+1; iprec++){
    Dsr = Dsr*0.1;
    for(int i=0; i<imax; i++){
      u += Dsr;
      Jacobi_elliptic(u,emmc,sn,cn,dn);
      printf(" %22.14e %22.14e %22.14e\n", u, sn, cn);
      if(cn < 0.0) break;
      if(i== imax-1) printf("Something wrong in setting Zolotarev\n");
    }
    u -= Dsr;
  }
  UK = u;

  Jacobi_elliptic(UK,emmc,sn,cn,dn);
  printf(" %22.14e %22.14e %22.14e\n", UK, sn, cn);

// Determination of c_l
  double FK = UK/(2.0*Np +1.0);

  for(int l=0; l<2*Np; l++){
    u = FK*(l+1.0);
    Jacobi_elliptic(u,emmc,sn,cn,dn);
    cl[l] = sn*sn/(1.0-sn*sn);
  }

// Determination of b_l
  double d0 = 1.0;
  for(int l=0; l<Np; l++) d0 *= (1.0+cl[2*l])/(1.0+cl[2*l+1]);

  for(int l=0; l<Np; l++){
    bl[l] = d0;
    for(int i=0; i<Np; i++){
      if(i < Np-1) bl[l] *= cl[2*i+1]-cl[2*l];
      if(i != l)   bl[l] /= cl[2*i  ]-cl[2*l];
    }
  }

// Correction
  int Nj = 10000;
  double Dx = (bmax-1.0)/Nj;

  double d1 = 0.0;
  double d2 = 2.0;

  for(int jx=0; jx <= Nj; jx++){
    double x = 1.0 + Dx*jx;
    double sgnx = sign_Zolotarev(x,bl,cl);
    if(fabs(sgnx) > d1) d1 = sgnx;
    if(fabs(sgnx) < d2) d2 = sgnx;
  }
  printf(" |sgnx|_up = %8.4e \n", d1-1.0);
  printf(" |sgnx|_dn = %8.4e \n", 1.0-d2);
  
  double dmax = d1-1.0;
  if (dmax < 1.0-d2) dmax = 1.0-d2;
  printf(" Amp(1-|sgn|) = %16.8e \n", dmax);

  return;
}

/* ****************************************************** */
void Jacobi_elliptic(double uu,double emmc,double &sn,double &cn,double &dn){
//  Return the Jacobi elliptic functions sn(u,kc),
//  cn(u,kc), and dn(u,kc), where kc = 1 - k^2, and
//  arguments:
//    uu: u, emmc: kc.
//
//  This routine is from
//  W.H.Press et al., Numerical Recipes in Fortran 2nd ed.
//  (Cambridge Univ. Press, 1986,1992).
//                         Typed by H. Matsufuru 1 May 2005

  double const CA = 0.00000001;
  //          (The accuracy is the square of CA.)
  double a,b,c,d,em[13],en[13];
  int bo;
  //  LOGICAL bo

  // main
  double emc = emmc;
  double u = uu;

  if(emc != 0.0){
    //    bo = (emc.LT.0.D0)   // how to do ?
    int bo = 0;
    if(emc < 0.0) bo = 1;

    if(bo != 0){
      d = 1.0 - emc;
      emc = - emc/d;
      d = sqrt(d);
      u = d*u;
    }
    a  = 1.0;
    dn = 1.0;

    int l;
    for(int i=0; i<13; i++){
      l = i;
      em[i] = a;
      emc = sqrt(emc);
      en[i] = emc;
      c = 0.5*(a + emc);
      if(fabs(a-emc) <= CA*a) goto continued;
      emc = a * emc;
      a = c;
    }
    
  continued:
    u = c*u;
    sn = sin(u);
    cn = cos(u);

    if(sn == 0.0) goto continued2;

    a = cn/sn;
    c = a * c;
    for(int ii=l; ii >= 0; --ii){
      b = em[ii];
      a = c * a;
      c = dn * c;
      dn = (en[ii] + a)/(b + a);
      a = c/b;
    }

    a = 1.0/sqrt(c*c + 1.0);
    if(sn < 0.0) sn = -a;
    else         sn = a;
    cn = c * sn;

  continued2:
    if(bo != 0){
      a = dn;
      dn = cn;
      cn = a;
      sn = sn/d;
    }

  }else{
    cn = 1.0/cosh(u);
    dn = cn;
    sn = tanh(u);
  }
}
/* ****************************************************** */

