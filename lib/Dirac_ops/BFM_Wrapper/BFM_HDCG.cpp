/*!
 * @file BFM_HDCG.cpp
 * @brief Declares classes for P. Boyle HDCG inverter
 * Time-stamp: <2014-08-07 16:57:08 neo>
 */

#include "BFM_HDCG.hpp"
#include "include/common_fields.hpp"
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <malloc.h>
#include <builtins.h> // builtin functions 
int MGreport;

#undef DEBUG_HDCG

int    ForceMultishiftContinue;

// Send buffers for SPI are shared across all instances
int BfmHDCGStatic::StaticInitialised;
std::vector< std::complex<double> * > BfmHDCGStatic::sendbufs_static; 
std::vector< std::complex<double> * > BfmHDCGStatic::recvbufs_static;

int AnalyseSpectrumInternal;

//////////////////////////////////////////////////
// Simple Chebyshev approximation
//////////////////////////////////////////////////
double Chebyshev::Evaluate(double x)
{
  double *c = Coeffs;

  double Tn;
  double Tnm;
  double Tnp;

  double y=( x-0.5*(hi+lo))/(0.5*(hi-lo));

  double T0=1;
  double T1=y;

  double sum;
  sum = 0.5*c[0]*T0;
  sum+= c[1]*T1;

  Tn =T1;
  Tnm=T0;
  for(int i=2;i<order;i++){
    Tnp=2*y*Tn-Tnm;
    Tnm=Tn;
    Tn =Tnp;
    sum+= Tn*c[i];
  }
  return sum;
}

void Chebyshev::Init(double _lo,double _hi,int _order, double (* func)(double))
{
  lo=_lo;
  hi=_hi;
  order=_order;
  
  if(order < 2) exit(-1);
  Coeffs = new double [order];
  for(int j=0;j<order;j++){
    double s=0;
    for(int k=0;k<order;k++){
      double y=cos(M_PI*(k+0.5)/order);
      double x=0.5*(y*(hi-lo)+(hi+lo));
      double f=func(x);
      s=s+f*cos( j*M_PI*(k+0.5)/order );
    }
    Coeffs[j] = s * 2.0/order;
  }
}
double PolynomialShift = 0.0;
double PolynomialShape (double x) { 
  return 1.0/(x+PolynomialShift);
}

// Use for analysing spectral decomposition
double Filter_0_1 (double x) { 
  if ( x < 0.1 ) return 1.0;
  else return 0.00;
}
double Filter_0_3 (double x) { 
  if ( (x>0.1) && (x < 0.3) ) return 1.0;
  else return 0.00;
}
double Filter_1_0 (double x) { 
  if ( (x>0.3) && (x < 1.0) ) return 1.0;
  else return 0.00;
}
double Filter_3_0 (double x) { 
  if ( (x>1.0) && (x < 3.0) ) return 1.0;
  else return 0.00;
}
double Filter_6_0 (double x) { 
  if ( (x>3.0) && (x < 6.0) ) return 1.0;
  else return 0.00;
}
double Filter_10_0 (double x) { 
  if ( (x>6.0) && (x < 10.0) ) return 1.0;
  else return 0.00;
}
double Filter_20_0 (double x) { 
  if ( (x>10.0) && (x < 20.0) ) return 1.0;
  else return 0.00;
}
double Filter_30_0 (double x) { 
  if ( (x>20.0) && (x < 30.0) ) return 1.0;
  else return 0.00;
}
double Filter_100_0 (double x) { 
  if ( (x>30.0) && (x < 100.0) ) return 1.0;
  else return 0.00;
}
template<class cFloat>   
template<class Float>   
void BfmHDCG<cFloat>::SpectralDecomposition(Fermion_t psi,bfm_internal<Float> *lop,const char *vname)
{
  Chebyshev Filter;
  const int nfilter = 9;
  double (*Functions[])(double)  = { Filter_0_1, Filter_0_3,Filter_1_0,Filter_3_0,Filter_6_0,Filter_10_0,Filter_20_0,Filter_30_0,Filter_100_0 };
  double FilterLo = 0;
  double FilterHi = 100;
  int    FilterOrder = 512;
  double spectrum[nfilter+1];

  Fermion_t tmp = lop->threadedAllocFermion();
  spectrum[nfilter] = lop->norm(psi);

  for ( int i=0;i<nfilter;i++ ) { 
    Filter.Init(FilterLo,FilterHi,FilterOrder,Functions[i]);
    PolyMdagMprec<Float>(lop,psi,tmp,Filter);
    spectrum[i] = lop->norm(tmp);
    Filter.End();
  }
  lop->threadedFreeFermion(tmp);
  double nn=spectrum[nfilter];
  lop->ThreadBossMessage("Vector %s Normsq = %le\n",vname,nn);
  lop->ThreadBossMessage("%s Power[0,  0.1 ] = %le\n",vname,spectrum[0]/nn);
  lop->ThreadBossMessage("%s Power[0.1,0.3 ] = %le\n",vname,spectrum[1]/nn);
  lop->ThreadBossMessage("%s Power[0.3,1.0 ] = %le\n",vname,spectrum[2]/nn);
  lop->ThreadBossMessage("%s Power[1.0,3.0 ] = %le\n",vname,spectrum[3]/nn);
  lop->ThreadBossMessage("%s Power[3.0,6.0 ] = %le\n",vname,spectrum[4]/nn);
  lop->ThreadBossMessage("%s Power[6, 10.0 ] = %le\n",vname,spectrum[5]/nn);
  lop->ThreadBossMessage("%s Power[10, 20  ] = %le\n",vname,spectrum[6]/nn);
  lop->ThreadBossMessage("%s Power[20, 30  ] = %le\n",vname,spectrum[7]/nn);
  lop->ThreadBossMessage("%s Power[30,100  ] = %le\n",vname,spectrum[8]/nn);
  lop->ThreadBossMessage("%s Power[0,inf   ] = %le\n",vname,spectrum[9]/nn);
  
}


//////////////////////////////
// Coarse space linalg
//////////////////////////////
template<class cFloat> 
void BfmHDCG<cFloat>::scale( std::vector<std::complex<cFloat> > &r,
			  std::vector<std::complex<cFloat> > &x,double a)
{
  axpby(r,x,x,a,0.0);
}
template<class cFloat> 
void BfmHDCG<cFloat>::zeroOut(std::vector<std::complex<cFloat> > &v) {
  axpby(v,v,v,0.0,0.0); // NAN propagation?
}

template<class cFloat> 
void BfmHDCG<cFloat>::zaxpy(  std::vector<std::complex<cFloat> > &r,
			   std::vector<std::complex<cFloat> > &x,
			   std::vector<std::complex<cFloat> > &y,std::complex<double> a)
{
  if ( r.size() != LocalNsubspace ) exit(0);
  if ( x.size() != LocalNsubspace ) exit(0);
  if ( y.size() != LocalNsubspace ) exit(0);
  int me,throff,thrlen;
  linop_d->thread_work(x.size(),me,thrlen,throff);
#ifdef USE_XLC_OPTIMISED_CODE
  if ( sizeof(cFloat) == sizeof(double) ){ 
    qpx_zaxpy(thrlen,
	      (double *)&r[throff],
	      (double *)&x[throff],
	      (double *)&y[throff],
	      (double *)&a);
  } else { 
    qpx_zaxpy_f(thrlen,
	      (float *)&r[throff],
	      (float *)&x[throff],
	      (float *)&y[throff],
	      (double *)&a);
  }
#else
  for(int i=throff;i<throff+thrlen;i++){
    r[i] = a*x[i]+y[i];
  }
#endif
  linop_d->thread_barrier();
}

template<class cFloat> 
void BfmHDCG<cFloat>::axpby(std::vector<std::complex<cFloat> > &r,
			   std::vector<std::complex<cFloat> > &x,
			   std::vector<std::complex<cFloat> > &y,double a,double b)
{
  int me,throff,thrlen;
  if ( r.size() != LocalNsubspace ) { fprintf(stderr,"axpby wrong vector length r %d/%d %lx %lx %lx\n",r.size(),LocalNsubspace,
					      __builtin_return_address(0),__builtin_return_address(1),__builtin_return_address(2));exit(0);}
  if ( x.size() != LocalNsubspace ) { fprintf(stderr,"axpby wrong vector length x %d/%d %lx %lx %lx\n",x.size(),LocalNsubspace,
					      __builtin_return_address(0),__builtin_return_address(1),__builtin_return_address(2));exit(0);}
  if ( y.size() != LocalNsubspace ) { fprintf(stderr,"axpby wrong vector length y %d/%d %lx %lx %lx\n",y.size(),LocalNsubspace,
					      __builtin_return_address(0),__builtin_return_address(1),__builtin_return_address(2));exit(0);}
  linop_d->thread_work(x.size(),me,thrlen,throff);
#ifdef USE_XLC_OPTIMISED_CODE
  if ( sizeof(cFloat) == sizeof(double ) ) { 
    qpx_axpby(thrlen,(double *)&r[throff],(double *)&x[throff],(double *)&y[throff],a,b);
  } else { 
    qpx_axpby_f(thrlen,(float *)&r[throff],(float *)&x[throff],(float *)&y[throff],a,b);
  }
#else
  for(int i=throff;i<throff+thrlen;i++){
    r[i] = a*x[i]+b*y[i];
  }
#endif
  linop_d->thread_barrier();
}
template<class cFloat>
void BfmHDCG<cFloat>::axpy(std::vector<std::complex<cFloat> > &r,
			std::vector<std::complex<cFloat> > &x,
			std::vector<std::complex<cFloat> > &y,double a)
{
  int me,throff,thrlen;
  if ( r.size() != LocalNsubspace ) { fprintf(stderr,"wrong vector length %d/%d %lx %lx %lx r\n",r.size(),LocalNsubspace,
					      __builtin_return_address(0),__builtin_return_address(1),__builtin_return_address(2));exit(0);}
  if ( x.size() != LocalNsubspace ) { fprintf(stderr,"wrong vector length %d/%d %lx %lx %lx x\n",x.size(),LocalNsubspace,
					      __builtin_return_address(0),__builtin_return_address(1),__builtin_return_address(2));exit(0);}
  if ( y.size() != LocalNsubspace ) { fprintf(stderr,"wrong vector length %d/%d %lx %lx %lx y\n",y.size(),LocalNsubspace,
					      __builtin_return_address(0),__builtin_return_address(1),__builtin_return_address(2));exit(0);}
  linop_d->thread_work(x.size(),me,thrlen,throff);
#ifdef USE_XLC_OPTIMISED_CODE
  if ( sizeof(cFloat) == sizeof(double) ) {
    qpx_axpy(thrlen,(double *)&r[throff],(double *)&x[throff],(double *)&y[throff],a);
  }else { 
    qpx_axpy_f(thrlen,(float *)&r[throff],(float *)&x[throff],(float *)&y[throff],a);
  }
#else
  for(int i=throff;i<throff+thrlen;i++){
    r[i] = a*x[i]+y[i];
  }
#endif
  linop_d->thread_barrier();
}
template<class cFloat> 
std::complex<double> BfmHDCG<cFloat>::innerProduct(std::vector<std::complex<cFloat> > &v1,
						std::vector<std::complex<cFloat> > &v2) 
{
  std::complex<double> ir = 0;
  int me,throff,thrlen;
  
  linop_d->thread_work(v1.size(),me,thrlen,throff);
  
  if ( v1.size() != LocalNsubspace ) exit(0);
  if ( v2.size() != LocalNsubspace ) exit(0);
    
#ifdef USE_XLC_OPTIMISED_CODE
  if ( sizeof(cFloat) == sizeof(double) ) { 
    qpx_inner((double *)&ir,thrlen,(double *)&v1[throff],(double *)&v2[throff]);
  } else { 
    qpx_inner_f((double *)&ir,thrlen,(float *)&v1[throff],(float *)&v2[throff]);
  }
#else
  for(int i=throff;i<throff+thrlen;i++){
    ir = ir + conj(v1[i])*v2[i];
  }
#endif
  
  double re = real(ir);
  double im = imag(ir);
  
  linop_d->thread_sum(re,me);
  linop_d->comm_gsum(re);
  
  linop_d->thread_sum(im,me);
  linop_d->comm_gsum(im);
  
  std::complex<double> a(re,im);
  return a;
}
template<class cFloat> 
double BfmHDCG<cFloat>::innerProductReal(std::vector<std::complex<cFloat> > &v1,
				      std::vector<std::complex<cFloat> > &v2)
{
  double ir = 0;
  int me,throff,thrlen;
  linop_d->thread_work(v1.size(),me,thrlen,throff);
#ifdef USE_XLC_OPTIMISED_CODE
  if ( sizeof(cFloat) == sizeof(double) ) { 
    ir = qpx_inner_real(thrlen,(double *)&v1[throff],(double *)&v2[throff]);
  } else {
    ir = qpx_inner_real_f(thrlen,(float *)&v1[throff],(float *)&v2[throff]);
  }
#else
  for(int i=throff;i<throff+thrlen;i++){
    ir = ir + real(v1[i])*real(v2[i]) + imag(v1[i])*imag(v2[i]);
  }
#endif
  linop_d->thread_sum(ir,me);
  linop_d->comm_gsum(ir);
  return ir;
}
template<class cFloat> 
double BfmHDCG<cFloat>::norm_vec(std::vector<std::complex<cFloat> > &v) {
  double nrm=0.0;
  int me,thrlen,throff;
  
  linop_d->thread_work(v.size(),me,thrlen,throff);
#ifdef USE_XLC_OPTIMISED_CODE
  if ( sizeof(cFloat) == sizeof(double) ) { 
    nrm = qpx_norm   (thrlen,(double *)&v[throff]);
  } else {
    nrm = qpx_norm_f (thrlen,(float *)&v[throff]);
  }
#else
  for(int i=throff;i<throff+thrlen;i++){
    nrm += real(v[i])*real(v[i])+imag(v[i])*imag(v[i]);
  }
#endif
  linop_d->thread_sum(nrm,me);
  linop_d->comm_gsum(nrm);
  return nrm;
}

// Pass a linop and solve in single
template<class cFloat> 
template<class Float>
void BfmHDCG<cFloat>::RationalSubspace(bfm_internal<Float> * rop,Fermion_t src,Fermion_t sol)
{
  const int nshift=4;
  Fermion_t sol_guess[nshift];

  double Lo = SubspaceRationalLo;
  //
  // [ 1/6(x+Lo)  - 1/2(x+2Lo) + 1/2(x+3Lo)  -1/6(x+4Lo) = Lo^3 /[ (x+1Lo)(x+2Lo)(x+3Lo)(x+4Lo) ]
  //
  // 1/(x+Lo)  - 1/(x+2 Lo)
  //
  double epsilon      = Lo/3;
  double alpha[4]     = {1.0,1.0,1.0,1.0};
  double shifts[4]    = {Lo,Lo+epsilon,Lo+2*epsilon,Lo+3*epsilon};
  //  double mresidual[4] = {3.0*SubspaceRationalResidual,SubspaceRationalResidual,SubspaceRationalResidual,SubspaceRationalResidual};
  double mresidual[4] = {SubspaceRationalResidual,SubspaceRationalResidual,SubspaceRationalResidual,SubspaceRationalResidual};

  int me = rop->thread_barrier();
  if ( !me )  SubspaceRationalLo *= 1.00;

  rop->ThreadBossLog("Rational 4th order low pass : Bandstop=%f epsilon=%f tol = %le\n",Lo,epsilon,SubspaceRationalResidual);

  for(int i=0;i<nshift;i++){
    sol_guess[i] = rop->threadedAllocFermion();
  }
  int power=1;
  for(int i=0;i<power;i++){
    int single = 0;
    rop->CGNE_prec_MdagM_multi_shift(sol_guess,
				     src,
				     shifts,
				     alpha,
				     nshift,
				     mresidual,
				     single,
				     ForceMultishiftContinue);
    

    rop->axpby(sol,sol_guess[0],sol_guess[1],1.0/6.0,-1.0/2.0);
    rop->axpy(sol,sol_guess[2],sol,1.0/2.0);
    rop->axpy(sol,sol_guess[3],sol,-1.0/6.0);

    rop->axpy(src,sol,sol,0.0);
  }

  for(int i=0;i<nshift;i++){
    rop->threadedFreeFermion(sol_guess[i]);
  }

}
template<class cFloat> 
void BfmHDCG<cFloat>::SinglePrecSubspace(void)
{
  int threads   = linop_f->threads;
  subspace_f = (Fermion_t *) malloc(Nvec*sizeof(Fermion_t));
  for(int i=0;i<Nvec;i++){
    subspace_f[i] = linop_f->allocFermion();
  }

#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      for(int i=0;i<Nvec;i++){
	linop_f->precisionChange(subspace_d[i],subspace_f[i],DoubleToSingle,1);
      }
    }
  }
}



template<class cFloat> 
template<class Float>
void BfmHDCG<cFloat>::RelaxSubspace(bfm_internal<Float> *rop)
{
  int Ls=N5;
  int rLs=rop->Ls;
  int Nvol = CommonPrms::instance()->Nvol();

  //  multi1d<LatticeFermion> gauss(rLs);
  //multi1d<LatticeFermion> solution(rLs);

  FermionField gauss(Nvol*Ls);
  FermionField solution(Nvol*Ls);


  if ( rop->SPIcomms() ) rop->comm_init();
  Fermion_t sol = rop->allocFermion();
  Fermion_t src = rop->allocFermion();

  subspace_d   = (Fermion_t *) malloc(2*Nvec*sizeof(Fermion_t));
  subspace_r   = (Fermion_t *) malloc(2*Nvec*sizeof(Fermion_t)); // Don't orthogonalise this guy

  for(int v=0;v<Nvec;v++) {
    subspace_d[v] = linop_d->allocFermion();
    subspace_r[v] = linop_d->allocFermion();

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop_d->threads;i++) {
	linop_d->set_zero(subspace_r[v]);
	linop_d->set_zero(subspace_d[v]);
      }
    }
  }

  uint64_t t_orthog=0;
  uint64_t t0;
  for(int v=0;v<Nvec;v++) {

    ////////////////////////////////
    // Prepare the source
    ////////////////////////////////
    rop->UniformRandom(src);
    rop->exportFermion(gauss,src,1);
    gauss[0]   = chiralProjectPlus(gauss[0]);
    gauss[rLs-1]= chiralProjectMinus(gauss[rLs-1]);
    for(int s=1;s<rLs-1;s++) gauss[s]=zero;

    //////////////////////////////////////////
    // Project out the existing subspace
    //////////////////////////////////////////
    if ( 1 ) { 
      Fermion_t tmp1= linop_d->allocFermion();    
      Fermion_t tmp2= linop_d->allocFermion();    

      linop_d->master_fill(tmp1,0.0);
      linop_d->importFermion(gauss[0],tmp1,1,0);
      linop_d->importFermion(gauss[rLs-1],tmp1,1,Ls-1);

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<rop->threads;i++) {
	ProjectToSubspace(tmp1,PleftProj);     
	PromoteFromSubspace(PleftProj,tmp2);  
	linop_d->axpy(tmp1,tmp2,tmp1,-1.0);
      }
    }

      for(int s=0;s<rLs/2;s++) {
	linop_d->exportFermion(gauss[s],tmp1,1,s);
	linop_d->exportFermion(gauss[rLs-1-s],tmp1,1,Ls-1-s);
      }

      linop_d->freeFermion(tmp1);
      linop_d->freeFermion(tmp2);
    }

    rop->importFermion(gauss,src,1);
    
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<rop->threads;i++) {
	RationalSubspace(rop,src,sol);
      }
    }

    for(int s=0;s<rLs;s++) gauss[s]=zero;

    rop->exportFermion(gauss,sol,1);    // + rest

    s_min=0; s_max=Ls-1;

    for(int s=0;s<rLs/2;s++) {
      linop_d->importFermion(gauss[s],subspace_d[v],1,s);
      linop_d->importFermion(gauss[rLs-1-s],subspace_d[v],1,Ls-1-s);
    }

    //    linop_d->residual = SubspaceRationalRefineResidual;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop_d->threads;i++) {
	linop_d->copy(subspace_r[v],subspace_d[v]);
      }
    }

    linop_d->BossLog("HDCG: Relax Subspace Orthogonalising %ld ms\n", t_orthog/MHz()/1000); fflush(stdout);
    t0 = GetTimeBase();
    OrthogonaliseSubspace();
    t_orthog+=GetTimeBase()-t0;
  }

  rop->freeFermion(src);
  rop->freeFermion(sol);
}



template void BfmHDCG<float>::RelaxSubspace<float>(bfm_internal<float> *rop);
template void BfmHDCG<float>::RationalSubspace<float>(bfm_internal<float> * rop,Fermion_t src,Fermion_t sol);
template void BfmHDCG<float>::RelaxSubspace<double>(bfm_internal<double> *rop);
template void BfmHDCG<float>::RationalSubspace<double>(bfm_internal<double> * rop,Fermion_t src,Fermion_t sol);

template void BfmHDCG<double>::RelaxSubspace<float>(bfm_internal<float> *rop);
template void BfmHDCG<double>::RationalSubspace<float>(bfm_internal<float> * rop,Fermion_t src,Fermion_t sol);
template void BfmHDCG<double>::RelaxSubspace<double>(bfm_internal<double> *rop);
template void BfmHDCG<double>::RationalSubspace<double>(bfm_internal<double> * rop,Fermion_t src,Fermion_t sol);

template<class cFloat> 
void BfmHDCG<cFloat>::FreeDoubleSubspace (void)
{
  for (int v=0;v<Nvec;v++){
    linop_d->freeFermion(subspace_d[v]);
    linop_d->freeFermion(subspace_r[v]);
  }
  free(subspace_r);
  free(subspace_d);
}

template<class cFloat> 
void BfmHDCG<cFloat>::FreeSingleSubspace (void)
{
  for (int v=0;v<Nvec;v++){
    linop_f->freeFermion(subspace_f[v]);
  }
  free(subspace_f);
}

template<class cFloat> 
template<class oFloat> 
void BfmHDCG<cFloat>::CloneSubspace(BfmHDCG<oFloat> &ref)
{
  ///////////////////////////////////////
  // Clone fine subspace vectors
  ///////////////////////////////////////
  subspace_d = ref.subspace_d;
  subspace_f = ref.subspace_f;
  s_min = ref.s_min;
  s_max = ref.s_max;

  ///////////////////////////////////////
  // Clone ldop representation
  ///////////////////////////////////////
  int sz = ref.Asparse.size();
  Asparse.resize(sz);
  AsparseSingle.resize(sz);
  for(int i=0;i<sz;i++){
    Asparse[i]=ref.Asparse[i];
    AsparseSingle[i]=ref.Asparse[i];
  }

  sz = ref.AsparseDiagInv.size();
  AsparseDiagInv.resize(sz);
  for(int i=0;i<sz;i++){
    AsparseDiagInv[i]=ref.AsparseDiagInv[i];
  }

  ///////////////////////////////////////
  // Clone deflation vectors
  ///////////////////////////////////////
  int nv  =  LdopDeflationBasisSize = ref.LdopDeflationBasisSize;
  LdopDeflationIsDiagonal = ref.LdopDeflationIsDiagonal;

  int nev=ref.LdopDeflationEvals.size();
  LdopDeflationEvals.resize(nev);
  LdopDeflationInvEvals.resize(nev);
  for(int i=0;i<nev;i++){
    LdopDeflationEvals[i]=ref.LdopDeflationEvals[i];
    LdopDeflationInvEvals[i]=ref.LdopDeflationInvEvals[i];
  }

  //resize
  LdopDeflationBasis.resize(LdopDeflationBasisSize);
  LdopDeflationAv.resize(LdopDeflationBasisSize);
  LdopDeflationMatrix.resize(nv*nv);
  LdopDeflationInverseMatrix.resize(nv*nv);
  
  for(int i=0;i<nv*nv;i++){
    LdopDeflationMatrix[i] = ref.LdopDeflationMatrix[i];
    LdopDeflationInverseMatrix[i] = ref.LdopDeflationInverseMatrix[i];
  }

  for(int i=0;i<nv;i++){
    LdopDeflationBasis[i].resize(LocalNsubspace);
    LdopDeflationAv[i].resize(LocalNsubspace);
    for(int j=0;j<LocalNsubspace;j++){
      LdopDeflationBasis[i][j] = ref.LdopDeflationBasis[i][j];
      LdopDeflationAv[i][j] = ref.LdopDeflationAv[i][j];
    }
  }
  // Stencil info
  for(int i=0;i<Nball;i++){
    StencilNonZero[i]=ref.StencilNonZero[i];
    StencilReord[i]  =ref.StencilReord[i];
    StencilDepth[i]  =ref.StencilDepth[i];
  }
  for(int deep=0;deep<5;deep++){
    StencilDepthCount[deep] = ref.StencilDepthCount[deep];
  }

  DeflKrylovProj.resize(nv);
  DeflKrylovMss.resize(nv);
  DeflKrylovGsum.resize(nv);

}

template void BfmHDCG<double>::CloneSubspace(BfmHDCG<double> &);
template void BfmHDCG<float>::CloneSubspace(BfmHDCG<double> &);
template void BfmHDCG<double>::CloneSubspace(BfmHDCG<float> &);
template void BfmHDCG<float>::CloneSubspace(BfmHDCG<float> &);

template void BfmHDCG<double>::FreeDoubleSubspace(void);
template void BfmHDCG<float>::FreeDoubleSubspace(void);
template void BfmHDCG<double>::FreeSingleSubspace(void);
template void BfmHDCG<float>::FreeSingleSubspace(void);

template<class cFloat> 
BfmHDCG<cFloat>::BfmHDCG(int _N5,
			   int _Nvec,
			   int _BlockSize[5],
			   int _QuadrantSize[4],
			   bfm_internal<double> * _linop_d,
			   bfm_internal<float> * _linop_f
			   ) 
{
  AnalyseSpectrum=0;
  Nvec=_Nvec;
  N5=_N5;
  LdopDeflationIsDiagonal=0;
  s_min=0;
  s_max=N5-1;

  LdopDeflationBasisSize = 0;

  linop_d=_linop_d;
  linop_f=_linop_f;

  for(int mu=0;mu<5;mu++){
    BlockSize[mu]=_BlockSize[mu];
    linop_d->BossLog("BlockSize[%d]=%d\n",mu,BlockSize[mu]);
  }
  for(int mu=0;mu<4;mu++){
    QuadrantSize[mu]=_QuadrantSize[mu];
  }
  linop_d=_linop_d;
  linop_f=_linop_f;

  // Block layout
  for(int mu=0;mu<4;mu++){
    BlockSize[mu] = _BlockSize[mu];
    glatt[mu]     = Layout::lattSize()[mu];
    global_nb[mu] = glatt[mu]/BlockSize[mu];
    if ( global_nb[mu]*BlockSize[mu] !=  Layout::lattSize()[mu] ) {
      printf("global block size error\n");
      exit(0);
    }
    local_nb[mu] = QDP::Layout::subgridLattSize()[mu]/BlockSize[mu];
    if ( local_nb[mu]*BlockSize[mu] !=  QDP::Layout::subgridLattSize()[mu] ) {
      printf("local block size error\n");
      exit(0);
    }
  }

  // Fifth dimension
  NblockS= N5/BlockSize[4];  
  global_nb[4]=NblockS;
  local_nb[4] =NblockS;
  if ( NblockS*BlockSize[4]!=N5 ) {
    QDP_error_exit("BlockSize err");
  }

  // Dimensions of subspace blocks (globally)
  Nblock4= global_nb[0]*global_nb[1]*global_nb[2]*global_nb[3];
  Nblock =NblockS*Nblock4;

  Nball=1;
  for(int d=0;d<5;d++){
    if ( global_nb[d] == 1 ) {
      stencil_lo[d] = 0;
      stencil_hi[d] = 0;
      stencil_size[d]= 1;
    } else if ( global_nb[d] == 2 ) {
      stencil_lo[d] = -1;
      stencil_hi[d] = 0;
      stencil_size[d]= 2;
    } else if ( global_nb[d] > 2 ) {
      stencil_lo[d] = -1;
      stencil_hi[d] =  1;
      stencil_size[d]= 3;
    }
    Nball=Nball*stencil_size[d];  // Keep matrix representation local
  }
  StencilNonZero.resize(Nball);
  StencilDepth.resize(Nball);
  StencilReord.resize(Nball);
  StencilDepthCount.resize(5);

  // Dimensions of subspace blocks (locally)
  LocalNblock4= local_nb[0]*local_nb[1]*local_nb[2]*local_nb[3];
  LocalNblock = NblockS*LocalNblock4;
  Nsubspace=Nblock*Nvec;
  LocalNsubspace=LocalNblock*Nvec;

  ProjectInnerProduct.resize(LocalNblock);
  PromoteCoeff.resize(LocalNblock);
  ComputeProj.resize(Nball);
  for(int b=0;b<Nball;b++){
    ComputeProj[b].resize(LocalNsubspace);
  }
  phases.resize(LocalNblock);

  slatt[0]=QDP::Layout::subgridLattSize()[0];
  slatt[1]=QDP::Layout::subgridLattSize()[1];
  slatt[2]=QDP::Layout::subgridLattSize()[2];
  slatt[3]=QDP::Layout::subgridLattSize()[3];
  slatt[4]=N5;

  ncoor[0] = QDP::Layout::nodeCoord()[0];
  ncoor[1] = QDP::Layout::nodeCoord()[1];
  ncoor[2] = QDP::Layout::nodeCoord()[2];
  ncoor[3] = QDP::Layout::nodeCoord()[3];
  ncoor[4] = 0;

  int vol   = slatt[0]*slatt[1]*slatt[2]*slatt[3]*slatt[4];
  int cbvol = vol/2;

  inner_reduce.resize(cbvol*Nvec);

  Asparse.resize      (LocalNblock*Nball*Nvec*Nvec);
  AsparseSingle.resize(LocalNblock*Nball*Nvec*Nvec);
  AsparseDiagInv.resize(LocalNblock*Nvec*Nvec);

  // Preallocate vectors needed in threaded routines
  PleftProj.resize(LocalNsubspace) ;
  PleftMss_proj.resize(LocalNsubspace) ;
  Krylov_p.resize(LocalNsubspace);
  Krylov_r.resize(LocalNsubspace);
  Krylov_mu.resize(LocalNsubspace);
  Krylov_Ap.resize(LocalNsubspace);
  Krylov_Ar.resize(LocalNsubspace);
  Krylov_Amu.resize(LocalNsubspace);
  Krylov_tmp.resize(LocalNsubspace);
  Krylov_Atmp.resize(LocalNsubspace*Nball);
  KrylovNest_p.resize(LocalNsubspace);
  KrylovNest_Ap.resize(LocalNsubspace);
  KrylovNest_r.resize(LocalNsubspace);

  block_id.resize(2);
  quadrant_id.resize(2);
  local_block_id.resize(2);
  for(int cb=0;cb<2;cb++){
    block_id[cb].resize(cbvol);
    local_block_id[cb].resize(cbvol);
    quadrant_id[cb].resize(cbvol);
  }

  std::vector<std::vector< std::vector<int> > > coordinate(5);
  for(int mu=0;mu<5;mu++){
    coordinate[mu].resize(2);
    coordinate[mu][0].resize(cbvol); 
    coordinate[mu][1].resize(cbvol);
  }
    

  int nq[4];
  int qs[4];
  for(int mu=0;mu<4;mu++){
    nq[mu] = BlockSize[mu]/QuadrantSize[mu];
    qs[mu] = QuadrantSize[mu];
    if( nq[mu]*QuadrantSize[mu] != BlockSize[mu] ) { 
      linop_d->Error("Quadrant does not divide block\n");
      nq[mu]++;
    }
  }
  Nquadrant = nq[0]*nq[1]*nq[2]*nq[3];

  for(int site=0;site<vol;site++){
      
    int ss=site;
    int x[5];
    x[0]=ss%slatt[0];    ss=ss/slatt[0];
    x[1]=ss%slatt[1];    ss=ss/slatt[1];
    x[2]=ss%slatt[2];    ss=ss/slatt[2];
    x[3]=ss%slatt[3];    ss=ss/slatt[3];
    x[4]=ss%slatt[4];
	
    int sp;
    if ( linop_d->precon_5d ) sp=x[4];
    else sp=0;
      
    int cb = ((x[0]+x[1]+x[2]+x[3]+sp)&0x1) ;
    int cbsite = linop_d->bagel_idx5d(x,x[4],0,0,1,1); cbsite=cbsite>>1; //wipes out reim factor of 2

    int b[5];
    int lb[5];
    int xb[5];

    for(int mu=0;mu<5;mu++){
      coordinate[mu][cb][cbsite] = x[mu]+ncoor[mu]*slatt[mu]; //Global coordinate
      b[mu]=coordinate[mu][cb][cbsite]/BlockSize[mu];         //Global block
      lb[mu]=x[mu]/BlockSize[mu];                         //Local block
      xb[mu]=x[mu]%BlockSize[mu];
    }

    // Needed for constructing little dirac op by picking subvectors
    block_id[cb][cbsite] = b[0]+global_nb[0]*(b[1]+global_nb[1]*(b[2]+global_nb[2]*(b[3]+global_nb[3]*b[4])));

    quadrant_id[cb][cbsite] =                   (xb[0]/QuadrantSize[0]) ;
    quadrant_id[cb][cbsite]+=             nq[0]*(xb[1]/QuadrantSize[1]) ;
    quadrant_id[cb][cbsite]+=       nq[1]*nq[0]*(xb[2]/QuadrantSize[2]) ;
    quadrant_id[cb][cbsite]+= nq[2]*nq[1]*nq[0]*(xb[3]/QuadrantSize[3]) ;

    //      ball_id_cb[cb][cbsite]= ball_id[cb][cbsite] + ball_cb[cb][cbsite]*32;

    local_block_id[cb][cbsite]=lb[0]+local_nb[0]*(lb[1]+local_nb[1]*(lb[2]+local_nb[2]*(lb[3]+local_nb[3]*lb[4])));
    
  }

  HaloInit();

  LdopM1control= atoi(getenv("LDOPM1_CONTROL"));
  LdopM1Lo   = atof(getenv("LDOPM1_LO"));
  LdopM1Hi   = atof(getenv("LDOPM1_HI"));
  LdopM1resid= atof(getenv("LDOPM1_RESID"));
  LdopM1iter = atoi(getenv("LDOPM1_ITER"));
  SloppyComms= atoi(getenv("SLOPPY_COMMS"));
  ForceMultishiftContinue=1;

  linop_d->BossLog("HDCG Initialised Little Dirac Operator\n");
}


template<class cFloat> 
int BfmHDCG<cFloat>::GetLocalBlockIdx(int d[5])
{
  return d[0] + local_nb[0]*(d[1]+local_nb[1]*(d[2]+local_nb[2]*(d[3]+local_nb[3]*d[4])));
}


template<class cFloat> 
void BfmHDCG<cFloat>::GlobalToLocalBlock(int gb[5],int &proc,int lb[5])
{
  multi1d<int> gcoor(4);
  for(int mu=0;mu<5;mu++)   lb[mu]=gb[mu]%local_nb[mu];
  for(int mu=0;mu<4;mu++)   gcoor[mu] = gb[mu]*BlockSize[mu];
  proc = QDP::Layout::nodeNumber(gcoor);  
}
template<class cFloat> 
int BfmHDCG<cFloat>::ReverseDirection(int mu)
{
  int d[5];
  int rd[5];
  d[0] = mu % stencil_size[0]; mu = mu/stencil_size[0];
  d[1] = mu % stencil_size[1]; mu = mu/stencil_size[1];
  d[2] = mu % stencil_size[2]; mu = mu/stencil_size[2];
  d[3] = mu % stencil_size[3]; mu = mu/stencil_size[3];
  d[4] = mu % stencil_size[4]; 
  for(int dir=0;dir<5;dir++){
    if ( stencil_size[dir]==2 )
      rd[dir] = stencil_size[dir]-d[dir];
    else 
      rd[dir] = d[dir];
  }
  mu = rd[4]+stencil_size[4]*(rd[3]+stencil_size[3]*(rd[2]+stencil_size[2]*(rd[1]+stencil_size[1]*rd[0])));
  return mu;
}

template<class cFloat> 
void BfmHDCG<cFloat>::GetDelta(int delta[5],int _ballidx,int rev)
{
  int d[5];
  int ballidx=0;

  for(d[0]=stencil_lo[0];d[0]<=stencil_hi[0];d[0]++){
    for(d[1]=stencil_lo[1];d[1]<=stencil_hi[1];d[1]++){
      for(d[2]=stencil_lo[2];d[2]<=stencil_hi[2];d[2]++){
	for(d[3]=stencil_lo[3];d[3]<=stencil_hi[3];d[3]++){
	  for(d[4]=stencil_lo[4];d[4]<=stencil_hi[4];d[4]++){
	    if ( _ballidx == ballidx ) {
	      if( rev ) 
		for(int mu=0;mu<5;mu++) {
		  if ( stencil_size[mu]==2 ) delta[mu] = -d[mu];
		  else delta[mu] = d[mu];
		}
	      else for(int mu=0;mu<5;mu++) delta[mu] = d[mu];
	      return;
	    }
	    ballidx++;
	  }
	}
      }
    }
  }
  exit(0);
}

template<class cFloat> 
int BfmHDCG<cFloat>::NeighbourBlockId(int myblock,int _mu, int &me,
				   int &from,int &xblock,int &to,int rev)
{
  int idx = myblock;
  
  int b[5];// Local block coord of this site
  b[0] = idx%local_nb[0]; idx = idx/local_nb[0];
  b[1] = idx%local_nb[1]; idx = idx/local_nb[1];
  b[2] = idx%local_nb[2]; idx = idx/local_nb[2];
  b[3] = idx%local_nb[3]; idx = idx/local_nb[3];
  b[4] = idx%local_nb[4];

  int gb[5]; // Global block coordinate
  for(int mu=0;mu<5;mu++)    gb[mu] = b[mu]+ncoor[mu]*local_nb[mu];

  int lwb[5]; 
  int seeker[5]; 

  int delta[5]; 
  GetDelta(delta,_mu,rev); 

  // Block plus delta
  int gbpd[5];
  for(int mu=0;mu<5;mu++){
    // Periodic wrap at global lattice boundaries
    gbpd[mu] = gb[mu]+delta[mu];
    if ( gbpd[mu] >= global_nb[mu] ) gbpd[mu]-=global_nb[mu];
    if ( gbpd[mu] <  0             ) gbpd[mu]+=global_nb[mu];

    // Locally wrap internal to node. Someone must be also seeking this point
    lwb[mu]= b[mu]+delta[mu];
    if ( lwb[mu] >= local_nb[mu] ) lwb[mu]-=local_nb[mu];
    if ( lwb[mu] <  0            ) lwb[mu]+=local_nb[mu];
    
    // Who would be looking for that
    seeker[mu] = lwb[mu]+ncoor[mu]*local_nb[mu] - delta[mu];
    if ( seeker[mu]>= global_nb[mu] ) seeker[mu]-=global_nb[mu];
    if ( seeker[mu]<0 )               seeker[mu]+=global_nb[mu];

  }

  int fb[5];
  int tb[5];
  int mb[5];
  GlobalToLocalBlock(gbpd,from,fb);
  GlobalToLocalBlock(seeker,to,tb);

  GlobalToLocalBlock(gb,me,mb);
  int mecheck= QDP::Layout::nodeNumber();
  if( me != mecheck) QDP_error_exit("oops");

  xblock   =GetLocalBlockIdx(fb);
  int lwblock =GetLocalBlockIdx(lwb);
  if ( xblock != lwblock ) QDP_error_exit("block oops");

  int meblock= GetLocalBlockIdx(mb);
  if ( meblock != myblock ) {
    cout << "logic error meblock = "<<meblock<<" by myblock="<< myblock<<endl;
    exit(0);
  }

}



template<class cFloat> 
void BfmHDCG<cFloat>::HaloEnd(void)
{

  /*Original QMP
  int Nmsg = tproc.size();
  for(int m=0;m<Nmsg;m++) {
    QMP_free_msghandle(xmit_msghandle[m]);
    QMP_free_msghandle(recv_msghandle[m]);
    QMP_free_msgmem(xmit_msgmem[m]);
    QMP_free_msgmem(recv_msgmem[m]);
#ifndef BFM_BGQ_SPI_DSLASH
    bfm_free(sendbufs[m]);
    bfm_free(recvbufs[m]);
#endif
  }
*/
  int Nmsg = tproc.size();
  for(int m=0;m<Nmsg;m++) {
#ifndef BFM_BGQ_SPI_DSLASH
    bfm_free(sendbufs[m]);
    bfm_free(recvbufs[m]);
#endif
  }

}
template<class cFloat> 
void BfmHDCG<cFloat>::HaloInit(void)
{
  ////////////////////////////////////////////////////////////
  // Map the hopping mu within ball to a depth (number hops)
  ////////////////////////////////////////////////////////////
  for(int b=0;b<Nball;b++){
    int    imom[5];
    GetDelta(imom,b,0);  
    int dist = abs(imom[0])+abs(imom[1])+abs(imom[2])+abs(imom[3])+abs(imom[4]);
    StencilDepth[b] = dist;
  }
  ////////////////////////////////////////////////////////////
  // Enumerate these in order of increasing stencil depth
  ////////////////////////////////////////////////////////////
  for(int deep=0;deep<5;deep++){

    if ( deep == 0 ) StencilDepthCount[0] = 0;
    else  StencilDepthCount[deep] = StencilDepthCount[deep-1];

    for(int mu=0;mu<Nball;mu++){
      if ( StencilDepth[mu] == deep ) {
	StencilReord[StencilDepthCount[deep]] = mu;
	StencilDepthCount[deep]++;
      }
    }
    linop_d->BossLog("%d-hop matrix has Nstencil %d\n",deep,StencilDepthCount[deep]);
  }

  tproc.resize(0);
  fproc.resize(0);

  sendbufs.resize(Nball);
  recvbufs.resize(Nball);
  
  sendbuf_bytes.resize(Nball);
  recvbuf_bytes.resize(Nball);
  for(int mu=0;mu<Nball;mu++){
    sendbuf_bytes[mu]=0;
    recvbuf_bytes[mu]=0;
  }

  sendbuf_depth_bytes.resize(5);
  sbuf_gather_depth_count.resize(5);
  for(int d=0;d<5;d++){
    sendbuf_depth_bytes[d].resize(Nball);
    sbuf_gather_depth_count[d] = 0;
  }

  ////////////////////////////////////////////////
  // Which buffer (or none) each block and dir maps to 
  ////////////////////////////////////////////////
  sbuf_id.resize (LocalNblock*Nball); 
  sbuf_idx.resize(LocalNblock*Nball); 
  ////////////////////////////////////////////////
  // Which buffer (or none) each block and dir maps to 
  ////////////////////////////////////////////////
  rbuf_id.resize (LocalNblock*Nball); 
  rbuf_idx.resize(LocalNblock*Nball); 
 
  //Build and Fill halo buffers
  // Enumerates in hop-0, hop-1, hop-2,...
  for(int mmu=0;mmu<Nball;mmu++){
      
    int mu = StencilReord[mmu];
    int delta[5];
    GetDelta(delta,mu,0);
    int deltasq = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]+delta[3]*delta[3]+delta[4]*delta[4];

    for(int idx=0;idx<LocalNblock;idx++){

      if ( deltasq < 5  ){

	int me,from,fromblock,to;
	int rev=0;
 
	NeighbourBlockId(idx,mu,me,from, fromblock,to,rev);
 
	rbuf_id [idx+LocalNblock*mmu] = -1; // Indicates local to node
	rbuf_idx[idx+LocalNblock*mmu] = fromblock*Nvec;
	
	int fidx=-1;
	if ( me != from ) {
	  for(int p=0;p<fproc.size();p++){
	    if (fproc[p]==from) fidx=p;
	  }
	  if(fidx==-1) { 
	    fproc.push_back(from);
	    fidx = fproc.size()-1;
	  }
	  rbuf_id [idx+LocalNblock*mmu] = fidx;
	  rbuf_idx[idx+LocalNblock*mmu] = recvbuf_bytes[fidx]/(2*sizeof(cFloat));
	  for(int i=0;i<Nvec;i++){
	    recvbuf_bytes[fidx]+=2*sizeof(cFloat);
	  }
	}
      
	int tidx=-1;
	if ( me != to ) {
	  for(int p=0;p<tproc.size();p++){
	    if (tproc[p]==to) tidx=p;
	  }
	  if(tidx==-1) {
	    std::vector<int> empty(0); // gather table for filling buffer
	    sbuf_gather.push_back(empty);
	    tproc.push_back(to);
	    tidx = tproc.size()-1;
	  }
	  sbuf_id [idx+LocalNblock*mmu] = tidx;
	  sbuf_idx[idx+LocalNblock*mmu] = sendbuf_bytes[tidx]/(2*sizeof(cFloat));

	  //////////////////////////////////////////////////////
	  // Information for this site and this direction in a flat table of everything
	  // to be gathered. Easier to thread.
	  //////////////////////////////////////////////////////
	  sbuf_gather_from_idx.push_back(fromblock);          //from index in array
	  sbuf_gather_buff_idx.push_back(sbuf_gather[tidx].size()); //buffer index in sendbuf
	  sbuf_gather_buff_id.push_back(tidx);                 //buffer id of send buffer

	  ///////////////////////////////
	  // Old format gather table
	  // Also maintains the size
	  ///////////////////////////////
	  sbuf_gather[tidx].push_back(fromblock);

	  for(int i=0;i<Nvec;i++){
	    sendbuf_bytes[tidx]+=2*sizeof(cFloat);
	  }      

	  ////////////////////////////////////////////////////////
	  // Track how much work we need to do for each depth
	  // Store depth0/1/2/3/4 sequentially in buffers
	  ////////////////////////////////////////////////////////
	  for(int d=0;d<5;d++){
	    if ( d >= deltasq ) {
	      sendbuf_depth_bytes[d][tidx]= sendbuf_bytes[tidx];
	      sbuf_gather_depth_count[d] = sbuf_gather_from_idx.size();
	    }
	  }

	}
      }
    } // idx
  } //mu


  if ( tproc.size() != fproc.size() ) {
    exit(-1);
  }

  HaloInitBuffers(this);
 
  // QMP init
  int Nmsg = tproc.size();
  xmit_msgmem.resize(Nmsg);
  recv_msgmem.resize(Nmsg);
  xmit_msghandle.resize(Nmsg);
  recv_msghandle.resize(Nmsg);

  for(int m=0;m<Nmsg;m++){
    xmit_msgmem[m] = QMP_declare_msgmem((void *) &sendbufs[m][0],sendbuf_bytes[m]) ;  // send buf
    recv_msgmem[m] = QMP_declare_msgmem((void *) &recvbufs[m][0],recvbuf_bytes[m]) ; // receive buf
    
    xmit_msghandle[m] = QMP_declare_send_to(xmit_msgmem[m],tproc[m],0);
    recv_msghandle[m] = QMP_declare_receive_from(recv_msgmem[m],fproc[m],0);
  }
  linop_d->BossLog("Initialised haloes\n");
  for(int d=0;d<5;d++){
    int bytes = 0;
    for(int m=0;m<Nmsg;m++) {
      bytes = bytes + sendbuf_depth_bytes[d][m];
    }
    linop_d->BossLog("Depth %d : %d gathers %d bytes\n",d,sbuf_gather_depth_count[d],bytes);
  }

}
template<class cFloat> 
std::complex<cFloat> * 
BfmHDCG<cFloat>::HaloGetStencilData(int _idx,int _ballidx,
				 std::vector<std::complex<cFloat> > &my_data)
{
  // Locate the neighbour
  int me, from, to;
  int fromblock;
  int mu;

  std::complex<cFloat> *nbr_data;
  mu =_ballidx;

  int bid = rbuf_id [_idx+LocalNblock*mu] ;
  int bidx= rbuf_idx [_idx+LocalNblock*mu] ;

  if (bid==-1){ 
    nbr_data = &my_data[bidx];
  } else { 
    nbr_data = &recvbufs[bid][bidx];
  }
  return nbr_data;
}

//
// Optimisation thoughts
// i)  Overlap commms and compute here. 2x in ldop
// ii) Drop fill in; add an infra-red shift use as smoother for 2nd level -> Adef2
//     Possible 4x
//      
// No particular saving possible for MdagM smoother from fill in dropping
// as must evaluate Moe Meo and its adjoint.
//

template<class cFloat> 
void BfmHDCG<cFloat>::HaloExchange(std::vector<std::complex<cFloat> >&my_data,int depth)
{
  int Nmsg = tproc.size();

  int me, thrlen, throff;
  static int pooh;
  static unsigned long totbytes=0;

  linop_d->thread_barrier();

  uint64_t t0 = GetTimeBase();
#if 1 
  //  int nwork = sbuf_gather_from_idx.size();
  int nwork = sbuf_gather_depth_count[depth];
  linop_d->thread_work(nwork,me,thrlen,throff);
  for(int w=throff;w<throff+thrlen;w++){
    int tidx = sbuf_gather_buff_id[w];
    int fidx = sbuf_gather_from_idx[w];
    int bidx = sbuf_gather_buff_idx[w];
    for(int i=0;i<Nvec;i++){
      sendbufs[tidx][bidx*Nvec+i] = my_data[i+fidx*Nvec];
    }
  }
#else 
  int nwork = sbuf_gather.size();
  linop_d->thread_work(nwork,me,thrlen,throff);
  for(int tidx=throff;tidx<throff+thrlen;tidx++){
#if 0
    for(int j=0;j<sbuf_gather[tidx].size();j++){
      for(int i=0;i<Nvec;i++){
	sendbufs[tidx][j*Nvec+i] = my_data[i+sbuf_gather[tidx][j]*Nvec];
      }
    }
#else
    if ( sizeof(cFloat) == sizeof(double) ) { 
      qpx_gather(Nvec,
		 &sbuf_gather[tidx][0],
		 sbuf_gather[tidx].size(),
		 (double *)&sendbufs[tidx][0],
		 (double *)&my_data[0]);
    } else { 
      qpx_gather_f(Nvec,
		 &sbuf_gather[tidx][0],
		 sbuf_gather[tidx].size(),
		 (float *)&sendbufs[tidx][0],
		 (float *)&my_data[0]);
    }
#endif
  }
#endif
  linop_d->thread_barrier();
  uint64_t t1 = GetTimeBase();
  HaloExchangeCommStart(this,depth);
  linop_d->thread_barrier();
  uint64_t tw = GetTimeBase();
  HaloExchangeCommComplete(this,depth);

#ifdef DEBUG_HDCG
  linop_d->thread_barrier();
  uint64_t t2 = GetTimeBase();
  if ( MGreport ) {
    linop_d->ThreadBossMessage("Halo: Gather %ld cyc Start %ld Wait %ld total %ld\n",t1-t0,tw-t1,t2-tw,t2-t0);
  }
  linop_d->thread_barrier();
#endif
}
#ifndef BFM_BGQ_SPI_DSLASH
  ////////////////////////////////////////////////
  // QMP implementation
  ////////////////////////////////////////////////
template<class cFloat> void HaloExchangeCommComplete(BfmHDCG<cFloat> * BMG,int depth)
{
  int Nmsg=BMG->tproc.size();
  int me = BMG->linop_d->thread_barrier();
  if ( me == 0){
    for(int m=0;m<Nmsg;m++) QMP_wait (BMG->xmit_msghandle[m]);
    for(int m=0;m<Nmsg;m++) QMP_wait (BMG->recv_msghandle[m]);
  }
  BMG->linop_d->thread_barrier();
}
template<class cFloat> void HaloExchangeCommStart(BfmHDCG<cFloat> * BMG,int depth)
{
  int Nmsg=BMG->tproc.size();
  int me = BMG->linop_d->thread_barrier();
  if ( me == 0){
    for(int m=0;m<Nmsg;m++) QMP_start(BMG->xmit_msghandle[m]);
    for(int m=0;m<Nmsg;m++) QMP_start(BMG->recv_msghandle[m]);
  }
  BMG->linop_d->thread_barrier();
}
template<class cFloat> void HaloInitBuffers(BfmHDCG<cFloat> * BMG)
{
  int Nmsg = BMG->tproc.size();

  for(int m=0;m<Nmsg;m++){
    BMG->sendbufs[m] = (std::complex<cFloat> *)bfm_alloc(BMG->sendbuf_bytes[m]);
    BMG->recvbufs[m] = (std::complex<cFloat> *)bfm_alloc(BMG->recvbuf_bytes[m]);
  }

  /* QMP 
  // QMP init
  BMG->xmit_msgmem.resize(Nmsg);
  BMG->recv_msgmem.resize(Nmsg);
  BMG->xmit_msghandle.resize(Nmsg);
  BMG->recv_msghandle.resize(Nmsg);

  for(int m=0;m<Nmsg;m++){
    BMG->xmit_msgmem[m] = QMP_declare_msgmem((void *) &BMG->sendbufs[m][0],BMG->sendbuf_bytes[m]) ;  // send buf
    BMG->recv_msgmem[m] = QMP_declare_msgmem((void *) &BMG->recvbufs[m][0],BMG->recvbuf_bytes[m]);
    
    BMG->xmit_msghandle[m] = QMP_declare_send_to(BMG->xmit_msgmem[m],BMG->tproc[m],0);
    BMG->recv_msghandle[m] = QMP_declare_receive_from(BMG->recv_msgmem[m],BMG->fproc[m],0);
  }
  */


  for(int m=0;m<Nmsg;m++){
    BGNET_SetSendBuffer((void *) &BMG->sendbufs[m][0], m+1 ,BMG->sendbuf_bytes[m]);
    BGNET_SetRecvBuffer((void *) &BMG->recvbufs[m][0], m+1, BMG->recvbuf_bytes[m]);
  }



  BMG->linop_d->BossLog("Initialised haloes\n ");
}
#endif



template<class cFloat> 
void BfmHDCG<cFloat>::OrthogonaliseSubspace(void)
{
  // Remove earlier (normalised) vectors
  std::vector<std::complex<double> > bip(LocalNblock);
  std::vector<double> bn(LocalNblock);
  int Nv=0;

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {

      for(int v=0;v<Nvec;v++){
	if ( linop_d->norm(subspace_d[v]) > 1.0e-12 ) { 
	  if ( t==0 ) Nv++;
	}
      }
    }
  }

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      for(int v=0;v<Nv;v++){
	for(int u=0;u<v;u++){
	  //Inner product & remove component
	  linop_d->block_inner(subspace_d[u],subspace_d[v],LocalNblock,(double *)&bip[0],(double *)&inner_reduce[0],&local_block_id[Odd][0]);
	  linop_d->block_zaxpy(subspace_d[v],subspace_d[u],subspace_d[v],-1,(double *)&bip[0],&local_block_id[Odd][0]);
	}
	// Normalise this vector
	linop_d->block_norm(subspace_d[v],LocalNblock,(double *)&bn[0],&local_block_id[Odd][0]);
	linop_d->block_normalise(subspace_d[v],LocalNblock,(double *)&bn[0],&local_block_id[Odd][0]);
      }
    }
  }


#if DEBUG_HDCG
  for(int u=0;u<Nv;u++){
    int v=u; 
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      linop_d->block_inner(subspace_d[u],subspace_d[v],LocalNblock,(double *)&bip[0],(double *)&inner_reduce[0],&local_block_id[Odd][0]);
    }
  }
      double c = (u==v);
      for(int i=0;i<bip.size();i++) { 
	if ( abs(bip[i]-c) > 1.0e-6 ) {  
	  cout << "Block "<< i<<" vec "<< u<< "," << v << " inner= " << bip[i] << " : " 
	       <<"Doesnt look orthogonal and should be "<<c<<endl;
	  exit(0);
	}
      }
  }

  for(int v=0;v<Nv;v++){
  for(int u=0;u<Nv;u++){
    
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      linop_d->block_inner(subspace_d[u],subspace_d[v],LocalNblock,(double *)&bip[0],(double *)&inner_reduce[0],&local_block_id[Odd][0]);
    }
  }
      double c = (u==v);
      for(int i=0;i<bip.size();i++) { 
	if ( abs(bip[i]-c) > 1.0e-6 ) {  
	  cout << "Block "<< i<<" vec "<< u<< "," << v << " inner= " << bip[i] << " : " 
	       <<"Doesnt look orthogonal and should be "<<c<<endl;
	  exit(0);
	}
      }
  }}
#endif
}

template<class cFloat> 
template<class lFloat>
void BfmHDCG<cFloat>::ComputeLittleMatrixColored(bfm_internal<lFloat> * lop, Fermion_t *subspace)
{
  if ( lop->SPIcomms() ) lop->comm_init();

  Fermion_t zero_t= lop->allocFermion();
  Fermion_t phi_t = lop->allocFermion();
  Fermion_t tmp_t = lop->allocFermion();
  Fermion_t mmp   = lop->allocFermion();
  Fermion_t mp    = lop->allocFermion();

  uint64_t t,t1,t2,t_mat, t_proj;
  t_mat = t_proj=0;

  t1=GetTimeBase();

  lop->BossLog("HDCG: ComputeColored: ");

#pragma omp parallel 
  {
#pragma omp for 
  for(int i=0;i<lop->threads;i++) {

  ////////////////////////////////////
  // Preallocate maximum sized matrix for FT.
  ////////////////////////////////////
  std::vector< std::vector< std::complex<double> > > imat; // Inverse matrix for this mu.
  imat.resize(3);
  for(int i=0;i<3;i++) imat[i].resize(3);
    
  int me, thrlen,throff;
  me = lop->thread_barrier();
  lop->fill(zero_t,0.0);

  for(int i=0;i<Nvec;i++){

    if ( lop->isBoss() && !me ) printf("."); fflush(stdout);
	
    // Problematic for threading
    for(int b=0;b<Nball;b++){  // Loop over momenta (Nball)

      /////////////////////////////////////////////////////
      // Stick a different phase on every block
      /////////////////////////////////////////////////////
      int lb[5];// Local Block
      int gb[5];// Global Block
      int    imom[5];
      double dmom[5];

      GetDelta(imom,b,0);  
      for(int mu=0;mu<5;mu++){
 	dmom[mu] = imom[mu]*2*M_PI/global_nb[mu];
      }
      
      lop->thread_barrier();
      if ( !me ) {
      for(lb[0]=0;lb[0]<local_nb[0];lb[0]++){ // Redundant computation alert
      for(lb[1]=0;lb[1]<local_nb[1];lb[1]++){
      for(lb[2]=0;lb[2]<local_nb[2];lb[2]++){
      for(lb[3]=0;lb[3]<local_nb[3];lb[3]++){
      for(lb[4]=0;lb[4]<local_nb[4];lb[4]++){

	double pdotx = 0.0;
	for(int mu=0;mu<5;mu++) gb[mu] = ncoor[mu]*local_nb[mu]+lb[mu];
	for(int mu=0;mu<5;mu++) pdotx+= gb[mu]*dmom[mu];

	std::complex<double> pha(cos(pdotx),sin(pdotx));

	int lbidx=GetLocalBlockIdx(lb);
	phases[lbidx] = pha;

      }}}}}
      }
      lop->thread_barrier();


      ///////////////////////////////////////////////////////
      // We apply the matrix Nball times using these phases
      ///////////////////////////////////////////////////////
      lop->block_zaxpy(phi_t,subspace[i],zero_t,1.0,(double *)&phases[0],&local_block_id[Odd][0]);
      t=GetTimeBase();
      lop->Mprec(phi_t,mp,tmp_t,0);
      lop->Mprec(mp,mmp,tmp_t,1); 
      if( me==0) t_mat+=GetTimeBase()-t;

      double mmp_norm= lop->norm(mmp);
      double phi_norm= lop->norm(phi_t);
      double ss_norm = lop->norm(subspace[i]);

      ///////////////////////////////////////////////////////
      // Project this into subspace
      ///////////////////////////////////////////////////////
      t=GetTimeBase();
      
      int smin = 0;
      int smax = lop->Ls-1;
      lop->template block_project<cFloat>(subspace,Nvec,mmp,LocalNblock,
			     (cFloat *)&ComputeProj[b][0],
  			     (double *)&inner_reduce[0],
			     &local_block_id[Odd][0],smin,smax);
      if( me==0)   t_proj+=GetTimeBase()-t;

      lop->thread_work(LocalNblock,me,thrlen,throff);
      for(int lbidx=throff;lbidx<throff+thrlen;lbidx++) { 
	for(int j=0;j<Nvec;j++){
	  std::complex<cFloat> pha(real(phases[lbidx]),imag(phases[lbidx]));
	  ComputeProj[b][lbidx*Nvec+j] = ComputeProj[b][lbidx*Nvec+j]*conj(pha);
	}
      }
      lop->thread_barrier();
    }


    //////////////////////////////////////////////////////////
    // Solve system of equations to get Aij
    //////////////////////////////////////////////////////////
    /*
     *     Here, k,l index which possible shift within the 3^Nd "ball" connected by MdagM.
     *
     *     conj(phases[block]) proj[k][ block*Nvec+j ] =  \sum_ball  e^{i q_k . delta} < phi_{block,j} | MdagM | phi_{(block+delta),i} > 
     *                                                 =  \sum_ball e^{iqk.delta} A_ji
     *
     *     Must invert matrix M_k,l = e^[i q_k . delta_l]
     *
     *     Where q_k = delta_k . (2*M_PI/global_nb[mu])
     */
    {
      
      lop->thread_work(LocalNblock,me,thrlen,throff);
      for(int lbidx=throff;lbidx<throff+thrlen;lbidx++){
	for(int j=0;j<Nvec;j++){

	for(int mu=0;mu<5;mu++){

	  double pmu=2*M_PI/global_nb[mu];      
	  std::complex<double> pha(cos(pmu),-sin(pmu));
	  std::vector<std::complex<double> > FT(Nball);
	  std::complex<double> a = pha+1.0;
	  std::complex<double> b = (pha-1.0)*(pha-1.0);
	  std::complex<double> ab=a*b;

	  ///////////////////////////////////////////////////////////////////////////////
	  // Stencil == 3
	  //  pha* 1 pha                 pha/a -pha        pha^2/a
	  //  1    1  1     -> inv ->    -pha  (1+pha^2)   -pha           / [ 1-pha]^2     where a=1+pha
 	  //  pha  1 pha*               pha^2/a -pha       pha/a
	  ///////////////////////////////////////////////////////////////////////////////
	  ///////////////////////////////////////////////////////////////////////////////
	  // Stencil == 2
	  // -1 1                  1/2   -1 1 
	  //  1 1      -> inv ->          1 1
 	  //                  
	  ///////////////////////////////////////////////////////////////////////////////
	  if ( stencil_size[mu]==3 ) {
	    imat[0][0] = pha/ab; 	   imat[0][1] =-pha/b;	                  imat[0][2] = pha*pha/ab;
	    imat[1][0] =-pha/b;	           imat[1][1] = (1.0+pha*pha)/b;	  imat[1][2] =-pha/b;
	    imat[2][0] = pha*pha/ab;       imat[2][1] =-pha/b;         	          imat[2][2] = pha/ab;
	  } else if (stencil_size[mu]==2) { 
	    imat[0][0] = -0.5; 	   imat[0][1] = 0.5;
	    imat[1][0] =  0.5;	   imat[1][1] = 0.5;
	  } else if (stencil_size[mu]==1) { 
	    imat[0][0]=1.0;
	  }

	  int d[5];
	  int nd[5];

	  for(d[0]=0;d[0]<stencil_size[0];d[0]++){
	  for(d[1]=0;d[1]<stencil_size[1];d[1]++){
	  for(d[2]=0;d[2]<stencil_size[2];d[2]++){
	  for(d[3]=0;d[3]<stencil_size[3];d[3]++){
	  for(d[4]=0;d[4]<stencil_size[4];d[4]++){

	    int  d_mu;
	    d_mu =   d[4]+stencil_size[4]*(d[3]+stencil_size[3]*(d[2]+stencil_size[2]*(d[1]+stencil_size[1]*d[0])));
	    for(int i=0;i<5;i++) nd[i]=d[i];
	    FT[d_mu] = 0;
	    
	    for(nd[mu]=0;nd[mu]<stencil_size[mu];nd[mu]++){
	      int nd_mu;
	      nd_mu = nd[4]+stencil_size[4]*(nd[3]+stencil_size[3]*(nd[2]+stencil_size[2]*(nd[1]+stencil_size[1]*nd[0])));
	      std::complex<double> CP(real(ComputeProj[nd_mu][lbidx*Nvec+j]),
				      imag(ComputeProj[nd_mu][lbidx*Nvec+j]));
	      FT[d_mu]+= imat[d[mu]][nd[mu]] * CP;
	      
	    }
	  }}}}}

	  for(int b=0;b<Nball;b++) ComputeProj[b][lbidx*Nvec+j] = FT[b];

	}//5-Directions FT
	////////////////////////
	// Copy back solution of system of eqns
	////////////////////////
	for(int b=0;b<Nball;b++)  Asparse[lbidx*Nball*Nvec*Nvec+b*Nvec*Nvec+Nvec*j+i] = ComputeProj[b][lbidx*Nvec+j];
	for(int b=0;b<Nball;b++)  AsparseSingle[lbidx*Nball*Nvec*Nvec+b*Nvec*Nvec+Nvec*j+i] = ComputeProj[b][lbidx*Nvec+j];
	
	}//j
      }// Block
    }// Scope
  }//Nvec
  }}

  if ( lop->isBoss() ) printf("\n"); fflush(stdout);
  int diag_ball_idx = StencilReord[0];
  if ( lop->isBoss() ) printf("Diagonal ball index is %d\n",diag_ball_idx); fflush(stdout);

  ///////////////////////////////////////////////////////////////
  // Not sure if lapack is reentrant so don't thread this.
  ///////////////////////////////////////////////////////////////
  for(int b=0;b<LocalNblock;b++){
    for(int j=0;j<Nvec*Nvec;j++){
      AsparseDiagInv[b*Nvec*Nvec+j] = Asparse[b*Nvec*Nvec*Nball +diag_ball_idx*Nvec*Nvec +j] ;
    }    
    if ( sizeof(cFloat) == sizeof(double) ) { 
      LapackHermitianInvert(Nvec,(double *)&AsparseDiagInv[b*Nvec*Nvec]);
    } else { 
      LapackHermitianInvert_f(Nvec,(float *)&AsparseDiagInv[b*Nvec*Nvec]);
    }
  }

  ////////////////////////////////////////////////////
  // Some elements are zero by symmetry. Drop these
  ////////////////////////////////////////////////////
  for(int b=0;b<Nball;b++){
    double nrm =0.0;
    for(int j=0;j<Nvec*Nvec;j++){
      std::complex<double> e = Asparse[b*Nvec*Nvec+j];
      nrm=nrm+real(e)*real(e)+imag(e)*imag(e);
    }
    if(nrm < 1.0e-7 ) {
#ifdef DEBUG_HDCG
      //      lop->BossLog("Stencil %d is zero %le distance %d\n",b,nrm,StencilDepth[b]);
#endif
      StencilNonZero[b] = 0;
    } else { 
#ifdef DEBUG_HDCG
      //      lop->BossLog("Stencil %d is nonzero %le distance %d \n",b,nrm,StencilDepth[b]);
#endif
      StencilNonZero[b] = 1;
    }
  }
  
  t2=GetTimeBase();
  lop->BossLog("ComputeColored<%d>: Total    : %ld ms\n",sizeof(lFloat),(t2-t1 )/MHz()/1000);
  lop->BossLog("ComputeColored<%d>: MdagM    : %ld ms\n",sizeof(lFloat),(t_mat )/MHz()/1000);
  lop->BossLog("ComputeColored<%d>: Proj     : %ld ms\n",sizeof(lFloat),(t_proj)/MHz()/1000);
  
  lop->freeFermion( zero_t ); 
  lop->freeFermion( phi_t ); 
  lop->freeFermion( tmp_t ); 
  lop->freeFermion( mmp   ); 
  lop->freeFermion( mp    );
  
}

template void BfmHDCG<float>::ComputeLittleMatrixColored(bfm_internal<float> * lop, Fermion_t *subspace);
template void BfmHDCG<float>::ComputeLittleMatrixColored(bfm_internal<double> * lop, Fermion_t *subspace);
template void BfmHDCG<double>::ComputeLittleMatrixColored(bfm_internal<float> * lop, Fermion_t *subspace);
template void BfmHDCG<double>::ComputeLittleMatrixColored(bfm_internal<double> * lop, Fermion_t *subspace);

template<class cFloat> 
void BfmHDCG<cFloat>::PromoteFromSubspace_f(std::vector<std::complex<cFloat> > &v,Fermion_t prom)
{
  linop_f->block_promote<cFloat>(prom,subspace_f,Nvec,(cFloat *)&v[0],&local_block_id[Odd][0],s_min,s_max);
}

// Threaded
template<class cFloat> 
void BfmHDCG<cFloat>::PromoteFromSubspace(std::vector<std::complex<cFloat> > &v,Fermion_t prom)
{
  linop_d->block_promote<cFloat>(prom,subspace_d,Nvec,
				 (cFloat *)&v[0],
				 &local_block_id[Odd][0],s_min,s_max);
}

template<class cFloat> 
void BfmHDCG<cFloat>::ProjectToSubspace(Fermion_t vec, std::vector<std::complex<cFloat> > &proj)
{
  linop_d->block_project<cFloat>(subspace_d,Nvec,vec,LocalNblock,
			     (cFloat *)&proj[0],
  			     (double *)&inner_reduce[0],
			     &local_block_id[Odd][0],s_min,s_max);
}
template<class cFloat> 
void BfmHDCG<cFloat>::ProjectToSubspace_f(Fermion_t vec, std::vector<std::complex<cFloat> > &proj)
{
  linop_f->block_project<cFloat>(subspace_f,Nvec,vec,LocalNblock,
			     (cFloat *)&proj[0],
  			     (double *)&inner_reduce[0],
			     &local_block_id[Odd][0],s_min,s_max);
}






template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverse(std::vector<std::complex<cFloat> > &v,std::vector<std::complex<cFloat> > &vinv)
{
  if ( LittleDopSolver == LittleDopSolverDeflCG ) {
    ApplyInverseDeflCG(v,vinv);
    return;
  }

  if ( LittleDopSolver == LittleDopSolverADef2 ) {
    ApplyInverseADef2(v,vinv);
    return;
  }

  if ( LittleDopSolver == LittleDopSolverMCR ) {
    ApplyInverseMCR(v,vinv);
    return;
  }  
  if ( LittleDopSolver == LittleDopSolverCG ) {
    ApplyInverseCG(v,vinv);
    return;
  }
  exit(0);
}

template<class cFloat> void BfmHDCG<cFloat>::LdopM1(std::vector<std::complex<cFloat> > &in,
						    std::vector<std::complex<cFloat> > &out,
						    Chebyshev & Approx,
						    LdopM1ControlType control)
{
  //Adef1
  //MP  + Q
  // [MP+Q] = M  [ 1 - A Q] + Q 
  //
  // Qin
  // 1-A Q in from precompute AQ.
  // PolyMdagMprec -- eliminates the extra mult. 20% speed up of ldop.
  // Further 
  static int call;
  int me=linop_d->thread_barrier();

  static std::vector<double> CGpoly;

#define LDOP_ADEF1
#ifdef  LDOP_ADEF1
  uint64_t t0 = GetTimeBase();
  LdopDeflationProject(in,DeflKrylovProj);
  uint64_t t1 = GetTimeBase();
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  uint64_t t2 = GetTimeBase();
  LdopDeflationMatPromote(DeflKrylovMss,Krylov_Amu);
  uint64_t t3 = GetTimeBase();
  LdopDeflationPromote(DeflKrylovMss,out);
  axpy(Krylov_Amu,Krylov_Amu,in,-1.0);
  uint64_t t4 = GetTimeBase();

  // 0=LdopM1Chebyshev -- Chebyshev
  // 1 -- Mirs CG (record poly)
  // 2 -- MirsPoly (replay poly)
  // 3 -- Zero hop exact invers preconditioner
  if (control==LdopM1Chebyshev) {        
    PolyMdagMprec(Krylov_Amu,Krylov_tmp,Approx);
  } else if (control==LdopM1Mirs || control==LdopM1MirsPolyRecord) {
    int halo=1;
    double lo=LdopM1Lo;
    double rr=LdopM1resid;
    int maxit=LdopM1iter;
    ApplyInverseCGopt(Krylov_Amu,Krylov_tmp,CGpoly,lo,rr,maxit,halo);
  } else if ( control==LdopM1MirsPoly ) {
    PolyMdagMprec(Krylov_Amu,Krylov_tmp,CGpoly);
  } else if ( control==LdopM1DiagInv) {
    ApplyDiagInv(Krylov_Amu,Krylov_tmp);
  }
  axpy(out,Krylov_tmp,out,1.0);
  uint64_t t5 = GetTimeBase();
#else
  //Adef2
  // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]

  // M is a Chebyshev polynomial approximate inverse

  uint64_t t0 = GetTimeBase();
  PolyMdagMprec(in,Krylov_tmp,Approx);
  uint64_t t1 = GetTimeBase();
  ApplyThread(Krylov_tmp,Krylov_Amu,1);
  uint64_t t2 = GetTimeBase();
  axpy(Krylov_Amu,Krylov_Amu,in,-1.0);          // Atmp  = in - A Min

  LdopDeflationProject(Krylov_Amu,DeflKrylovProj);
  uint64_t t3 = GetTimeBase();
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  uint64_t t4 = GetTimeBase();
  LdopDeflationPromote(DeflKrylovMss,Krylov_Amu);

  axpy(out,Krylov_tmp,Krylov_Amu,1.0); // Min+tmp

  uint64_t t5 = GetTimeBase();
#endif 

#ifdef DEBUG_HDCG
  linop_d->thread_barrier();
  if ( call == 17 ) { 
    linop_d->ThreadBossLog("LdopM1 : %d cyc\n",t1-t0);
    linop_d->ThreadBossLog("LdopM1 : %d cyc\n",t2-t1);
    linop_d->ThreadBossLog("LdopM1 : %d cyc\n",t3-t2);
    linop_d->ThreadBossLog("LdopM1 : %d cyc\n",t4-t3);
    linop_d->ThreadBossLog("LdopM1 : %d cyc\n",t5-t4);
    linop_d->ThreadBossLog("LdopM1 : total %d cyc\n",t5-t0);
  }
  linop_d->thread_barrier();
  if ( !me ) call++;
  linop_d->thread_barrier();
#endif

}

template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverseADef2(std::vector<std::complex<cFloat> > &src,
					std::vector<std::complex<cFloat> > &x, double shift )
{
  double rtzp,rtz,a,d,b;
  uint64_t t0,t1;
  uint64_t mat0,mat1;
  uint64_t prec0,prec1;

  int me = linop_d->thread_barrier();
  Chebyshev Approx;  
  Approx.Init(LdopM1Lo,LdopM1Hi,LdopM1iter,PolynomialShape);


  double tmp =  norm_vec(src); // Forces sync of threads & nodes
  t0 = GetTimeBase();

  // Deflated precond CG (3.6 in Saad umsi-98-97)
  ///////////////////////////////////
  // x_{-1} = src
  ///////////////////////////////////

  axpy(x,src,src,0.0);
  ApplyThread(x,Krylov_Ap);
  axpy(Krylov_r, Krylov_Ap, src,-1.0);        // r_{-1} = src - A x

  ///////////////////////////////////
  // Choose x_0 such that 
  // x_0 = guess +  (A_ss^inv) r_s
  ///////////////////////////////////
  LdopDeflationProject(Krylov_r,DeflKrylovProj);
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  LdopDeflationPromote(DeflKrylovMss,Krylov_Ap);
  axpy(x,x,Krylov_Ap,1.0);

  //////////////////////////////////
  // Recomputes r=src-x0
  //////////////////////////////////
  ApplyThread(x,Krylov_Ap);
  axpy (Krylov_r,Krylov_Ap, src,-1.0);  

  LdopDeflationProject(Krylov_r,DeflKrylovProj);

  double nv = norm_vec(DeflKrylovProj);
  
  ///////////////////////////////////////
  // Compute z = M1 x
  ///////////////////////////////////////
  LdopM1ControlType control = LdopM1control;
  if ( (LdopM1control == LdopM1MirsPoly) ) { // rerecord the polynomial
    control = LdopM1MirsPolyRecord;
  }

  LdopM1(Krylov_r,Krylov_mu,Approx,control);
  rtzp =innerProductReal(Krylov_r,Krylov_mu);

  // M2 trivial for all except Def2
  // M2 p = z
  axpy(Krylov_p,Krylov_mu,Krylov_mu,0.0);

  double ssq =  norm_vec(src); // Forces sync of threads & nodes
  double rsq =  LittleDopSolverResidual*LittleDopSolverResidual*ssq;
  int singleprecision = 0;
  if ( LittleDopSolverResidual > 1.0e-6 ) singleprecision=1;

  for (int k=1;k<=10000;k++){

    rtz=rtzp;
    double pn = norm_vec(Krylov_p);

    MGreport=0;
    //    if(linop_d->isBoss() && (k==3) ) MGreport=1;

    mat0=GetTimeBase();
    ApplyThreadOpt(Krylov_p,Krylov_Ap,singleprecision);// Full depth halo
    mat1=GetTimeBase();
    double Apn = norm_vec(Krylov_Ap);
    d  = innerProductReal(Krylov_p,Krylov_Ap);
    a = rtz/d;


    axpy(x,Krylov_p,x,a);
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a);
    
    double rn =norm_vec(Krylov_r);

    prec0=GetTimeBase();

    control=LdopM1control;
    if ( (LdopM1control == LdopM1MirsPoly) && ( (k%5) == 1 ) ) { // rerecord the polynomial
	control = LdopM1MirsPolyRecord;
    }
    LdopM1(Krylov_r,Krylov_mu,Approx,control);
    prec1=GetTimeBase();
    rtzp =innerProductReal(Krylov_r,Krylov_mu);
    b = rtzp/rtz;

    axpy(Krylov_p,Krylov_p,Krylov_mu,b);     // Overlap with comms if lift
#ifdef DEBUG_HDCG
    linop_d->ThreadBossDebug("ApplyInverseADef2: k= %d residual = %le bad %le %le %le pnApn %le %le \n",k,sqrt(rn/ssq),b,a,d,pn,Apn);
#endif

    // Stopping condition
    if ( rn <= rsq ) { 
      t1 = GetTimeBase();

      /*
      ApplyThread(x,Krylov_Ap);
      axpy(Krylov_r,src,Krylov_Ap,-1.0);
      double tmpnorm=sqrt(norm_vec(Krylov_r));
      double srcnorm=sqrt(ssq);
      double true_residual = tmpnorm/srcnorm;

      linop_d->ThreadBossLog("ApplyInverseAdef2<%d>: true residual is %le, %d iterations %f millisecs\n",
			     sizeof(cFloat),
			     true_residual,k,
			     (t1-t0)/MHz()/1000.0);
      */

      linop_d->ThreadBossLog("ApplyInverseAdef2<%d>: Mat %f Prec %f millisecs %d iters cp=%le\n",
			     sizeof(cFloat),
			     k*(mat1-mat0)/MHz()/1000.0,
			     k*(prec1-prec0)/MHz()/1000.0,
			     k,
			     sqrt(rn/rsq)
			     );

      Approx.End();


      return;
    }

  }
  linop_d->Error("ApplyInverseADef2: CG not converged \n"); fflush(stdout);
  exit(0);

}

template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverseDeflCG(std::vector<std::complex<cFloat> > &src,
					      std::vector<std::complex<cFloat> > &x, double shift )
{
  double rtzp,rtz,a,d,b;

  int me = linop_d->thread_barrier();

  double ssq =  norm_vec(src); // Forces sync of threads & nodes

  // Deflated precond CG (3.6 in Saad umsi-98-97)
  ///////////////////////////////////
  // x_{-1} = src
  ///////////////////////////////////

  axpy(x,src,src,0.0);
  ApplyThread(x,Krylov_Ap);
  axpy(Krylov_r, Krylov_Ap, src,-1.0);        // r_{-1} = src - A x

  ///////////////////////////////////
  // Choose x_0 such that 
  // x_0 = guess +  (A_ss^inv) r_s
  // 
  // W^T (src - A x_0) = src_s - A guess_s - r_s  
  //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
  //                   = 0 
  ///////////////////////////////////

  LdopDeflationProject(Krylov_r,DeflKrylovProj);
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  LdopDeflationPromote(DeflKrylovMss,Krylov_Ap);
  axpy(x,x,Krylov_Ap,1.0);

  // Recomputes r=src-x0
  ApplyThread(x,Krylov_Ap);
  axpy (Krylov_r,Krylov_Ap, src,-1.0);  

  rtzp =innerProductReal(Krylov_r,Krylov_r);

  // Check orthogonal
  LdopDeflationProject(Krylov_r,DeflKrylovProj);     
  double nv=norm_vec(DeflKrylovProj);
  
  ///////////////////////////////////////
  // Solve for Mss mu = P A r and set p = r-mu
  ///////////////////////////////////////
  ApplyThread(Krylov_r,Krylov_Ap);
  LdopDeflationProject(Krylov_Ap,DeflKrylovProj);
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  LdopDeflationPromote   (DeflKrylovMss,Krylov_mu);
  LdopDeflationMatPromote(DeflKrylovMss,Krylov_Amu);

  axpy (Krylov_p,Krylov_mu,Krylov_r,-1.0);
  axpy (Krylov_Ap,Krylov_Amu,Krylov_Ap,-1.0);

  double rsq =  LittleDopSolverResidual*LittleDopSolverResidual*ssq;
  int singleprecision = 0;
  if ( LittleDopSolverResidual > 1.0e-6 ) singleprecision=1;

  uint64_t t0= GetTimeBase();

  for (int k=1;k<=10000;k++){

    uint64_t mat_0= GetTimeBase();

    MGreport=0;
    //    if(linop_d->isBoss() && (k==100) ) MGreport=1;

    rtz=rtzp;
    d  = innerProductReal(Krylov_p,Krylov_Ap);
    a = rtz/d;
    axpy(x,Krylov_p,x,a);
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a);

    rtzp =norm_vec(Krylov_r);
    b = rtzp/rtz;

    // WAW mu = W A z
    uint64_t mat_2= GetTimeBase();
    // Could overlap with norm_vec, axpy but that is NOT much overhead
    ApplyThread(Krylov_r,Krylov_Ar,singleprecision);
    uint64_t mat_3= GetTimeBase();

    uint64_t defl_0= GetTimeBase();
    LdopDeflationProject(Krylov_Ar,DeflKrylovProj);
    uint64_t defl_1= GetTimeBase();
    LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
    uint64_t defl_2= GetTimeBase();
    LdopDeflationPromote   (DeflKrylovMss,Krylov_mu);
    LdopDeflationMatPromote(DeflKrylovMss,Krylov_Amu);
    uint64_t defl_3= GetTimeBase();


    //=>    Krylov_p = b*Krylov_p + Krylov_r -  Krylov_mu
    axpy(Krylov_p,Krylov_p,Krylov_r,b);     // Overlap with comms if lift
    axpy(Krylov_p,Krylov_mu,Krylov_p,-1.0);

    //=>    Krylov_Ap = bKrylov_p + Krylov_Ar +  Krylov_Amu
    axpy(Krylov_Ap,Krylov_Ap,Krylov_Ar,b);
    axpy(Krylov_Ap,Krylov_Amu,Krylov_Ap,-1.0);

    if ( k%200 == 0 ){
      linop_d->ThreadBossDebug("ApplyInverseDelfCG: k= %d residual = %le \n",k,sqrt(rtzp/ssq));
    }

    // Stopping condition
    if ( rtzp <= rsq ) { 
      uint64_t t1= GetTimeBase();
      ApplyThread(x,Krylov_Ap);
      axpy(Krylov_r,src,Krylov_Ap,-1.0);
      double tmpnorm=sqrt(norm_vec(Krylov_r));
      double srcnorm=sqrt(ssq);
      double true_residual = tmpnorm/srcnorm;

      linop_d->ThreadBossLog("ApplyInverseDeflCG<%d>: true residual is %le after %d iterations %f milliseconds \n",sizeof(cFloat),
			     true_residual,k,(t1-t0)/MHz()/1000.0);
      if (0){
	//	linop_d->ThreadBossMessage("*** Last Iter Stats: *** \n");
	linop_d->ThreadBossMessage("Cycles total %ld\n",t1-mat_0);
	linop_d->ThreadBossMessage("Cycles mat %ld precision %d\n",mat_3-mat_2,singleprecision);
	linop_d->ThreadBossMessage("Cycles defl1 %ld\n",defl_1-defl_0);
	linop_d->ThreadBossMessage("Cycles defl2 %ld\n",defl_2-defl_1);
	linop_d->ThreadBossMessage("Cycles defl3 %ld\n",defl_3-defl_2);
      }


      return;
    }

  }
  linop_d->Error("ApplyInverseDeflCG: CG not converged \n"); fflush(stdout);
  exit(0);
}

template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverseMCR(std::vector<std::complex<cFloat> > &v,std::vector<std::complex<cFloat> > &vinv)
{
  int me, thrlen,throff;
  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;

  double rAr;
  double rArp;
  double pAAp;

  // Conjugate gradient on A. A is hermitian posdef.
  // Need these to be shared across threads
  me = linop_d->thread_barrier();

  linop_d->thread_work(LocalNsubspace,me,thrlen,throff);
  for(int i=throff;i<throff+thrlen;i++){
    vinv[i]=0;
  }
  linop_d->thread_barrier();
  axpy(Krylov_r,v,v,0.0);
  axpy(Krylov_p,v,v,0.0);

  ApplyThread(Krylov_p,Krylov_Ap); //Threaded
  ApplyThread(Krylov_r,Krylov_Ar); //Threaded
  rAr = innerProductReal(Krylov_r,Krylov_Ar);
  pAAp= norm_vec(Krylov_Ap);

  cp =norm_vec(Krylov_r);
  ssq=norm_vec(v);
  rsq = LittleDopSolverResidual*LittleDopSolverResidual*ssq;

  uint64_t t1 = GetTimeBase();
  for(int k=1;k<10000;k++){


    MGreport=0;
    //    if(linop_d->isBoss() && (k==10) ) MGreport=1;
    a   = rAr/pAAp;

    uint64_t axpy_0 = GetTimeBase();
    axpy(vinv,Krylov_p,vinv,a);           // overlap with comms
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a);
    uint64_t axpy_1 = GetTimeBase();


    uint64_t mat_0 = GetTimeBase();
    ApplyThread(Krylov_r,Krylov_Ar); //          Ar_i+1
    uint64_t mat_1 = GetTimeBase();

    rArp=rAr;
    uint64_t inner_0 = GetTimeBase();
    rAr = innerProductReal(Krylov_r,Krylov_Ar); // Fuse into ApplyThread
    uint64_t inner_1 = GetTimeBase();
    b = rAr/rArp;
 
    uint64_t axpy1_0 = GetTimeBase();
    axpy(Krylov_p,Krylov_p,Krylov_r,b);
    axpy(Krylov_Ap,Krylov_Ap,Krylov_Ar,b); // Axpy_norm
    uint64_t axpy1_1 = GetTimeBase();

    uint64_t norm_0 = GetTimeBase();
    pAAp= norm_vec(Krylov_Ap);// Fuse into axpy norm
    cp  = norm_vec(Krylov_r); // Fuse into earlier ApplyThread
    uint64_t norm_1 = GetTimeBase();

    if(cp<rsq) {
      ApplyThread(vinv,Krylov_Ap);
      axpy(Krylov_Ap,v,Krylov_Ap,-1.0);

      double nv = norm_vec(Krylov_Ap);
      uint64_t t2 = GetTimeBase();


      return;
    }
  }
  linop_d->Error("Little Dirac Matrix inversion failed\n");
  fflush(stdout);
  exit(0);
}
// Threaded
template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverseMultiShift(std::vector<std::complex<cFloat> >& src,
					     std::vector<std::vector<std::complex<cFloat> > >& psi,
					     std::vector<std::vector<std::complex<cFloat> > >& ps,
					     int nshift,
					     double *mass,
					     double *alpha,
					     double *residual )
{
  int me = linop_d->thread_barrier();
  int max_iter = 10000;

  // Per shift fields
  double    bs [nshift];
  double    rsq[nshift];
  double    z[nshift][2];
  int       converged[nshift];

  const int       primary =0;

  //Primary shift fields CG iteration
  double a,b,c,d;
  double cp,bp; //prev

  struct timeval start,stop;

  // Check lightest mass
  for(int s=0;s<nshift;s++){
    if ( mass[s] < mass[primary] ) {
      linop_d->Error("First shift not lightest - oops\n");
      exit(-1);
    }
  }

  //Allocate per shift fields
  for(int i=0;i<nshift;i++){
    converged[i]=0;
  }

  // Wire guess to zero
  // Residuals "r" are src
  // First search direction "p" is also src
  cp = norm_vec(src);
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * residual[s] * residual[s];
    linop_d->ThreadBossDebug("ApplyInverseMultiShift: shift %d target resid %le %le %le\n",
			     s,rsq[s],cp,residual[s]);
    copy(ps[s],src);
  }
  // r and p for primary
  scale(Krylov_r,src,1.0);
  scale(Krylov_p,src,1.0);
  
  //MdagM+m[0]

  ApplyThread(Krylov_p,Krylov_Ap); //Threaded
  axpy(Krylov_Ap,Krylov_p,Krylov_Ap,mass[0]);
  d=innerProductReal(Krylov_p,Krylov_Ap);
  b = -cp /d;

  // Set up the various shift variables
  int       iz=0;
  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; 
  }
  
  // c= norm(r)
  axpy(Krylov_r,Krylov_Ap,Krylov_r,b);
  c = norm_vec(Krylov_r);

  for(int s=0;s<nshift;s++) {
    axpby(psi[s],src,src,0.,-bs[s]*alpha[s]);
  }

  int continue_all_vecs=1;
  // Iteration loop
  int k;
  gettimeofday(&start,NULL);
  for (k=1;k<=max_iter;k++){

    a = c /cp;
    axpy(Krylov_p,Krylov_p,Krylov_r,a);

    for(int s=0;s<nshift;s++){
      if ( (!converged[s]) || continue_all_vecs ) { 
        if (s==0){
          axpy(ps[s],ps[s],Krylov_r,a);
        } else{
	  double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  axpby(ps[s],Krylov_r,ps[s],z[s][iz],as);
        }
      }
    }

    cp=c;
    
    ApplyThread(Krylov_p,Krylov_Ap); //Threaded
    axpy(Krylov_Ap,Krylov_p,Krylov_Ap,mass[0]);
    d=innerProductReal(Krylov_p,Krylov_Ap);

    bp=b;
    b=-cp/d;
    
    axpy(Krylov_r,Krylov_Ap,Krylov_r,b);
    c=norm_vec(Krylov_r);

    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;
    for(int s=1;s<nshift;s++){
      if((!converged[s]) || continue_all_vecs){
	double z0 = z[s][1-iz];
	double z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    for(int s=0;s<nshift;s++){
      int ss = s;
      if((!converged[s])||continue_all_vecs)
	axpy(psi[ss],ps[s],psi[ss],-bs[s]*alpha[s]);
    }

    // Convergence checks
    if ( (k%500)==0) {
      linop_d->ThreadBossMessage("ApplyInverseMultiShift: k=%d c=%g\n",k,c);
    }

    int all_converged = 1;
    for(int s=0;s<nshift;s++){
      double css  = c * z[s][iz]* z[s][iz];
      if(css<rsq[s]){
	converged[s]=1;
	linop_d->ThreadBossDebug("ApplyInverseMultiShift: k=%d Shift %d has converged\n",k,s);
      } else {
	all_converged=0;
      }
    }

    if ( all_converged ){

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);
      double t = diff.tv_sec*1.0E6 + diff.tv_usec;

      
      linop_d->ThreadBossMessage("ApplyInverseMultiShift: k=%d All shifts have converged after %le s\n",k,t*1.0e-6);
      linop_d->ThreadBossDebug("ApplyInverseMultiShift: k=%d Checking solutions\n",k);
      // Check answers 
      for(int s=0; s < nshift; s++) { 
	ApplyThread(psi[s],Krylov_Ap); //Threaded
	axpy(Krylov_Ap,psi[s],Krylov_Ap,mass[s]);
	axpy(Krylov_r,Krylov_Ap,src,-1);
	double rn = norm_vec(Krylov_r);
	double cn = norm_vec(src);
	linop_d->ThreadBossMessage("ApplyInverseMultiShift: shift[%d] true residual %le \n",s,sqrt(rn/cn));
      }

      return;
    }
  }

  linop_d->ThreadBossDebug("ApplyInverseMultiShift: CG not converged after %d iterations\n",k);

  return ;
}

// Threaded
template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverseCG(std::vector<std::complex<cFloat> > &v,std::vector<std::complex<cFloat> > &vinv, 
				     double shift,double resid, int maxit)
{
  int me, thrlen,throff;

  // Conjugate gradient on A. A is hermitian posdef.
  // Need these to be shared across threads
  me = linop_d->thread_barrier();

  axpby(vinv,v,v,0.0,0.0);
  axpy(Krylov_r,v,v,0.0);
  axpy(Krylov_p,v,v,0.0);

  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;

  ssq=norm_vec(v);
  a  =norm_vec(Krylov_p);
  cp =norm_vec(Krylov_r);

  // replace default args
  if ( resid == 0 ) resid = LittleDopSolverResidual;
  if ( maxit == 0 ) maxit = 10000;
  rsq = resid*resid*ssq;

  uint64_t t0 = GetTimeBase();
  for(int k=1;k<maxit;k++){
    c=cp;

    MGreport=0;
    //if(linop_d->isBoss() && (k==10) ) MGreport=1;
    
    uint64_t mat_0 = GetTimeBase();
    ApplyThread(Krylov_p,Krylov_Ap); //Threaded
    if ( shift != 0.0 ) { 
      axpy(Krylov_Ap,Krylov_p,Krylov_Ap,shift);
    }
    uint64_t mat_1 = GetTimeBase();

    uint64_t inner_0 = GetTimeBase();
    d=innerProductReal(Krylov_p,Krylov_Ap);
    a=c/d;
    uint64_t inner_1 = GetTimeBase();

    uint64_t axpy_norm_0 = GetTimeBase();
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a); // Fuse Axpy-Norm
    cp=norm_vec(Krylov_r);
    uint64_t axpy_norm_1 = GetTimeBase();
    b=cp/c;


    if ( (k%200==0) ) { 
      linop_d->ThreadBossDebug("ApplyInverseCG:: k=%d cp= %le b = %le d=%le a=%le\n",k,cp,b,d,a);
    }
    uint64_t axpy_0 = GetTimeBase();
    axpy(vinv,Krylov_p,vinv,a);
    axpy(Krylov_p,Krylov_p,Krylov_r,b);
    uint64_t axpy_1 = GetTimeBase();


    if(cp<rsq) {
      uint64_t t1 = GetTimeBase();
      ApplyThread(vinv,Krylov_Ap);
      if ( shift != 0.0 ) { 
	axpy(Krylov_Ap,vinv,Krylov_Ap,shift);
      }
      linop_d->thread_work(LocalNsubspace,me,thrlen,throff);
      axpy(Krylov_Ap,v,Krylov_Ap,-1.0);

      double nv = norm_vec(Krylov_Ap);
      linop_d->ThreadBossMessage("HDCG LittleDiracInversion via CG converged : %d iterations , residual %g, ssq %g, millisec %f shift %f\n",k,
	       sqrt(nv/ssq),sqrt(ssq),(t1-t0)/MHz()/1000.0, shift);
      
      return;
    }
  }

  printf("Little Dirac Matrix inversion failed \n");
  exit(0);
}

template<class cFloat> 
void BfmHDCG<cFloat>::ApplyInverseCGopt(std::vector<std::complex<cFloat> > &v,std::vector<std::complex<cFloat> > &vinv, 
					std::vector<double> & polynomial,
					double shift,double resid, int maxit, int halo)
{
  int me, thrlen,throff;

  // Conjugate gradient on A. A is hermitian posdef.
  // Need these to be shared across threads
  me = linop_d->thread_barrier();

  
  axpby(vinv,v,v,0.0,0.0);
  axpy(KrylovNest_r,v,v,0.0);
  axpy(KrylovNest_p,v,v,0.0);

  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;
  int single = 1;

  static int count;
  if ( !me ) { 
    count ++;
  }
  ssq=norm_vec(v);
  a  =norm_vec(KrylovNest_p);
  cp =norm_vec(KrylovNest_r);

  // replace default args
  if ( resid == 0 ) resid = LittleDopSolverResidual;
  if ( maxit == 0 ) maxit = 10000;
  if ( resid < 1.0e-6 ) single =0;
  rsq = resid*resid*ssq;

  std::vector<double> poly_p(1);
  std::vector<double> poly_r(1);
  std::vector<double> poly_Ap;

  if ( !me ) { 
    poly_p[0]=1.0;
    poly_r[0]=1.0;
    polynomial.resize(0);
  }


  uint64_t t0 = GetTimeBase();
  for(int k=1;k<maxit;k++){
    c=cp;

    MGreport=0;
    //if(linop_d->isBoss() && (k==10) ) MGreport=1;
    
    uint64_t mat_0 = GetTimeBase();
    ApplyThreadOpt(KrylovNest_p,KrylovNest_Ap,single,halo); //Threaded
    if ( shift != 0.0 ) { 
      axpy(KrylovNest_Ap,KrylovNest_p,KrylovNest_Ap,shift);
    }
    uint64_t mat_1 = GetTimeBase();

    uint64_t inner_0 = GetTimeBase();
    d=innerProductReal(KrylovNest_p,KrylovNest_Ap);
    a=c/d;
    uint64_t inner_1 = GetTimeBase();

    uint64_t axpy_norm_0 = GetTimeBase();
    axpy(KrylovNest_r,KrylovNest_Ap,KrylovNest_r,-a); // Fuse Axpy-Norm
    cp=norm_vec(KrylovNest_r);
    uint64_t axpy_norm_1 = GetTimeBase();
    b=cp/c;


    if(!me){
      //  Ap= right_shift(p)
      poly_Ap.resize(k+1);
      poly_Ap[0]=0.0;
      for(int i=0;i<k;i++){
	poly_Ap[i+1]=poly_p[i];
      }

      //  x = x + a p
      polynomial.resize(k);
      polynomial[k-1]=0.0;
      for(int i=0;i<k;i++){
	polynomial[i] = polynomial[i] + a * poly_p[i];
      }
    
      //  r = r - a Ap
      //  p = r + b p
      poly_r.resize(k+1);
      poly_p.resize(k+1);
      poly_r[k] = poly_p[k] = 0.0;
      for(int i=0;i<k+1;i++){
	poly_r[i] = poly_r[i] - a * poly_Ap[i];
	poly_p[i] = poly_r[i] + b * poly_p[i];
      }
    }


    if ( (k%200==0) ) { 
      linop_d->ThreadBossDebug("ApplyInverseCG:: k=%d cp= %le b = %le d=%le a=%le\n",k,cp,b,d,a);
    }
    uint64_t axpy_0 = GetTimeBase();
    axpy(vinv,KrylovNest_p,vinv,a);
    axpy(KrylovNest_p,KrylovNest_p,KrylovNest_r,b);
    uint64_t axpy_1 = GetTimeBase();


    if(cp<rsq) {
      uint64_t t1 = GetTimeBase();
      ApplyThreadOpt(vinv,KrylovNest_Ap,single,halo);
      if ( shift != 0.0 ) { 
	axpy(KrylovNest_Ap,vinv,KrylovNest_Ap,shift);
      }
      linop_d->thread_work(LocalNsubspace,me,thrlen,throff);
      axpy(KrylovNest_Ap,v,KrylovNest_Ap,-1.0);
#ifdef DEBUG_HDCG
      double nv = norm_vec(KrylovNest_Ap);
      if ( count < 50 ) {
	linop_d->ThreadBossMessage("HDCG LittleDiracInversion via CGopt<%d> converged : %d iterations , residual %g, ssq %g, millisec %f shift %f\n",
				   halo,k,
				   sqrt(nv/ssq),sqrt(ssq),(t1-t0)/MHz()/1000.0, shift);
      }
      if(1){
      if ( !me && linop_d->isBoss() ){
	printf("f(x) = ");
	for(int i=0;i<polynomial.size();i++){
	  if ( i!= 0 ) printf("+");
	  printf("(%.16f) * pow(x,%d)",polynomial[i],i);
	}
	printf("\n");
      }
      }
#endif
      return;
    }
  }

#ifdef DEBUG_HDCG
  if(1){
  if ( !me && linop_d->isBoss() ){
    printf("f(x) = ");
    for(int i=0;i<polynomial.size();i++){
      if ( i!= 0 ) printf("+");
      printf("(%.16f) * pow(x,%d)",polynomial[i],i);
    }
    printf("\n");
  }
  }
  if ( count < 50 ) {
    linop_d->ThreadBossMessage("HDCG LittleDiracInversion via CGopt<%d> not converged : residual %g\n",
			       halo, sqrt(cp/ssq));
  }
#endif
}

void mycgemv(int N,float * __restrict A, double * __restrict in,double c, double * __restrict out)
{
  if (((uint64_t) A) &0x1f) { printf("unaligned A\n"); fflush(stdout); while(1){};}
  if (((uint64_t)in) &0x1f) { printf("unaligned in\n"); fflush(stdout);while(1){};}
  if (((uint64_t)out)&0x1f) { printf("unaligned out\n"); fflush(stdout);while(1){};}
  qpx_cgemv(N,A,in,c,out);
}
void mycgemv_f(int N,float * __restrict A, float * __restrict in,float c, float * __restrict out)
{
  if (((uint64_t) A) &0x1f) { printf("unaligned A\n"); fflush(stdout); while(1){};}
  if (((uint64_t)in) &0x1f) { printf("unaligned in\n"); fflush(stdout);while(1){};}
  if (((uint64_t)out)&0x1f) { printf("unaligned out\n"); fflush(stdout);while(1){};}
  qpx_cgemv_f(N,A,in,c,out);
}  

void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out)
{
  if (((uint64_t) A) &0x1f) { printf("unaligned A\n"); fflush(stdout); while(1){};}
  if (((uint64_t)in) &0x1f) { printf("unaligned in\n"); fflush(stdout);while(1){};}
  if (((uint64_t)out)&0x1f) { printf("unaligned out\n"); fflush(stdout);while(1){};}
  qpx_zgemv(N,A,in,c,out);
}


int gcd(int a, int b)
{
  int c = a % b;
  while(c != 0)
    {
      a = b;
      b = c;
      c = a % b;
    }
  return b;
}

template<class cFloat> 
void BfmHDCG<cFloat>::ApplyDiagInv(std::vector<std::complex<cFloat> > &in,
				   std::vector<std::complex<cFloat> > &out)
{
  int me, thrlen,throff;
  int nwork     = LocalNblock;
  linop_d->thread_work(nwork,me,thrlen,throff);
  for(int w=throff;w<throff+thrlen;w++){

    double  coeff=0;
    if ( sizeof(cFloat)==sizeof(double) ) {
      double *oo= (double *)&out[w*Nvec*Nvec];
      double *ii= (double *)&in[w*Nvec*Nvec];
      double *Aps=(double *)&AsparseDiagInv[w*Nvec*Nvec];
      myzgemv(Nvec,Aps,ii,coeff,oo);
    } else { 
      float  *oo= (float  *)&out[w*Nvec*Nvec];
      float  *ii= (float  *)&in[w*Nvec*Nvec];
      float  *Aps=(float  *)&AsparseDiagInv[w*Nvec*Nvec];
      mycgemv_f(Nvec,Aps,ii,coeff,oo);
    }
  }
  linop_d->thread_barrier();
}

template<class cFloat> 
void BfmHDCG<cFloat>::ApplyThreadOpt(std::vector<std::complex<cFloat> > &in,
				     std::vector<std::complex<cFloat> > &out,
				     int single ,int depth)
{
  int me, thrlen,throff;

  double nA=0.0;
  double nH=0.0;

  linop_d->thread_barrier();

  uint64_t t0 = GetTimeBase();
  HaloExchange(in,depth);                
  uint64_t t1 = GetTimeBase();

  int NballDeep = StencilDepthCount[depth];
  
  int nwork     = LocalNblock*NballDeep;

  linop_d->thread_work(nwork,me,thrlen,throff);

  // Thread across both blocks and directions jointly to maximally load balance.
  for(int w=throff;w<throff+thrlen;w++){

    int  b = w/NballDeep;    // Which block
    int mmu= w%(NballDeep);  // Which neigbour
    int mu = StencilReord[mmu];

    cFloat *oo= (cFloat *)&Krylov_Atmp[mmu*LocalNsubspace+b*Nvec];
    cFloat *nbr_data = (cFloat *)HaloGetStencilData(b,mmu,in);
    double coeff=0.0;
    if ( sizeof(cFloat)==sizeof(double) ) { 
      if ( single ) { 
	float  *Aps=(float *) &AsparseSingle[b*Nball*Nvec*Nvec+mu*Nvec*Nvec];
	if ( StencilNonZero[mu] && (StencilDepth[mu]<=depth) ) {
	  mycgemv(Nvec,Aps,(double *)nbr_data,coeff,(double *)oo);
	}
      } else { 
	double  *Aps=(double *) &Asparse[b*Nball*Nvec*Nvec+mu*Nvec*Nvec];
	if ( StencilNonZero[mu] && (StencilDepth[mu]<=depth) ) {
	  myzgemv(Nvec,Aps,(double *)nbr_data,coeff,(double *)oo);
	}
      }
    } else { 
      float  *Aps=(float *) &Asparse[b*Nball*Nvec*Nvec+mu*Nvec*Nvec];
      if ( StencilNonZero[mu] && (StencilDepth[mu]<=depth) ) {
	mycgemv_f(Nvec,Aps,(float *)nbr_data,coeff,(float *)oo);
      }
    }
  }
  
  linop_d->thread_barrier();
  uint64_t t2 = GetTimeBase();

  // Reduce the result vectors down.
  linop_d->thread_work(LocalNsubspace,me,thrlen,throff);
  double coeff = 0.0;
  for(int mmu=0;mmu<NballDeep;mmu++){
    if ( sizeof(cFloat) == sizeof(double) ) {
      qpx_axpby(thrlen,
		(double *)&out[throff],
		(double *)&out[throff],
		(double *)&Krylov_Atmp[mmu*LocalNsubspace+throff],
		coeff,1.0);
    }else { 
      qpx_axpby_f(thrlen,
		  (float *)&out[throff],
		  (float *)&out[throff],
		  (float *)&Krylov_Atmp[mmu*LocalNsubspace+throff],
		  coeff,1.0);
    }
    coeff = 1.0;
  }
  linop_d->thread_barrier();
#ifdef DEBUG_HDCG

  uint64_t t3 = GetTimeBase();
  if (  MGreport  ) {
    linop_d->ThreadBossMessage("ApplyThreadOpt -- Halo    : depth %d t1-t0 cycles %ld, %f us\n",depth,t1-t0,1.0*(t1-t0)/MHz());
    linop_d->ThreadBossMessage("ApplyThreadOpt -- MatMul  : reduce cycles %ld, mult cycles %ld single precision %d\n",t3-t2,t2-t1,single);
    linop_d->thread_barrier();
  }
#endif
}
template<class cFloat> 
void BfmHDCG<cFloat>::ApplyThread(std::vector<std::complex<cFloat> > &in,
				  std::vector<std::complex<cFloat> > &out,
				  int single,
				  int depth)
{
  int me, thrlen,throff;
  double nA=0.0;
  double nH=0.0;
  me = linop_d->thread_barrier();


  // Going async comms  here => 2x in ldop.
  // Only way to avoid the blocking here is to

  // Phase-1 
  //      start the comms for half (A) of the data, 
  //      compute multiply for half (B) of the data.
  //
  // Phase-2
  //      start the comms for half (B) of the data, 
  //      compute multiply for half (A) of the data.
  // 
  // Phase-3 
  //      add half B of data.
  //
  // If done correctly there is a factor of two speed up in the Ldop.
  // This is 2/3+ of Ldop CG time.
  // Conclude that we get around 25% and 190->140s. 
  //
  //
  uint64_t t0 = GetTimeBase();
  HaloExchange(in,depth);                
  uint64_t th = GetTimeBase();

  double flops = Nvec*Nvec*Nball*8*LocalNblock;
  int nwork = LocalNblock*Nball;

  int threads   = linop_d->threads;
  //  int block_threads= gcd(LocalNblock,threads);
  //  int ball_threads = threads/block_threads;
  int ball_threads = 8;
  if ( depth < 4 )  ball_threads=4;

  int block_threads = threads/ball_threads;
  int block_me =me/ball_threads;
  int block_len;
  int block_off;
  {
    ThreadModelSingle WorkAllocator; 
    WorkAllocator.nthread=block_threads;
    WorkAllocator.thread_work_nobarrier(LocalNblock,block_me,block_len,block_off);
  }
  
  int ball_me  =me%ball_threads;
  int ball_len;
  int ball_off;
  {
    ThreadModelSingle WorkAllocator; 
    WorkAllocator.nthread=ball_threads;
    WorkAllocator.thread_work_nobarrier(Nball,ball_me,ball_len,ball_off);
  }
  linop_d->thread_barrier();

  uint64_t t1 = GetTimeBase();
  uint64_t tt;
  uint64_t t3;

  //  double nn=norm_vec(in);
  linop_d->thread_barrier();
  tt = GetTimeBase();
  for(int b =block_off;b<block_off+block_len;b++){

    int w = ball_me*LocalNblock+b;
    cFloat *oo= (cFloat *)&Krylov_Atmp[w*Nvec];
    
    double coeff=0.0;

    for(int mmu= ball_off;mmu<ball_off+ball_len;mmu++){
      
      int mu = StencilReord[mmu];
      cFloat *nbr_data = (cFloat *)HaloGetStencilData(b,mmu,in);
      if ( sizeof(cFloat)==sizeof(double) ) { 
	if ( single ) { 
	  float  *Aps=(float *) &AsparseSingle[b*Nball*Nvec*Nvec+mu*Nvec*Nvec];
	  if ( StencilNonZero[mu] && (StencilDepth[mu]<=depth) ) {
	    mycgemv(Nvec,Aps,(double *)nbr_data,coeff,(double *)oo);
	    coeff=1.0;
	  }
	} else { 
	  double  *Aps=(double *) &Asparse[b*Nball*Nvec*Nvec+mu*Nvec*Nvec];
	  if ( StencilNonZero[mu] && (StencilDepth[mu]<=depth) ) {
	    myzgemv(Nvec,Aps,(double *)nbr_data,coeff,(double *)oo);
	    coeff=1.0;
	  }
	}
      } else { 
	float  *Aps=(float *) &Asparse[b*Nball*Nvec*Nvec+mu*Nvec*Nvec];
	if ( StencilNonZero[mu] && (StencilDepth[mu]<=depth) ) {
	  mycgemv_f(Nvec,Aps,(float *)nbr_data,coeff,(float *)oo);
	  coeff=1.0;
	}
      }
    }
  }
  linop_d->thread_barrier();
  t3 = GetTimeBase();

  nwork = LocalNblock;
  linop_d->thread_work(nwork,me,thrlen,throff);
  for(int b=throff;b<throff+thrlen;b++){
    std::complex<double> ctmp;
    for(int i=0;i<Nvec;i++){
      int o = b*Nvec+i;
      ctmp= Krylov_Atmp[o];
      for(int mu=1;mu<ball_threads;mu++){
	ctmp+= Krylov_Atmp[mu*Nvec*LocalNblock+o];
      }
      out[o]=ctmp;
    }
  }
  linop_d->thread_barrier();
#ifdef DEBUG_HDCG
  uint64_t t2 = GetTimeBase();
  if (  MGreport  ) {
    linop_d->ThreadBossMessage("ApplyThread -- Halo<%d>: th-t0 cycles %ld, %f us\n",depth,th-t0,1.0*(th-t0)/MHz());
    linop_d->ThreadBossMessage("ApplyThread -- MatMul  : t3-tt cycles %ld, %f Mflop/s %f us\n",t3-tt,flops*MHz()/(t3-tt),1.0*(t3-tt)/MHz());
    linop_d->ThreadBossMessage("ApplyThread -- MatMul  : t2-t1 cycles %ld, %f Mflop/s %f us\n",t2-t1,flops*MHz()/(t2-t1),1.0*(t2-t1)/MHz());
    linop_d->ThreadBossMessage("ApplyThread -- MatMul  : reduce cycles %ld, mult cycles %ld single precision %d\n",t2-t3,t3-t1,single);
    linop_d->thread_barrier();
  }
#endif
}

/////////////////////////////////////////
// Deflated Solver support
/////////////////////////////////////////

// Flaws here:
// -- Getting expensive on 48^3 and 4^4 blocks. 6h on 256 nodes unacceptable. 
// -- Initial tests [16^3] indicate reducing precision a lot [10^-7] doesn't help
// -- Ideally go to an intermediate grid and block the vectors we generate here. But challenging
// -- Reduce precision. Use earlier vectors to deflate later vectors => accelerate generation => 5x-10x ??
// -- Use single precision matrix storage in the generation => 2x 
// -- Multishift rational solve => ~ 3-6x.
//
// Note: could use fewer vectors, block them, and store as dense matrix still. 1 block per node??
//
// Plan ahead:
// i)   Identify limits of precision on LittleDopInverse and SubspaceGeneration for finest level.
//      Present seqence of jobs should do this.
//      1.0e-8  on inverse iteration.
//
//    a)      Single prec mat storage OK.
//             Done
//
//    b)      Could single project the vectors in ldop
//
//    c)      Could use single in preconditioning inverse [M+a]
//
//    d)      Could optimise the multishift CG.
//
//    e)      Profile code and assess whether I can grow subspace.
//            Grow this from residual in some way?
//            Inverse iteration from residual?
//
//    f)      Use smaller Ls for subspace generation?
//            Even for a portion of the space??
//
//    
// ii)  Use 16^3 to investigate trade off of 2^4,3^4 quadrants and fewer vectors. 
//      Q. Can I decrease vectors for free?
// 16^3
//           2^2x4^2 and  9 vecs => 140 iterations
//           4^4     and 36 vecs =>  35 iterations.
//           2^2x4^2 and 32 vecs =>  31 iterations  - factor of 4 in blockvol, small factor in convergence
//
// 6^4      36 Deep 1.0: => 486
// 4,4,6,6  36 Deep 1.0: => 399              
// 4,4,4,4, 36 Deep 1.0: => 317    - factor of 4 in block vol <=> 5/3 in iters.
//
//      A. No. More vectors better than subdivision of same vectors.
//                  
//      Information also of use to predict behaviour of next level where I use 128 vecs
//    
// iii) Reduce cost of next level subspace generation. [subroutine below]
//  Y   - Reduce precision 
//      1.0e-8 ok. Not sure how far this can be pushed. Should test this.
//
//
//  Y   - using float (single prec) rep of matrix.     
//      - No effect. Use it.
//
//  Y    - Create 12 vecs and use these to deflate rest of subspace gen
//         - Make the DeflationBasisInit "grow" the size and then use the deflated CG in this.
//      - Worked on 16^3.
//
//      - Multishift rational filter
//
//      - reduce vecs
//      - proper RNGs
//
// iv) Blocking next level subspace => fewer vectors [32?]. => bigger matrix to treat ~ num nodes x num nodes.
//     Blocking flexibly, beyond level of node may make sense.
//
template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationBasisTrivial(void)
{
  int No= LdopDeflationBasisSize;
  int N = No+Nvec;
  int Nbasis=N;
  int Ns= LocalNsubspace;

  std::vector<std::vector<std::complex<cFloat> > > DeflationBasis = LdopDeflationBasis;
  DeflationBasis.resize(N);
  for(int v=No;v<N;v++){
    DeflationBasis[v].resize(LocalNsubspace);
  }

  linop_d->BossLog("LdopDeflationBasisTrivial: Adding trivially obtained vectors for little dirac op %d -> %d\n",No,N);

  // Orthogonalise subspace
  // Remove earlier (normalised) vectors
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      for(int i=0;i<Nvec;i++){
	ProjectToSubspace(subspace_r[i],DeflationBasis[i+No]);
      }
    }
  }

  linop_d->BossLog("LdopDeflationBasisTrivial: Projected to subspace\n");

  LdopDeflationBasisSize = Nbasis;
  LdopDeflationMatrix.resize(N*N);
  LdopDeflationInverseMatrix.resize(N*N);
  LdopDeflationAv.resize(N);
  LdopDeflationBasis.resize(N);
  for(int v=No;v<N;v++){
    LdopDeflationBasis[v].resize(LocalNsubspace);
#pragma omp parallel 
    {
#pragma omp for 
      for(int t=0;t<linop_d->nthread;t++) {
	axpy(LdopDeflationBasis[v],DeflationBasis[v],DeflationBasis[v],0.0);
      }
    }
    LdopDeflationAv[v].resize(LocalNsubspace);
  }  

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      std::complex<double> c;
      double n,a;

      for(int v=0;v<N;v++){
	for(int u=0;u<v;u++){
	  linop_d->thread_barrier();

	  //Inner product & remove component
	  c = innerProduct(LdopDeflationBasis[u],LdopDeflationBasis[v]);
	  zaxpy(LdopDeflationBasis[v],LdopDeflationBasis[u],LdopDeflationBasis[v],-c);
	}
	// Normalise this vector
	n = norm_vec(LdopDeflationBasis[v]);
	a = 1.0/sqrt(n);
	axpby(LdopDeflationBasis[v],LdopDeflationBasis[v],LdopDeflationBasis[v],a,0.0);
      }
    }
  }

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      // Compute matrix elements
      for(int v=0;v<N;v++){
	linop_d->thread_barrier();
	ApplyThread(LdopDeflationBasis[v],LdopDeflationAv[v]);
	for(int u=0;u<N;u++){
	  std::complex<double> val          = innerProduct(LdopDeflationBasis[u],LdopDeflationAv[v]);
	  LdopDeflationMatrix[N*u+v]        =val;
	  LdopDeflationInverseMatrix[N*u+v] =val;
	}
      }
    }
  }

  linop_d->BossLog("LdopDeflationBasisInit: Calling out to LAPACK for dense matrix inversion\n");
  // Precompute inverse matrix rep on subspace
  // This could be made very efficient by LDU, and change basis as only N zaxpy's required to apply

  LapackHermitianInvert(N,(double *)&LdopDeflationInverseMatrix[0]);


}

template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationBasisInit(int Nbasis)
{
  int No= LdopDeflationBasisSize;
  int N = Nbasis;
  int Ns= LocalNsubspace;

  std::vector<std::vector<std::complex<cFloat> > > DeflationBasis = LdopDeflationBasis;
  DeflationBasis.resize(N);
  for(int v=No;v<N;v++){
    DeflationBasis[v].resize(LocalNsubspace);
  }

  linop_d->BossLog("LdopDeflationBasisInit: Building deflation space for little dirac op\n");
  //  std::normal_distribution<double> dist;
  //  std::mt19937 eng;
  for(int v=No;v<N;v++){ 

    //////////////////////////////////////////////////////
    // Gaussian noise???
    //////////////////////////////////////////////////////
    std::vector<std::complex<cFloat> > src(LocalNsubspace);
    std::vector<std::complex<cFloat> > src_p(LocalNsubspace);
    for(int i=0;i<LocalNsubspace;i++){
      src[i] = std::complex<cFloat>(drand48()-0.5,drand48()-0.5);
      src[i] = src[i]*((cFloat)1.0/LocalNsubspace);
    }

    ////////////////////////////////////////////////////////////
    // Remove existing subspace from source.
    ////////////////////////////////////////////////////////////
    if(LdopDeflationBasisSize >0){
      std::vector<std::complex<cFloat> > ss(LdopDeflationBasisSize);
      if ( linop_d->SPIcomms() ) linop_d->comm_init();
#pragma omp parallel 
      {
#pragma omp for 
	for(int t=0;t<linop_d->nthread;t++) {
	  LdopDeflationProject      (src,ss);
	  LdopDeflationPromote      (ss,src_p);
	  axpy(src,src_p,src,-1.0);
	}
      }
    }

    if ( LittleDopSubspaceRational == 0 ) { 

    // Try this, and also force the solver to use the single precision matrix???
    // Can also make the residual tight only on last solve.
#pragma omp parallel 
    {
#pragma omp for 
      for(int t=0;t<linop_d->nthread;t++) {
	for(int i=0;i<3;i++){

	  LittleDopSolverResidual= LittleDopSolverResidualSubspace;

	  linop_d->thread_barrier();
	  if ( No == 0 ) {
	    ApplyInverseCG(src,DeflationBasis[v],0.0e-3);
	  } else { 
	    ApplyInverseDeflCG(src,DeflationBasis[v]);
	  }
	  double nn=norm_vec(DeflationBasis[v]);
	  double scale=1.0/sqrt(nn);
	  axpby(DeflationBasis[v],DeflationBasis[v],DeflationBasis[v],scale,0.0);
	  axpy(src,DeflationBasis[v],DeflationBasis[v],0.0);
	}
      }
    }

    } else { 

      int Nrat=4;
      std::vector<std::vector<std::complex<cFloat> >  > ps(Nrat);
      std::vector<std::vector<std::complex<cFloat> >  > sols(Nrat);
      for(int r=0;r<Nrat;r++){
	ps[r].resize(LocalNsubspace);
	sols[r].resize(LocalNsubspace);
      }
      double Lo =SubspaceRationalLo/2;
      double epsilon      = Lo/3;
      double residuals[4]     = {LittleDopSolverResidualSubspace,
				 LittleDopSolverResidualSubspace,
				 LittleDopSolverResidualSubspace,
				 LittleDopSolverResidualSubspace
      };
      double alpha[4]     = {1.0,1.0,1.0,1.0};
      double shifts[4]    = {Lo,Lo+epsilon,Lo+2*epsilon,Lo+3*epsilon};

#pragma omp parallel 
      {
#pragma omp for 
	for(int t=0;t<linop_d->nthread;t++) {
	  
	  for(int r=0;r<Nrat;r++){
	    zeroOut(ps[r]);
	    zeroOut(sols[r]);
	  }
	  
	  ApplyInverseMultiShift(src,sols,ps,Nrat,(double *)shifts,(double *)alpha,(double *)residuals);
	  axpby(DeflationBasis[v],sols[0],sols[1],1.0/6.0,-1.0/2.0);
	  axpy(DeflationBasis[v],sols[2],DeflationBasis[v],1.0/2.0);
	  axpy(DeflationBasis[v],sols[3],DeflationBasis[v],-1.0/6.0);
	}
      }
    }
  }

  LdopDeflationBasisSize = Nbasis;
  LdopDeflationMatrix.resize(N*N);
  LdopDeflationInverseMatrix.resize(N*N);
  LdopDeflationAv.resize(N);
  LdopDeflationBasis.resize(N);
  for(int v=No;v<N;v++){
    LdopDeflationBasis[v].resize(LocalNsubspace);
#pragma omp parallel 
    {
#pragma omp for 
      for(int t=0;t<linop_d->nthread;t++) {
	axpy(LdopDeflationBasis[v],DeflationBasis[v],DeflationBasis[v],0.0);
      }
    }
    LdopDeflationAv[v].resize(LocalNsubspace);
  }  

  linop_d->BossLog("LdopDeflationBasisInit: Orthogonalising subspace\n");

  // Orthogonalise subspace
  // Remove earlier (normalised) vectors
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {

      std::complex<double> c;
      double n,a;

      for(int v=0;v<N;v++){
	for(int u=0;u<v;u++){
	  linop_d->thread_barrier();

	  //Inner product & remove component
	  c = innerProduct(LdopDeflationBasis[u],LdopDeflationBasis[v]);
	  zaxpy(LdopDeflationBasis[v],LdopDeflationBasis[u],LdopDeflationBasis[v],-c);
	}
	// Normalise this vector
	n = norm_vec(LdopDeflationBasis[v]);
	a = 1.0/sqrt(n);
	axpby(LdopDeflationBasis[v],LdopDeflationBasis[v],LdopDeflationBasis[v],a,0.0);
      }
    }
  }

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      // Compute matrix elements
      for(int v=0;v<N;v++){
	linop_d->thread_barrier();
	ApplyThread(LdopDeflationBasis[v],LdopDeflationAv[v]);
	for(int u=0;u<N;u++){
	  std::complex<double> val          = innerProduct(LdopDeflationBasis[u],LdopDeflationAv[v]);
	  LdopDeflationMatrix[N*u+v]        =val;
	  LdopDeflationInverseMatrix[N*u+v] =val;
	}
      }
    }
  }

  linop_d->BossLog("LdopDeflationBasisInit: Calling out to LAPACK for dense matrix inversion\n");
  // Precompute inverse matrix rep on subspace
  // This could be made very efficient by LDU, and change basis as only N zaxpy's required to apply

  LapackHermitianInvert(N,(double *)&LdopDeflationInverseMatrix[0]);

  // Check the inverse
  double maxerr=0.0;
#if 1
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    std::complex<double> dot=0.0;
    for(int k=0;k<N;k++){
      dot += LdopDeflationInverseMatrix[N*i+k]* LdopDeflationMatrix[k*N+j];
    }
    std::complex<double> expect;
    if(i==j) expect=1.0;
    else     expect=0.0;

    double abserr=abs(expect - dot);

    if ( abserr > maxerr ) maxerr=abserr;

    if ( abs(expect - dot) > 1.0e-10 ) { 
      printf("Oops inverse test failed\n");
      exit(0);
    }
  }}
  linop_d->BossLog("LdopDeflationBasisInit: Max inverse error is %le\n",maxerr);

#endif

  DeflKrylovProj.resize(N);
  DeflKrylovMss.resize(N);
  DeflKrylovGsum.resize(N);

}

template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationBasisDiagonalise(int Nbasis)
{
  int No=  LdopDeflationBasisSize;

  std::vector<std::complex<cFloat> > mat  =LdopDeflationMatrix;
  std::vector<std::complex<cFloat> > evecs=LdopDeflationMatrix;
  std::vector<cFloat > evals(No);

  linop_d->BossLog("Calling out to LAPACK for matrix diagonalisation\n");
  if ( sizeof(cFloat)==sizeof(double) ) {
    LapackHermitianDiagonalise(No,(double *)&mat[0],(double *)&evecs[0],(double *)&evals[0]);
  } else { 
    LapackHermitianDiagonalise_f(No,(float *)&mat[0],(float *)&evecs[0],(float *)&evals[0]);
  }
  LdopDeflationEvals.resize(No);
  LdopDeflationInvEvals.resize(No);
  for(int i=0;i<No;i++){
    LdopDeflationEvals[i]=evals[i];
    LdopDeflationInvEvals[i]=1.0/evals[i];
  }

  // Check the diagonalisation
  if ( linop_d->isBoss() ) { 
    
    for(int k=0;k<No;k++){
      std::vector<std::complex<cFloat> > Mpsik(No);
      for(int i=0;i<No;i++){
	Mpsik[i]=0;
      }

      std::complex<cFloat> defect = 0;
      for(int i=0;i<No;i++){
	for(int j=0;j<No;j++){
	  Mpsik[i] += LdopDeflationMatrix[j*No+i]*evecs[j+k*No]; // Left evec, not right; right evec is adjoint
	}
      }

      for(int i=0;i<No;i++){
	Mpsik[i] = Mpsik[i] - evals[k]*evecs[i+k*No];
	defect+=conj(Mpsik[i])*Mpsik[i];
      }
      linop_d->NodeMessage("Evec[%d] : lambda %le : %le: M psi_k - lambda_k psi_k = %le,%le\n",k,evals[k],1.0/evals[k],real(defect),imag(defect)); 
    }
  }

  // Truncate the basis
  int N= Nbasis;

  std::vector< std::vector<std::complex<cFloat > > > RotBasis(N);
  for(int v=0;v<N;v++){
    RotBasis[v].resize(LocalNsubspace);
  }  

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      for(int i=0;i<N;i++){
	axpby(RotBasis[i],LdopDeflationBasis[0],LdopDeflationBasis[0],0.0,0.0);
	for(int j=0;j<No;j++){
	  std::complex<cFloat> c = conj(evecs[j+No*i]);
	  zaxpy(RotBasis[i],LdopDeflationBasis[j],RotBasis[i],c);
	}
      }
    }
  }

  LdopDeflationBasis = RotBasis;
  LdopDeflationAv    = RotBasis;
  LdopDeflationBasisSize = Nbasis;
  LdopDeflationMatrix.resize(N*N);
  LdopDeflationInverseMatrix.resize(N*N);
  DeflKrylovProj.resize(N);
  DeflKrylovMss.resize(N);
  DeflKrylovGsum.resize(N);

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      // Compute matrix elements
      for(int v=0;v<N;v++){
	linop_d->thread_barrier();
	ApplyThread(LdopDeflationBasis[v],LdopDeflationAv[v]);
	for(int u=0;u<N;u++){
	  std::complex<double> val          = innerProduct(LdopDeflationBasis[u],LdopDeflationAv[v]);
	  LdopDeflationMatrix[N*u+v]        =val;
	  LdopDeflationInverseMatrix[N*u+v] =val;
	}
      }
    }
  }

  if ( sizeof(cFloat) == sizeof(double) ) { 
    LapackHermitianInvert(N,(double *)&LdopDeflationInverseMatrix[0]);
  } else { 
    LapackHermitianInvert_f(N,(float *)&LdopDeflationInverseMatrix[0]);
  }
  for(int v=0;v<N;v++){
    LdopDeflationEvals[v]=real(LdopDeflationMatrix[N*v+v]);
  }  
  LdopDeflationIsDiagonal=1;
}



template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationProject(std::vector<std::complex<cFloat> > &localsubspace,
						std::vector<std::complex<cFloat> > &deflvec)
{
  if( localsubspace.size() != LocalNsubspace ) { 
    exit(0);
  }
  if( deflvec.size() != LdopDeflationBasisSize ) { 
    exit(0);
  }
  int nwork,thrlen,throff,me;

  nwork = LdopDeflationBasisSize;
  linop_d->thread_work(nwork,me,thrlen,throff);
  uint64_t t0=GetTimeBase();
  for(int v=throff;v<throff+thrlen;v++){
    float *fp=(float *)&deflvec[v];
    if ( sizeof(cFloat)==sizeof(double) ) { 
      qpx_inner((double *)&DeflKrylovGsum[v],LocalNsubspace,
		(double *)&LdopDeflationBasis[v][0],
		(double *)&localsubspace[0]);
    } else { 
      qpx_inner_f((double *)&DeflKrylovGsum[v],LocalNsubspace,
		  (float *)&LdopDeflationBasis[v][0],
		  (float *)&localsubspace[0]);
    }
  }
  linop_d->thread_barrier();
  uint64_t t1=GetTimeBase();
  linop_d->comm_gsum((double *)&DeflKrylovGsum[0],2*LdopDeflationBasisSize);

  for(int v=throff;v<throff+thrlen;v++){
    deflvec[v] = DeflKrylovGsum[v];
  }

  uint64_t t2=GetTimeBase();
  linop_d->thread_barrier();
  static int printed=0;
  if ( !printed ) { 
    linop_d->ThreadBossMessage("DeflationProject: cycle breakdown  %ld %ld\n",t1-t0,t2-t1);
    printed=1;
  }
}
template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationPromote(std::vector<std::complex<cFloat> >&deflvec,
					std::vector<std::complex<cFloat> >&localsubspace)
{
  if( localsubspace.size() != LocalNsubspace ) { 
    exit(0);
  }
  if( deflvec.size() != LdopDeflationBasisSize ) { 
    exit(0);
  }
  int me,thrlen,throff;
#ifdef USE_XLC_OPTIMISED_CODE
  // Could optimise further with BAGEL kernel, about 10% overhead in LdopM1, and perhaps 5% overall.
  // ~5s saving in 115s on 48^3 and 1k rack.
  linop_d->thread_work(LocalNsubspace,me,thrlen,throff);
  if ( sizeof(cFloat)==sizeof(double) ) {
    for(int v=0;v<LdopDeflationBasisSize;v++){
      if(v==0){
	qpx_axpby(thrlen,(double *)&localsubspace[throff],
		  (double *)&LdopDeflationBasis[v][throff],
		  (double *)&LdopDeflationBasis[v][throff],0.0,0.0);
      }
      qpx_zaxpy(thrlen,(double *)&localsubspace[throff],
		(double *)&LdopDeflationBasis[v][throff],
		(double *)&localsubspace[throff],
		(double *)&deflvec[v]);
    }

  } else { 
    for(int v=0;v<LdopDeflationBasisSize;v++){
      if(v==0){
	qpx_axpby_f(thrlen,
		    (float *)&localsubspace[throff],
		    (float *)&LdopDeflationBasis[v][throff],
		    (float *)&LdopDeflationBasis[v][throff],0.0,0.0);
      }
      std::complex<double> zz= deflvec[v];
      qpx_zaxpy_f(thrlen,
		  (float *)&localsubspace[throff],
		  (float *)&LdopDeflationBasis[v][throff],
		  (float *)&localsubspace[throff],
		  (double *)&zz);
    }
  }
  linop_d->thread_barrier();
#else
  zeroOut(localsubspace);
  for(int v=0;v<LdopDeflationBasisSize;v++){
    zaxpy(localsubspace,LdopDeflationBasis[v],localsubspace,deflvec[v]);
  }
#endif 
}
template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationMatPromote(std::vector<std::complex<cFloat> >&deflvec,
					   std::vector<std::complex<cFloat> >&localsubspace)
{
  if( localsubspace.size() != LocalNsubspace ) { 
    exit(0);
  }
  if( deflvec.size() != LdopDeflationBasisSize ) { 
    exit(0);
  }
  int me,thrlen,throff;
#ifdef USE_XLC_OPTIMISED_CODE
  linop_d->thread_work(LocalNsubspace,me,thrlen,throff);
  if ( sizeof(cFloat)==sizeof(double) ) { 
    for(int v=0;v<LdopDeflationBasisSize;v++){
      if(v==0){
	qpx_axpby(thrlen,
		  (double *)&localsubspace[throff],
		  (double *)&LdopDeflationAv[v][throff],
		  (double *)&LdopDeflationAv[v][throff],0.0,0.0);
      }
      qpx_zaxpy(thrlen,
		(double *)&localsubspace[throff],
		(double *)&LdopDeflationAv[v][throff],
		(double *)&localsubspace[throff],
		(double *)&deflvec[v]);
    }
  } else { 
    for(int v=0;v<LdopDeflationBasisSize;v++){
      if(v==0){
	qpx_axpby_f(thrlen,
		  (float *)&localsubspace[throff],
		  (float *)&LdopDeflationAv[v][throff],
		  (float *)&LdopDeflationAv[v][throff],0.0,0.0);
      }
      std::complex<double> zz = deflvec[v];
      qpx_zaxpy_f(thrlen,
		(float *)&localsubspace[throff],
		(float *)&LdopDeflationAv[v][throff],
		(float *)&localsubspace[throff],
		(double *)&zz);
    }
  }
  linop_d->thread_barrier();
#else
  zeroOut(localsubspace);
  for(int v=0;v<LdopDeflationBasisSize;v++){
    zaxpy(localsubspace,LdopDeflationAv[v],localsubspace,deflvec[v]);
  }
#endif 
}

template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationMatrixInverseMult(std::vector<std::complex<cFloat> > &in,std::vector<std::complex<cFloat> >&out)
{
  int me, thrlen, throff;
  if(in.size() != LdopDeflationBasisSize ) exit(0);
  if(out.size()!= LdopDeflationBasisSize ) exit(0);
  // Could thread these => 32x speed up
  // Single thread -- could multithread if dominant
  // Block decomposition of evecs would speed it.
#if 1
  if ( LdopDeflationIsDiagonal ) { 
    linop_d->thread_work(LdopDeflationBasisSize,me,thrlen,throff);
    cFloat Lambda;
    for(int i=throff;i<throff+thrlen;i++){
      Lambda = LdopDeflationInvEvals[i];
      out[i] = in[i]*Lambda;
    }
  } else { 
    double coeff=0.0;
    me = linop_d->thread_barrier();
    if (me==0){
      if ( sizeof(cFloat)==sizeof(double) ) { 
	myzgemv(LdopDeflationBasisSize,
		(double *)&LdopDeflationInverseMatrix[0],
		(double *)&in[0],coeff,
		(double *)&out[0]);
      } else { 
	mycgemv_f(LdopDeflationBasisSize,
		  (float *)&LdopDeflationInverseMatrix[0],
		  (float *)&in[0],coeff,
		  (float *)&out[0]);
      }
    }
  }
#else
  me = linop_d->thread_barrier();
  if (me==0){
    for(int i=0;i<LdopDeflationBasisSize;i++){
      out[i] = 0.0;
      int o = i*LdopDeflationBasisSize;
      for(int j=0;j<LdopDeflationBasisSize;j++){
	out[i] += LdopDeflationInverseMatrix[o+j]*in[j];
      }
    }
  }
#endif
  linop_d->thread_barrier();
}
template<class cFloat> 
void BfmHDCG<cFloat>::LdopDeflationMatrixMult(std::vector<std::complex<cFloat> > &in,std::vector<std::complex<cFloat> >&out)
{
  int me,thrlen,throff;
  if(in.size() != LdopDeflationBasisSize ) exit(0);
  if(out.size() != LdopDeflationBasisSize ) exit(0);
  // Single thread -- could multithread if dominant
#if 1
  if ( LdopDeflationIsDiagonal ) { 
    cFloat Lambda;
    linop_d->thread_work(LdopDeflationBasisSize,me,thrlen,throff);
    for(int i=throff;i<throff+thrlen;i++){
      Lambda = LdopDeflationEvals[i];
      out[i] = in[i]*Lambda; 
    }
  } else {
    me = linop_d->thread_barrier();
    if ( me == 0 ) { 
      double coeff=0.0;
      if ( sizeof(cFloat) == sizeof(double) ) {
	myzgemv(LdopDeflationBasisSize,
		(double *)&LdopDeflationMatrix[0],
		(double *)&in[0],coeff,
		(double *)&out[0]);
      } else { 
	mycgemv_f(LdopDeflationBasisSize,
		  (float *)&LdopDeflationMatrix[0],
		  (float *)&in[0],coeff,
		  (float *)&out[0]);
      }
    }
  }
#else
  me = linop_d->thread_barrier();
  if ( me == 0 ) { 
    for(int i=0;i<LdopDeflationBasisSize;i++){
      out[i] = 0.0;
      int o = i*LdopDeflationBasisSize;
      for(int j=0;j<LdopDeflationBasisSize;j++){
	out[i] += LdopDeflationMatrix[o+j]*in[j];
      }
    }
  }
#endif

  linop_d->thread_barrier();
}

template<class cFloat> 
template<class Float> 
void BfmHDCG<cFloat>::PolyMdagMprec(bfm_internal<Float> *lop,Fermion_t in,Fermion_t out,std::vector<double> &coeffs )
{

  Fermion_t Mtmp= lop->threadedAllocFermion();
  lop->precision_test=SloppyComms;
  Fermion_t tmp = lop->threadedAllocFermion();
  Fermion_t AtoN= lop->threadedAllocFermion();

  uint64_t t1=GetTimeBase();
  lop->axpby(AtoN,in,in,1.0,0.0);
  lop->axpby(out,AtoN,AtoN,coeffs[0],0.0);
  for(int i=1;i<coeffs.size();i++){
    lop->Mprec(AtoN,tmp,Mtmp,DaggerNo);  
    lop->Mprec(tmp,AtoN,Mtmp,DaggerYes);    
    lop->axpy(out,AtoN,out,coeffs[i]);
  }
  uint64_t t2=GetTimeBase();

  int order = coeffs.size();
  double flops   = lop->mprecFlops()*2.0*(order-1);
  flops =    flops   +  order* 24*lop->cbSites()*(2); // axpy

  double microseconds = 1.0*(t2-t1)/MHz();
  lop->ThreadBossMessage("Polynomial of order %d applied in %le s %le mflop/s\n",order,microseconds*1.0e-6,flops/microseconds);

#ifdef DEBUG_HDCG
  if ( 0 ) { 
    lop->Mprec(out,tmp,Mtmp,DaggerNo);
    lop->Mprec(tmp,AtoN,Mtmp,DaggerYes);
    double ns=lop->norm(in);
    double nn=lop->axpy_norm(AtoN,AtoN,in,-1.0);
    lop->ThreadBossMessage("Polynomial residual %le\n",sqrt(nn/ns));
  }
#endif
  lop->threadedFreeFermion(Mtmp);
  lop->threadedFreeFermion(tmp);
  lop->precision_test=0;
  lop->threadedFreeFermion(AtoN);
}

template<class cFloat> 
template<class Float> 
void BfmHDCG<cFloat>::PolyMdagMprec(bfm_internal<Float> *lop,Fermion_t in,Fermion_t out,Chebyshev & approx )
{
  double *coeffs = approx.Coeffs;
  double lo      = approx.lo;
  double hi      = approx.hi;
  int order      = approx.order;
  uint64_t t0,t1,t2;
  double n;
  int me = lop->thread_barrier();

  t0=GetTimeBase();
  Fermion_t Mtmp= lop->threadedAllocFermion();
  Fermion_t Tnm = lop->threadedAllocFermion();
  Fermion_t Tn  = lop->threadedAllocFermion();
  Fermion_t Tnp = lop->threadedAllocFermion();
  Fermion_t y   = lop->threadedAllocFermion();
  
  double xscale = 2.0/(hi-lo);
  double mscale = -(hi+lo)/(hi-lo);

  t1=GetTimeBase();
  Fermion_t T0=Tnm;
  Fermion_t T1=Tn;
  // Tnm=T0=in
  lop->axpy(T0,in,in,0.0);              

  // Tn=T1 = (xscale MdagM + mscale)in
  Fermion_t tmp=Tnp;;
  lop->Mprec(T0 ,tmp,Mtmp,DaggerNo);  
  lop->Mprec(tmp, y,Mtmp,DaggerYes);    
  lop->axpby(T1 , y,in,xscale,mscale);    
  
  // sum = .5c[0]T0 + c[1]T1
  lop->axpby(out,T0,T1,coeffs[0]/2.0,coeffs[1]); 

  for(int i=2;i<order;i++){

    Fermion_t tmp=Tnp;;

    lop->Mprec(Tn,tmp,Mtmp,DaggerNo);  
    lop->Mprec(tmp, y,Mtmp,DaggerYes);    

    lop->axpby(y,y,Tn,xscale,mscale);    // y Tn = [xscale MdagM+mscale] Tn
    lop->axpby(Tnp,y,Tnm,2.0,-1.0);      // Tnp=2yTn - Tnm

    lop->axpy(out,Tnp,out,coeffs[i]);    //Accumulate

    Fermion_t swizzle;
    swizzle=Tnm;
    Tnm    =Tn;
    Tn     =Tnp;
    Tnp    =swizzle;
  }
  t2=GetTimeBase();
#ifdef DEBUG_HDCG
  double flops   = lop->mprecFlops()*2.0*(order-1);
  flops =    flops   +  (order-2)* 24*lop->cbSites()*(3+3+2); // axpby,axpby,axpy
  flops =    flops   +             24*lop->cbSites()*(3+3+2); // axpby,axpby,axpy in entry

  double microseconds = 1.0*(t2-t1)/MHz();
  double wastemicroseconds = 1.0*(t1-t0)/MHz();
  lop->ThreadBossMessage("Chebyshev Polynomial of order %d applied in %le s %le mflop/s\n",order,microseconds*1.0e-6,flops/microseconds);

  if ( 0 ) { 
    lop->Mprec(out,tmp,Mtmp,DaggerNo);
    lop->Mprec(tmp,y,Mtmp,DaggerYes);
    double ns=lop->norm(in);
    double nn=lop->axpy_norm(y,y,in,-1.0);
    lop->ThreadBossMessage("Chebyshev Polynomial residual %le\n",sqrt(nn/ns));
  }
#endif
  lop->threadedFreeFermion(y);
  lop->threadedFreeFermion(Mtmp);
  lop->threadedFreeFermion(Tnm);
  lop->threadedFreeFermion(Tn);
  lop->threadedFreeFermion(Tnp);


};
template<class cFloat>
template<class datum>
void BfmHDCG<cFloat>::ThreadedFreeVector(std::vector<datum> *vv)
{
  int me = linop_d->thread_barrier();
  if (!me ) { 
    delete vv;
  }
  linop_d->thread_barrier();
  return;
}

template<class cFloat>
template<class datum> 
std::vector<datum> *BfmHDCG<cFloat>::ThreadedAllocVector(int nvec)
{
  typedef std::vector<datum> * F_t ;
  int me = linop_d->thread_barrier();
  F_t tmp;
  if ( ! me ) {
    tmp = new std::vector<datum>; 
    tmp->resize(nvec);
  }
  tmp = (F_t)linop_d->thread_bcast(me,(void *)tmp);
  return tmp;
}


//
// Polynomial of the Ldop
//
template<class cFloat> 
void BfmHDCG<cFloat>::PolyMdagMprec(std::vector<std::complex<cFloat> > &in,
				    std::vector<std::complex<cFloat> > &out,
				    std::vector<double> &coeffs )
{
  typedef std::vector<std::complex<cFloat> > *F_t ;

  int me = linop_d->thread_barrier();

  int len = in.size();
  int depth=1;

  F_t tmp = ThreadedAllocVector<std::complex<cFloat> >(len);
  F_t AtoN= ThreadedAllocVector<std::complex<cFloat> >(len);

  uint64_t t1=GetTimeBase();
  axpby(*AtoN,in,in,1.0,0.0);
  axpby(out,*AtoN,*AtoN,coeffs[0],0.0);
  for(int i=1;i<coeffs.size();i++){
    F_t swizzle=tmp;
    tmp = AtoN;
    AtoN= swizzle;

    ApplyThreadOpt(*tmp,*AtoN,1,depth);
    axpy(out,*AtoN,out,coeffs[i]);
  }
  uint64_t t2=GetTimeBase();

#ifdef DEBUG_HDCG
  int order = coeffs.size();
  double microseconds = 1.0*(t2-t1)/MHz();
  linop_d->ThreadBossMessage("Ldop Polynomial<%d> order %d applied in %le s \n",depth,order,microseconds*1.0e-6);

  static int print;
  if ( print < 10 ) { 
    linop_d->ThreadBossMessage("LdopPolynomial of order %d depth %d applied in %le s\n",order,depth,microseconds*1.0e-6);

    ApplyThreadOpt(out,*tmp,1,depth); // reduced depth halo
    axpy(*tmp,*tmp,in,-1.0);
    double ns=norm_vec(in);
    double nn=norm_vec(*tmp);
    linop_d->ThreadBossMessage("LdopPolynomial residual(1) %le\n",sqrt(nn/ns));

    ApplyThreadOpt(out,*tmp,1,4);// Full depth halo
    axpy(*tmp,*tmp,in,-1.0);
    ns=norm_vec(in);
    nn=norm_vec(*tmp);
    linop_d->ThreadBossMessage("LdopPolynomia(l residual(4) %le\n",sqrt(nn/ns));

    if ( !me ) print ++;
  }
#endif
  linop_d->thread_barrier();

  ThreadedFreeVector<std::complex<cFloat> >(AtoN);
  ThreadedFreeVector<std::complex<cFloat> >(tmp);

}

template<class cFloat> 
void BfmHDCG<cFloat>::PolyMdagMprec(std::vector<std::complex<cFloat> > &in,
				    std::vector<std::complex<cFloat> > &out,
				    Chebyshev & approx )
{
  typedef std::vector<std::complex<cFloat> > *F_t ;

  double *coeffs = approx.Coeffs;
  double lo      = approx.lo;
  double hi      = approx.hi;
  int order      = approx.order;
  int depth = 1;
  double n;
  uint64_t t0,t1,t2;

  int me = linop_d->thread_barrier();
  t0=GetTimeBase();
  
  int len = in.size();
  F_t Mtmp= ThreadedAllocVector<std::complex<cFloat> >(len);
  F_t Tnm = ThreadedAllocVector<std::complex<cFloat> >(len);
  F_t Tn  = ThreadedAllocVector<std::complex<cFloat> >(len);
  F_t Tnp = ThreadedAllocVector<std::complex<cFloat> >(len);
  F_t y   = ThreadedAllocVector<std::complex<cFloat> >(len);
  
  double xscale = 2.0/(hi-lo);
  double mscale = -(hi+lo)/(hi-lo);

  t1=GetTimeBase();

  F_t T0=Tnm;
  F_t T1=Tn;
  // Tnm=T0=in
  axpy(*T0,in,in,0.0);              

  // Tn=T1 = (xscale MdagM + mscale)in

  ApplyThreadOpt(*T0,*y,1,depth); 
  axpby(*T1 , *y,in,xscale,mscale);    
  
  // sum = .5c[0]T0 + c[1]T1
  axpby(out,*T0,*T1,coeffs[0]/2.0,coeffs[1]); 

  for(int i=2;i<order;i++){

    F_t tmp=Tnp;

    ApplyThreadOpt(*Tn,*y,1,depth);

    axpby(*y,*y,*Tn,xscale,mscale);    // y Tn = [xscale MdagM+mscale] Tn
    axpby(*Tnp,*y,*Tnm,2.0,-1.0);      // Tnp=2yTn - Tnm
    axpy(out,*Tnp,out,coeffs[i])  ;    //Accumulate

    F_t swizzle;
    swizzle=Tnm;
    Tnm    =Tn;
    Tn     =Tnp;
    Tnp    =swizzle;
  }
  t2=GetTimeBase();

#ifdef DEBUG_HDCG
  double microseconds = 1.0*(t2-t1)/MHz();
  double wastemicroseconds = 1.0*(t1-t0)/MHz();

  static int print;
  if ( print < 10 ) { 
    linop_d->ThreadBossMessage("Chebyshev LdopPolynomial of order %d depth %d applied in %le s\n",order,depth,microseconds*1.0e-6);

    ApplyThreadOpt(out,*y,1,depth); // reduced depth halo
    axpy(*y,*y,in,-1.0);
    double ns=norm_vec(in);
    double nn=norm_vec(*y);
    linop_d->ThreadBossMessage("Chebyshev LdopPolynomial residual(1) %le\n",sqrt(nn/ns));

    ApplyThreadOpt(out,*y,1,4);// Full depth halo
    axpy(*y,*y,in,-1.0);
    ns=norm_vec(in);
    nn=norm_vec(*y);
    linop_d->ThreadBossMessage("Chebyshev LdopPolynomia(l residual(4) %le\n",sqrt(nn/ns));

    if ( !me ) print ++;
  }
#endif
  linop_d->thread_barrier();

  ThreadedFreeVector<std::complex<cFloat> >(y);
  ThreadedFreeVector<std::complex<cFloat> >(Mtmp);
  ThreadedFreeVector<std::complex<cFloat> >(Tnm);
  ThreadedFreeVector<std::complex<cFloat> >(Tnp);
  ThreadedFreeVector<std::complex<cFloat> >(Tn);
};


  /*
   * Compared to Tang-2009:  P=Pleft. P^T = PRight Q=MssInv. 
   * Script A = SolverMatrix 
   * Script P = Preconditioner
   *
   * Deflation methods considered
   *      -- Solve P A x = P b        [ like Luscher ]
   * DEF-1        M P A x = M P b     [i.e. left precon]
   * DEF-2        P^T M A x = P^T M b
   * ADEF-1       Preconditioner = M P + Q      [ Q + M + M A Q]
   * ADEF-2       Preconditioner = P^T M + Q
   * BNN          Preconditioner = P^T M P + Q
   * BNN2         Preconditioner = M P + P^TM +Q - M P A M 
   * 
   * Implement ADEF-2
   *
   * Vstart = P^Tx + Qb
   * M1 = P^TM + Q
   * M2=M3=1
   * Vout = x
   */

template<class cFloat> 
template<class Float> 
void BfmHDCG<cFloat>::PcgMlinop(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_internal<Float> *lop,double shift,std::vector<double>&poly)
{
  int me = lop->thread_barrier();
  double restore_rsd  = lop->residual;
  double restore_iter=  lop->max_iter;

  lop->thread_barrier();
  lop->residual = PreconditionerKrylovResidual;
  lop->max_iter = PreconditionerKrylovIterMax;
  lop->precision_test=SloppyComms;
  lop->thread_barrier();

  uint64_t t1 = GetTimeBase();
  lop->CGNE_single_shift(out,in,shift,poly); // Records the CG polynomial

  if (! me ) { 
    PcgMprec+=GetTimeBase()-t1;
  }
  
  lop->thread_barrier();
  lop->precision_test=0;
  lop->residual=restore_rsd;
  lop->max_iter=restore_iter;
  lop->thread_barrier();
  return;
}

template<class cFloat> 
void BfmHDCG<cFloat>::PcgM_f(Fermion_t in,Fermion_t out,Fermion_t tmp)
{
  PcgMlinop<float>(in,out,tmp,linop_f,PreconditionerKrylovShift,CGpolynomial);
  return;
}
template<class cFloat> 
void BfmHDCG<cFloat>::PcgM(Fermion_t in,Fermion_t out,Fermion_t tmp)
{
  PcgMlinop<double>(in,out,tmp,linop_d,PreconditionerKrylovShift,CGpolynomial);
  return;
}

template<class cFloat> 
void BfmHDCG<cFloat>::PcgM1(Fermion_t in, Fermion_t out,Fermion_t tmp,Fermion_t mp,int Smoother)
{
  uint64_t t1;
  int me = linop_d->thread_barrier();

  Fermion_t mgtmp;
  Fermion_t in_f;
  Fermion_t Mtmp= linop_d->threadedAllocFermion();
  Fermion_t Min = linop_d->threadedAllocFermion();
  int block;

  switch(PcgType) { 

  case PcgPrec:
  case PcgDef1:
  case PcgDef2:
    //M
    PcgM(in,out,tmp);

    break;
  case PcgMssDef:

    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop_d->axpy(out,in,tmp,1.0);
    break;

  case PcgAD:
    // M+Q
    PcgM(in,out,tmp);

    t1=GetTimeBase();
    ProjectToSubspace(in,PleftProj);     
    
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  


    linop_d->axpy(out,out,tmp,1.0);
    break;

  case PcgADef1:
    //MP  + Q
    Pleft(in,mp);
    PcgM(mp,out,tmp);
    
    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop_d->axpy(out,out,tmp,1.0);
    break;

  case PcgAdef2fSingleShift:

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    // A and Q are shifted by PcgSingleShift

    t1=GetTimeBase();
    linop_d->thread_barrier();
    if ( linop_f->SPIcomms() && !me )  linop_f->comm_init();
    linop_d->thread_barrier();
    in_f =linop_f->threadedAllocFermion();
    linop_f->precisionChange(in,in_f,DoubleToSingle,1);
    if ( !me ) PcgPrecChange+=GetTimeBase()-t1;

    PcgM_f(in_f,Min,tmp); // Shifted solver as preconditioner
    t1= GetTimeBase();
    linop_f->Mprec(Min,tmp,Mtmp,DaggerNo);
    linop_f->Mprec(tmp,out,Mtmp,DaggerYes);  // out  = A Min
    linop_f->axpy(out,Min,out,PcgSingleShift);
    linop_f->axpy(tmp,out,in_f,-1.0);        // tmp  = in - A Min
    if ( !me ) PcgMprec+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ProjectToSubspace_f(tmp,PleftProj);     
    if ( !me ) PcgProj+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ApplyInverseCG(PleftProj,PleftMss_proj,PcgSingleShift); // Ass^{-1} [in - A Min]_s
    if ( !me ) PcgMssInv+=GetTimeBase()-t1;

    t1=GetTimeBase();
    PromoteFromSubspace_f(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    if ( !me ) PcgProm+=GetTimeBase()-t1;

    t1=GetTimeBase();
    linop_f->axpy(Mtmp,Min,tmp,1.0); // Min+tmp
    if ( !me ) PcgLinalg+=GetTimeBase()-t1;


    t1=GetTimeBase();
    linop_f->precisionChange(Mtmp,out,SingleToDouble,1);
    linop_d->thread_barrier();
    if ( linop_d->SPIcomms() && !me )  linop_d->comm_init();
    linop_d->thread_barrier();
  
    linop_f->threadedFreeFermion(in_f);
    if ( !me ) PcgPrecChange+=GetTimeBase()-t1;

    break;
  case PcgV11f:
    t1=GetTimeBase();
    linop_d->thread_barrier();
    if ( linop_f->SPIcomms() && !me )  linop_f->comm_init();
    linop_d->thread_barrier();
    in_f =linop_f->threadedAllocFermion();
    mgtmp =linop_f->threadedAllocFermion();
    
    linop_f->precisionChange(in,in_f,DoubleToSingle,1);

    if ( !me ) PcgPrecChange+=GetTimeBase()-t1;
    
    ////////////////////////////////
    // M_{IRS} preconditioner
    ////////////////////////////////
    if ( Smoother == SmootherMirs ) { 
      PcgMlinop<float>(in_f,Min,tmp,linop_f,PreconditionerKrylovShift,CGpolynomial);
    } else if ( Smoother==SmootherMirsPoly) { 
      uint64_t t0=GetTimeBase();
      PolyMdagMprec<float>(linop_f,in_f,Min,CGpolynomial);
      if(!me) PcgMprec+=GetTimeBase()-t0;
    } else { 
      exit(-1);
    }

    t1= GetTimeBase();
    linop_f->Mprec(Min,tmp,Mtmp,DaggerNo);
    linop_f->Mprec(tmp,out,Mtmp,DaggerYes);  // out  = A Min
    linop_f->axpy(tmp,out,in_f,-1.0);        // tmp  = in - A Min
    if ( !me ) PcgMprec+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ProjectToSubspace_f(tmp,PleftProj);     
    if ( !me ) PcgProj+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    if ( !me ) PcgMssInv+=GetTimeBase()-t1;

    t1=GetTimeBase();
    PromoteFromSubspace_f(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    if ( !me ) PcgProm+=GetTimeBase()-t1;

    t1=GetTimeBase();
    linop_f->axpy(Mtmp,Min,tmp,1.0); // Min+tmp
    if ( !me ) PcgLinalg+=GetTimeBase()-t1;

    // Three preconditioner smoothing -- hermitian if C3 = C1
    // Compute error
    linop_f->Mprec(Mtmp,mgtmp,tmp,DaggerNo);
    linop_f->Mprec(mgtmp,out,tmp,DaggerYes);
    linop_f->axpy(out,out,in_f,-1.0);

    // Reapply smoother
    if ( Smoother==SmootherMirsPoly) { 
      uint64_t t0=GetTimeBase();
      PolyMdagMprec<float>(linop_f,out,Min,SparePolynomial);
      if(!me) PcgMprec+=GetTimeBase()-t0;
    } else if ( Smoother == SmootherMirs ) { 
      PcgMlinop<float>(out,Min,tmp,linop_f,PreconditionerKrylovShift,SparePolynomial);
    } else { 
      exit(-1);
    }
    linop_f->axpy(Mtmp,Min,Mtmp,1.0); // Mtmp+Min

    t1=GetTimeBase();
    linop_f->precisionChange(Mtmp,out,SingleToDouble,1);
    linop_d->thread_barrier();
    if ( linop_d->SPIcomms() && !me )  linop_d->comm_init();
    linop_d->thread_barrier();
  
    linop_f->threadedFreeFermion(in_f);
    linop_f->threadedFreeFermion(mgtmp);

    if ( !me ) PcgPrecChange+=GetTimeBase()-t1;

    break;

    break;
  case PcgAdef2f:
    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]

    ///////////////////////////////////////////
    // Move the project to single into here
    // New matrix
    // ADEF2f -- single precision
    // - Use single precision linop
    // - new routine: PcgM_float
    // - new routine: ProjectToSubspace_float
    // - new routine: PromoteFromSubspace_float
    // 
    ///////////////////////////////////////////
    AnalyseSpectrumInternal = 0;

    t1=GetTimeBase();
    linop_d->thread_barrier();
    if ( linop_f->SPIcomms() && !me )  linop_f->comm_init();
    linop_d->thread_barrier();
    in_f =linop_f->threadedAllocFermion();
    mgtmp =linop_f->threadedAllocFermion();
    
    
    linop_f->precisionChange(in,in_f,DoubleToSingle,1);

    if ( !me ) PcgPrecChange+=GetTimeBase()-t1;
    
    if ( Smoother != SmootherNone
	 &&Smoother != SmootherChebyshev
	 &&Smoother != SmootherMirs
	 &&Smoother != SmootherMirsPoly ) { 
      exit(-1);
    }
    if ( AnalyseSpectrumInternal ) {
      SpectralDecomposition<float>(in_f,linop_f,"in_f");
    }

    ////////////////////////////////
    // No preconditioning
    ////////////////////////////////
    if ( Smoother==SmootherNone) { 
      linop_f->axpy(Min,in_f,in_f,0.0);
    } 

    ////////////////////////////////
    // Chebyshev preconditioner
    ////////////////////////////////
    if ( Smoother==SmootherChebyshev ) { 
      Chebyshev Approx;  Approx.Init(0.3,30.0,7,PolynomialShape);
      uint64_t t0=GetTimeBase();
      PolyMdagMprec<float>(linop_f,in_f,Min,Approx);
      if(!me) PcgMprec+=GetTimeBase()-t0;
      Approx.End();
    } 
    
    ////////////////////////////////
    // M_{IRS} preconditioner
    ////////////////////////////////
    if ( Smoother==SmootherMirs) {   
      PcgMlinop<float>(in_f,Min,tmp,linop_f,PreconditionerKrylovShift,CGpolynomial);
    } 

    ////////////////////////////////
    // M_{IRS} preconditioner recorded polynomial
    ////////////////////////////////
    if ( Smoother==SmootherMirsPoly) { 
      uint64_t t0=GetTimeBase();
      PolyMdagMprec<float>(linop_f,in_f,Min,CGpolynomial);
      if(!me) PcgMprec+=GetTimeBase()-t0;

      if ( AnalyseSpectrumInternal ) {
	SpectralDecomposition<float>(Min,linop_f," Mirs in");
      }
    }

    t1= GetTimeBase();
    linop_f->Mprec(Min,tmp,Mtmp,DaggerNo);
    linop_f->Mprec(tmp,out,Mtmp,DaggerYes);  // out  = A Min
    linop_f->axpy(tmp,out,in_f,-1.0);        // tmp  = in - A Min
    if ( !me ) PcgMprec+=GetTimeBase()-t1;
    if ( AnalyseSpectrumInternal ) {
      SpectralDecomposition<float>(tmp,linop_f,"[1- A Mirs] in");
    }

    if ( AnalyseSpectrumInternal ) {
      linop_f->axpy(out,Min,tmp,-PreconditionerKrylovShift);        // tmp  = in - (A+shift) Min
      SpectralDecomposition<float>(out,linop_f,"[1- (A+shift) Mirs] in");
    }
    

    t1=GetTimeBase();
    ProjectToSubspace_f(tmp,PleftProj);     
    if ( !me ) PcgProj+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    if ( !me ) PcgMssInv+=GetTimeBase()-t1;

    t1=GetTimeBase();
    PromoteFromSubspace_f(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    if ( !me ) PcgProm+=GetTimeBase()-t1;

    if ( AnalyseSpectrumInternal ) {
      SpectralDecomposition<float>(tmp,linop_f," Q [1 - A Min] r");
    }

    t1=GetTimeBase();
    linop_f->axpy(Mtmp,Min,tmp,1.0); // Min+tmp
    if ( !me ) PcgLinalg+=GetTimeBase()-t1;


    if ( AnalyseSpectrumInternal ) {
      SpectralDecomposition<float>(Mtmp,linop_f,"new search z");
    }

    if ( AnalyseSpectrumInternal ) {
      linop_f->Mprec(Mtmp,mgtmp,tmp,DaggerNo);
      linop_f->Mprec(mgtmp,out,tmp,DaggerYes);
      SpectralDecomposition<float>(out,linop_f,"A z");
      linop_f->axpy(out,out,in_f,-1.0);
      SpectralDecomposition<float>(out,linop_f,"A z - r");
    }

    t1=GetTimeBase();
    linop_f->precisionChange(Mtmp,out,SingleToDouble,1);
    linop_d->thread_barrier();
    if ( linop_d->SPIcomms() && !me )  linop_d->comm_init();
    linop_d->thread_barrier();
  
    linop_f->threadedFreeFermion(in_f);
    linop_f->threadedFreeFermion(mgtmp);

    if ( !me ) PcgPrecChange+=GetTimeBase()-t1;

    break;
  case PcgAdef2:

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]


    PcgM(in,Min,tmp); // Shifted solver
    t1= GetTimeBase();
    linop_d->Mprec(Min,tmp,Mtmp,DaggerNo);
    linop_d->Mprec(tmp,out,Mtmp,DaggerYes);  // out  = A Min
    linop_d->axpy(tmp,out,in,-1.0);          // tmp  = in - A Min
    if ( !me ) PcgMprecO+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ProjectToSubspace(tmp,PleftProj);     
    if ( !me ) PcgProj+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    if ( !me ) PcgMssInv+=GetTimeBase()-t1;

    t1=GetTimeBase();
    PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    if ( !me ) PcgProm+=GetTimeBase()-t1;

    linop_d->axpy(out,Min,tmp,1.0); // Min+tmp
    

    // PT M + Q
    /*
    PcgM(in,tmp,out);
    Pright(tmp,out);

    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop_d->axpy(out,out,tmp,1.0);
    */
    break;

  case PcgBNN:
    // PT M P +Q 
    Pleft(in,out);
    PcgM(out,mp,tmp);
    Pright(mp,out);
    
    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop_d->axpy(out,out,tmp,1.0);
    break;
  default :
    exit(0);
  }
  linop_d->threadedFreeFermion(Mtmp);
  linop_d->threadedFreeFermion(Min);
}
template<class cFloat> 
void BfmHDCG<cFloat>::PcgM2(Fermion_t in, Fermion_t out)
{
  int me = linop_d->thread_barrier();
  uint64_t t1;
  switch(PcgType) { 
  case PcgPrec:
  case PcgAD:
  case PcgDef1:
  case PcgADef1:
  case PcgAdef2:
  case PcgAdef2f:
  case PcgV11f:
  case PcgAdef2fSingleShift:
  case PcgBNN:
  case PcgMssDef:
    t1=GetTimeBase();
    linop_d->axpy(out,in,in,0.0);
    if (!me) PcgLinalg+=GetTimeBase()-t1;
    break;
  case PcgDef2:
    Pright(in,out);
    break;
  default :
    exit(0);
  }
}

template<class cFloat> 
double BfmHDCG<cFloat>::PcgM3(Fermion_t p, Fermion_t mp,Fermion_t mmp, Fermion_t tmp)
{
  double d;
  uint64_t t1;
 
  int me = linop_d->thread_barrier();

  switch(PcgType) { 
  case PcgPrec:
  case PcgAD:
  case PcgDef2:
  case PcgADef1:
  case PcgAdef2:
  case PcgAdef2f:
  case PcgV11f:
  case PcgBNN:
  case PcgMssDef:
    t1=GetTimeBase();
    d=linop_d->Mprec(p,mp,tmp,0,1);// Dag no
      linop_d->Mprec(mp,mmp,tmp,1);// Dag yes
    if(!me) PcgMprecO+=GetTimeBase()-t1;
      
    break;
  case PcgAdef2fSingleShift:
    t1=GetTimeBase();
    linop_d->Mprec(p,mp,tmp,0,1);// Dag no
    linop_d->Mprec(mp,mmp,tmp,1);// Dag yes
    linop_d->axpy(mmp,p,mmp,PcgSingleShift);
    d=real(linop_d->inner(p,mmp));
    if(!me) PcgMprecO+=GetTimeBase()-t1;
      
    break;

  case PcgDef1:
    d=linop_d->Mprec(p,mmp,tmp,0,1);// Dag no
      linop_d->Mprec(mmp,mp,tmp,1);// Dag yes
    Pleft(mp,mmp);
    d=real(linop_d->inner(p,mmp));
    break;
  default :
    exit(0);
  }
  return d;
}

template<class cFloat> 
void BfmHDCG<cFloat>::PcgVstart(Fermion_t x, Fermion_t src,Fermion_t r,
			        Fermion_t mp, Fermion_t mmp, Fermion_t tmp)
{
  uint64_t t1;
  double nn;

  LittleDopSolverResidual=LittleDopSolverResidualVstart;

  int me = linop_d->thread_barrier();
  switch(PcgType) { 
  case PcgPrec:
  case PcgMssDef:
  case PcgAD:
  case PcgDef1:
  case PcgADef1:
  case PcgBNN:

    break;

  case PcgDef2:
  case PcgAdef2:
  case PcgAdef2f:
  case PcgV11f:
    ///////////////////////////////////
    // Choose x_0 such that 
    // x_0 = guess +  (A_ss^inv) r_s = guess + Ass_inv [src -Aguess]
    //                               = [1 - Ass_inv A] Guess + Assinv src
    //                               = P^T guess + Assinv src 
    //                               = Vstart  [Tang notation]
    // This gives:
    // W^T (src - A x_0) = src_s - A guess_s - r_s
    //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
    //                   = 0 
    ///////////////////////////////////
    linop_d->Mprec(x,mp,tmp,DaggerNo);
    linop_d->Mprec(mp,mmp,tmp,DaggerYes);

    linop_d->axpy (r, mmp, src,-1.0);        // r_{-1} = src - A x

    ProjectToSubspace(r,PleftProj);     
    ApplyInverseCG(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,mp);  

    linop_d->axpy(x,x,mp,1.0);
    break;
  case PcgAdef2fSingleShift:
    ///////////////////////////////////
    // Choose x_0 such that 
    // x_0 = guess +  (A_ss^inv) r_s = guess + Ass_inv [src -Aguess]
    //                               = [1 - Ass_inv A] Guess + Assinv src
    //                               = P^T guess + Assinv src 
    //                               = Vstart  [Tang notation]
    // This gives:
    // W^T (src - A x_0) = src_s - A guess_s - r_s
    //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
    //                   = 0 
    ///////////////////////////////////
    linop_d->Mprec(x,mp,tmp,DaggerNo);
    linop_d->Mprec(mp,mmp,tmp,DaggerYes);
    linop_d->axpy(mmp,x,mmp,PcgSingleShift);
    linop_d->axpy (r, mmp, src,-1.0);        // r_{-1} = src - A x

    nn = linop_d->norm(r);
    linop_d->ThreadBossDebug("Shifted Vstart r %le\n",nn);
    
    ProjectToSubspace(r,PleftProj);     
    ApplyInverseCG(PleftProj,PleftMss_proj,PcgSingleShift); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,mp);  

    nn = linop_d->norm(mp);
    linop_d->ThreadBossDebug("Shifted Vstart subspace %le\n",nn);

    linop_d->axpy(x,x,mp,1.0);
    break;
  default :
    exit(0);
  }
}

template<class cFloat> 
void BfmHDCG<cFloat>::PcgVout  (Fermion_t in, Fermion_t out,Fermion_t src,Fermion_t tmp)
{
  switch(PcgType) { 
  case PcgDef1:
    //Qb + PT x
    ProjectToSubspace(src,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  
    
    Pright(in,out);
    
    linop_d->axpy(out,tmp,out,1.0);
    break;

  case PcgPrec:
  case PcgAD:
  case PcgMssDef:
  case PcgDef2:
  case PcgADef1:
  case PcgAdef2:
  case PcgAdef2f:
  case PcgV11f:
  case PcgAdef2fSingleShift:
  case PcgBNN:
     linop_d->axpy(out,in,in,0.0);
    break;
  default :
    exit(0);
  }
}
template<class cFloat> 
void BfmHDCG<cFloat>::PcgReport (void)
{
  int me = linop_d->thread_barrier();

  LittleDopSolverResidual=LittleDopSolverResidualInner;
  // Timing in cycles
  PcgMssInv=0; 
  PcgMprecO=0;
  PcgMprec=0;
  PcgProj=0;
  PcgProm=0;
  PcgLinalg=0;
  PcgPrecChange=0;

  linop_d->ThreadBossLog("HDCG: ******************** PARAMS ***************\n");
  switch(PcgType){
  case PcgDef1:
    linop_d->ThreadBossLog("HDCG: algorithm DEF1\n");
    break;
  case PcgPrec:
    linop_d->ThreadBossLog("HDCG: algorithm Preconditioned CG\n");
    break;
  case PcgAD:
    linop_d->ThreadBossLog("HDCG: algorithm AD\n");
    break;
  case PcgDef2:
    linop_d->ThreadBossLog("HDCG: algorithm DEF2\n");
    break;
  case PcgADef1:
    linop_d->ThreadBossLog("HDCG: algorithm A-DEF1\n");
    break;
  case PcgAdef2:
    linop_d->ThreadBossLog("HDCG: algorithm A-DEF2\n");
    break;
  case PcgAdef2f:
    linop_d->ThreadBossLog("HDCG: algorithm A-DEF2(single prec accel)\n");
    break;
  case PcgV11f:
    linop_d->ThreadBossLog("HDCG: algorithm V(1,1) multigrid(single prec accel)\n");
    break;
  case PcgAdef2fSingleShift:
    linop_d->ThreadBossLog("HDCG: shifted algorithm A-DEF2(single prec accel)\n");
    break;
  case PcgBNN:
    linop_d->ThreadBossLog("HDCG: algorithm BNN\n");
    break;
  case PcgMssDef:
    linop_d->ThreadBossLog("HDCG: algorithm MssDef [equiv AD with M=1]\n");
      break;
  default:
    exit(0);
  }
  linop_d->ThreadBossLog("HDCG: Ldop inner solver residual %le \n",LittleDopSolverResidualInner);
  linop_d->ThreadBossLog("HDCG: PreconditionerKrylovResidual %le \n", PreconditionerKrylovResidual);
  linop_d->ThreadBossLog("HDCG: PreconditionerKrylovIterMax    %d \n",PreconditionerKrylovIterMax);
  linop_d->ThreadBossLog("HDCG: Preconditioner shift %le \n",PreconditionerKrylovShift);

}

template<class cFloat> 
int BfmHDCG<cFloat>::Pcg(BfmHDCG<double> * guess_solver,Fermion_t psi, Fermion_t src,Fermion_t resid)
{
  double f;
  double rtzp,rtz,a,d,b;
  double rptzp;

  int me = linop_d->thread_barrier();

  if ( linop_d->SPIcomms() && !me )  linop_d->comm_init();

  PcgReport();

  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterEnter();
  }

  Fermion_t x   = psi;
  Fermion_t p   = linop_d->threadedAllocFermion(); 
  Fermion_t z   = linop_d->threadedAllocFermion(); 
  Fermion_t tmp = linop_d->threadedAllocFermion(); 
  Fermion_t mp  = linop_d->threadedAllocFermion(); 
  Fermion_t mmp = linop_d->threadedAllocFermion(); 
  Fermion_t r   = linop_d->threadedAllocFermion(); 
  Fermion_t mu   = linop_d->threadedAllocFermion(); 
  Fermion_t rp = linop_d->threadedAllocFermion(); 

  //Initial residual computation & set up
  double guess = linop_d->norm(psi);
  double tn;

  //////////////////////////
  // x0 = Vstart -- possibly modify guess
  //////////////////////////
  PcgVstart(x,src,r,mp,mmp,tmp);

  // r0 = b -A x0
  linop_d->Mprec(x,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,mmp,tmp,DaggerYes);
  if ( PcgType == PcgAdef2fSingleShift ) { 
    linop_d->axpy(mmp,x,mmp,PcgSingleShift);
  }

  linop_d->axpy (r, mmp, src,-1.0);    // Recomputes r=src-x0
  linop_d->axpy(rp,r,r,0.0);

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  linop_d->ThreadBossMessage("HDCG: subspace residual for x0 = %le\n",nv);

  //////////////////////////////////
  // Compute z = M1 x
  //////////////////////////////////
  PcgM1(r,z,tmp,mp,SmootherMirs);
  rtzp =real(linop_d->inner(r,z));

  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  // Def2: p = 1 - Q Az = Pright z 
  // Other algos M2 is trivial
  ///////////////////////////////////////
  PcgM2(z,p);

  double ssq =  linop_d->norm(src);
  double rsq =  linop_d->residual* linop_d->residual*ssq;

  linop_d->ThreadBossLog("HDCG: k=0 residual %le rsq %le\n",rtzp,rsq);

  uint64_t t_start=GetTimeBase();
  uint64_t t1;
  LittleDopSolverResidual=LittleDopSolverResidualInner;
  linop_d->thread_barrier();

  for (int k=1;k<=linop_d->max_iter;k++){
    

    if ( k==33 && AnalyseSpectrum) { 
      AnalyseSpectrumInternal=1;
    } else { 
      AnalyseSpectrumInternal=0;
    }

    rtz=rtzp;
    d= PcgM3(p,mp,mmp,tmp);
    a = rtz/d;

    t1 = GetTimeBase();
    linop_d->axpy(x,p,x,a);
    double rn = linop_d->axpy_norm(r,mmp,r,-a);
    if ( !me )     PcgLinalg+=GetTimeBase()-t1;

    // Compute z = M x
    if ( ParamSmoother == SmootherMirsPoly ) {
      if ( (k < 5)  ||  ((k%10) == 0) ) { 
	PcgM1(r,z,tmp,mp,SmootherMirs);        // Record a polynomial 1,2,3,4, 10, 20, 30,...
      } else { 
	PcgM1(r,z,tmp,mp,SmootherMirsPoly);    // Replay a polynomial elsewise
      }
    } else {  // Freeze out the polynomials
      PcgM1(r,z,tmp,mp,ParamSmoother);
    }

    t1 = GetTimeBase();
    rtzp =real(linop_d->inner(r,z));

    int ipcg=1; // almost free inexact preconditioned CG
    if (ipcg) {
      rptzp =real(linop_d->inner(rp,z));
    } else {
      rptzp =0;
    }
    b = (rtzp-rptzp)/rtz;

    if ( !me )     PcgLinalg+=GetTimeBase()-t1;

    PcgM2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

    t1 = GetTimeBase();
    linop_d->axpy(p,p,mu,b);  // mu = A r
    if ( !me )     PcgLinalg+=GetTimeBase()-t1;

    double rrn=sqrt(rn/ssq);
    double rtn=sqrt(rtz/ssq);
    linop_d->ThreadBossLog("HDCG: Pcg k= %d residual = %le %le\n",k,rrn,rtn);

    if ( ipcg ) {
      linop_d->axpy(rp,r,r,0.0);
    }

    // Stopping condition
    if ( rn <= rsq ) { 

      t1=GetTimeBase();

	if ( AnalyseSpectrum ) { 
	  SpectralDecomposition<double>(r,linop_d,"accumulated residual");
	}

	linop_d->ThreadBossLog("HDCG: Pcg converged in %d iterations %f s\n",k,1.0e-6*(t1-t_start)/MHz());
	linop_d->ThreadBossMessage("HDCG: Mprec  %f s\n",1.0e-6*(PcgMprec)/MHz());
	linop_d->ThreadBossMessage("HDCG: MprecO %f s\n",1.0e-6*(PcgMprecO)/MHz());
	linop_d->ThreadBossMessage("HDCG: MssInv %f s\n",1.0e-6*(PcgMssInv)/MHz());
	linop_d->ThreadBossMessage("HDCG: Proj   %f s\n",1.0e-6*(PcgProj)/MHz());
	linop_d->ThreadBossMessage("HDCG: Prom   %f s\n",1.0e-6*(PcgProm)/MHz());
	linop_d->ThreadBossMessage("HDCG: Linalg %f s\n",1.0e-6*(PcgLinalg)/MHz());
	linop_d->ThreadBossMessage("HDCG: Prec   %f s\n",1.0e-6*(PcgPrecChange)/MHz());

	
	linop_d->Mprec(x,mp,tmp,0);
	linop_d->Mprec(mp,mmp,tmp,1); 
	if ( PcgType == PcgAdef2fSingleShift ) { 
	  linop_d->axpy(mmp,x,mmp,PcgSingleShift); 
	}
	linop_d->axpy(tmp,src,mmp,-1.0);

	if ( AnalyseSpectrum ) { 
	  SpectralDecomposition<double>(tmp,linop_d,"true residual");
	  AnalyseSpectrum=0;
	}

	double  mpnorm = sqrt(linop_d->norm(mp));
	double mmpnorm = sqrt(linop_d->norm(mmp));
	double psinorm = sqrt(linop_d->norm(x));
	double srcnorm = sqrt(linop_d->norm(src));
	double tmpnorm = sqrt(linop_d->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	linop_d->ThreadBossLog("HDCG: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	linop_d->ThreadBossMessage("HDCG: target residual was %le \n",linop_d->residual);
	
	
	linop_d->threadedFreeFermion(tmp);
	linop_d->threadedFreeFermion(p);
	linop_d->threadedFreeFermion(z);
	linop_d->threadedFreeFermion(mu);
	linop_d->threadedFreeFermion(mp);
	linop_d->threadedFreeFermion(mmp);
	linop_d->threadedFreeFermion(r);
	linop_d->threadedFreeFermion(rp);

	if ( linop_d->isBoss() && (!me) ) { 
	  linop_d->InverterExit();
	}
	return k;
    }

  }
  linop_d->ThreadBossLog("HDCG: Pcg not converged \n");
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(p);
  linop_d->threadedFreeFermion(mp);
  linop_d->threadedFreeFermion(mmp);
  linop_d->threadedFreeFermion(z);
  linop_d->threadedFreeFermion(mu);
  linop_d->threadedFreeFermion(r);

  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterExit();
  }

  return -1;
}
template<class cFloat> 
int BfmHDCG<cFloat>::gmres(bfm_fermion psi, bfm_fermion src)
{

  const int MaxGMRES=   LdopDeflVecs;
  
  if ( linop_d->SPIcomms() )  linop_d->comm_init();

  bfm_fermion p[MaxGMRES];
  bfm_fermion Ap[MaxGMRES];
  bfm_fermion v[MaxGMRES+1];

  bfm_fermion r;
  bfm_fermion ss;
  bfm_fermion tmp;

  std::vector<std::complex<double> > GmresMatrix(MaxGMRES*MaxGMRES);
  std::vector<std::complex<double> > GmresMatrixInv(MaxGMRES*MaxGMRES);
  std::vector<std::complex<double> > Beta(MaxGMRES);
  std::vector<std::complex<double> > Alpha(MaxGMRES);
  std::vector<double> Apn(MaxGMRES); // normsq of vectors
  double sn,rn,rr;

  for(int cb=0;cb<2;cb++){
    ss[cb] =linop_d->allocFermion();
    r[cb]  =linop_d->allocFermion();
    tmp[cb]=linop_d->allocFermion();
    for(int i=0;i<MaxGMRES+1;i++) { 
      v[i][cb]=linop_d->allocFermion();
    }
    for(int i=0;i<MaxGMRES;i++) { 
      p[i][cb]   = v[i][cb];
      Ap[i][cb]  = v[i+1][cb];
    }
  }

  linop_d->BossLog("restarted GMRES\n");

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      // Source prec. 
      linop_d->G5R(src[0],ss[0]);
      linop_d->G5R(src[1],ss[1]);
      sn =  linop_d->norm(ss);
      // Zero guess
      linop_d->axpby(psi,ss,ss,0.0,0.0);
    }
  }
  linop_d->BossLog("restarted GMRES source norm %le\n",sn);
  int iter = 0;

  do { 
    // Plan:
    //
    // Build up Krylov space K(r,MaxGMRES) = { A^n r ; n=1... MaxGmres} 
    // p[n] = A^n r
    // Ap[n] = A^{n+1} r
    //
    // Update psi to leave residual orthogonal to K(r,MaxGMRES)
    //
    //  r = eta - A psi   => eta - A (psi + Ainv r ) = eta - A [psi + Ainv (eta - A psi) ] = 0
    //
    // Beta[n] = <Ap[n]|r>
    //
    // Mat[n,m] = <Ap[n]|Ap[m]>
    //
    // Alpha[m] = MatInv * Beta
    // 
    // dPsi = Alpha[m] p[m]
    // 
    // A dPsi = Alpha[m] Ap[m]
    // 
    // This implies 
    // <Ap[m] |A dPsi-r> =  <Ap[m]|Ap[n]> Alpha[n] - Beta[m]
    //                   =  Mat[m,n] MatInv[n,k] Beta[k]  - Beta[m]
    //                   =  0
    //
    // Thus removing all components of residual parallel to Krylov space.
    //
    // Later refinement -- try a Blocking scheme and use Ldop to do an approximate Krylov
    // inversion.
    iter ++;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {

	// p0 = r = eta - A psi 

	linop_d->Munprec(psi,tmp,p[0][0],0);
	linop_d->G5R(tmp[0],p[0][0]);
	linop_d->G5R(tmp[1],p[0][1]);
	linop_d->axpy(p[0],p[0],ss,-1.0);

	rn = linop_d->norm(p[0]);
	rr = sqrt(rn/sn);

	linop_d->ThreadBossLog("%d residual %le (res %le / src %le)\n",iter,rr,sqrt(rn),sqrt(sn));
	if ( iter >= 2 ) { 
	  for(int i=0;i<MaxGMRES;i++){
	    std::complex<double> beta=linop_d->inner(Ap[i][0],p[0][0])
                                     +linop_d->inner(Ap[i][1],p[0][1]);
	    double Apnn=linop_d->norm(Ap[i]);
	    linop_d->ThreadBossLog("<residual|Apold > %le \n",abs(beta)/sqrt(Apnn));	    
	  }
	}

	// Solve G5 R5 M psi = G5 R5 eta
	// Minimise residual over Krylov space spanned by A^n r
	
	for(int i=0;i<MaxGMRES;i++) { 
	  
	  // Apply matrix ; note alias p[i+1] and Ap[i].
	  linop_d->Munprec(p[i],tmp,Ap[i][0],0);
	  linop_d->G5R(tmp[0],Ap[i][0]);
	  linop_d->G5R(tmp[1],Ap[i][1]);
	  Apn[i] = sqrt(linop_d->norm(Ap[i]));

	  //	  linop_d->ThreadBossLog("%d Ap %le \n",iter,Apn[i]);
	  
	}
	
	for(int i=0;i<MaxGMRES;i++) { 
	  for(int j=i;j<MaxGMRES;j++) { // 2x too much work here
	    std::complex<double> mm=(linop_d->inner(Ap[i][0],Ap[j][0])
 	 	  	            +linop_d->inner(Ap[i][1],Ap[j][1]))/(Apn[i]*Apn[j]);
	    if ( iter <= 2 ) {
	      //	      linop_d->ThreadBossLog("Mat[%d,%d] = (%le,%le) \n",i,j,real(mm),imag(mm));
	    }
	    if ( i == j ) { 
	      GmresMatrix[j+MaxGMRES*i] =real(mm);
	    } else {
	      GmresMatrix[j+MaxGMRES*i] =mm;
	      GmresMatrix[i+MaxGMRES*j] =conj(mm);
	    }
	    GmresMatrixInv[j+MaxGMRES*i] =GmresMatrix[j+MaxGMRES*i];
	    GmresMatrixInv[i+MaxGMRES*j] =GmresMatrix[i+MaxGMRES*j];
	  }}
	for(int i=0;i<MaxGMRES;i++) { 
	  std::complex<double> beta=(linop_d->inner(Ap[i][0],p[0][0])+linop_d->inner(Ap[i][1],p[0][1]))/Apn[i];
	  Beta[i] = beta;
	  //	  linop_d->ThreadBossLog("Beta[%d] = (%le,%le) \n",i,real(beta),imag(beta));
	}
    }}

#if 0
  // This is a really crappy inverse apparently as test below fails frequently
  LapackHermitianInvert(MaxGMRES,(double *)&GmresMatrixInv[0]);

#if 1
  int N = MaxGMRES;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    std::complex<double> dot=0.0;
    for(int k=0;k<N;k++){
      dot += GmresMatrixInv[N*i+k]*GmresMatrix[k*N+j];
    }
    std::complex<double> expect;
    if(i==j) expect=1.0;
    else     expect=0.0;
    double re = real(expect) - real(dot);
    double im = imag(expect) - imag(dot);
    double ndiff = sqrt(re*re+im*im);
    if ( ndiff > 1.0e-4 ) { 
      printf("Oops GMRES inverse %d %d test failed (%le,%le) (%le,%le) : %le %le\n",i,j,real(dot),imag(dot),real(expect),imag(expect),re,im);
    }
  }}
#endif

  for(int i=0;i<MaxGMRES;i++) { 
    for(int j=0;j<MaxGMRES;j++) { // 2x too much work here
      std::complex<double> mm=GmresMatrix[i+MaxGMRES*j];
      //      linop_d->BossLog("Mat[%d,%d] = (%le,%le) \n",i,j,real(mm),imag(mm));
    }
  }

  for(int i=0;i<MaxGMRES;i++) { 
    Alpha[i]=0.0;
    for(int j=0;j<MaxGMRES;j++) { 
      Alpha[i]=Alpha[i]+ GmresMatrixInv[j+MaxGMRES*i]*Beta[j];
    }
    std::complex<double> mm = Alpha[i];
    linop_d->BossLog("Alpha[%d] = (%le,%le) \n",i,real(mm),imag(mm));
  }
#else 
  // Apply CGNE to solve alpha = MatInv beta
  gmres_solve(MaxGMRES,Alpha,Beta,GmresMatrix);
#endif 


 
  //  dpsi ~ Ainv r 
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop_d->nthread;t++) {
      for(int i=0;i<MaxGMRES;i++) { 
	for(int cb=0;cb<2;cb++){
	  linop_d->caxpy(psi[cb],p[i][cb],psi[cb],real(Alpha[i])/Apn[i],imag(Alpha[i])/Apn[i]);
	}
      }
    }
  }

  } while ( rr>1.0e-8 ) ;


  for(int cb=0;cb<2;cb++){
    linop_d->freeFermion(r[cb]);
    linop_d->freeFermion(ss[cb]);
    linop_d->freeFermion(tmp[cb]);
    for(int i=0;i<MaxGMRES+1;i++) { 
      linop_d->freeFermion(v[i][cb]);
    }
  }
}

template<class cFloat> 
int BfmHDCG<cFloat>::gmres_solve(int N, std::vector<std::complex<double> > &sol,std::vector<std::complex<double> > &src, std::vector<std::complex<double> > &A)
{
  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;

  std::vector<std::complex<double> > ssrc(N);
  std::vector<std::complex<double> > p(N);
  std::vector<std::complex<double> > Ap(N);
  std::vector<std::complex<double> > r(N);
  std::vector<std::complex<double> > tmp(N);

  // For cgne need to solve A^2 sol = Asrc
  // ssrc = A src
  for(int i=0;i<N;i++) {
    ssrc[i]=0;
    for(int j=0;j<N;j++) {
      ssrc[i]+=conj(A[j*N+i])*src[j];
    }
  }

  // guess 0
  // res = ssrc
  // p = res
  for(int i=0;i<N;i++) sol[i]=0;
  for(int i=0;i<N;i++)   r[i]=ssrc[i];
  for(int i=0;i<N;i++)   p[i]=r[i];

  ssq=a=cp=0;
  for(int i=0;i<N;i++){
    ssq+=real(conj(ssrc[i])*ssrc[i]);
    a  +=real(conj(p[i])*p[i]);
    cp +=real(conj(r[i])*r[i]);
  }

  rsq = 1.0e-24;

  for(int k=1;k<10000;k++){

    c=cp;

    double nnn=0;
    // A^2 to keep posdef for herm indef A
    for(int i=0;i<N;i++){
      tmp[i]=0;
      for(int j=0;j<N;j++){
	tmp[i]+=A[i*N+j]*p[j];  // A mult
      }
      nnn+=real(conj(tmp[i])*tmp[i]);
    }
    for(int i=0;i<N;i++){
      Ap[i]=0;
      for(int j=0;j<N;j++){
	Ap[i]+=conj(A[j*N+i])*tmp[j]; // A mult
      }
    }

    d=0;
    for(int i=0;i<N;i++) d+=real(conj(p[i])*Ap[i]);
    a=c/d;
    
    for(int i=0;i<N;i++)   r[i]=r[i]  -a*Ap[i];
    for(int i=0;i<N;i++) sol[i]=sol[i]+a*p[i];

    cp=0; for(int i=0;i<N;i++) cp+=real(conj(r[i])*r[i]);
    b=cp/c;
    for(int i=0;i<N;i++)    p[i]  =  r[i]+b*p[i];

    linop_d->BossLog("gmres_solve:: k=%d cp= %le b = %le d=%le a=%le nnn=%le\n",k,cp,b,d,a,nnn);


      for(int i=0;i<N;i++){
	tmp[i]=0;
	for(int j=0;j<N;j++){
	  tmp[i]+=A[i*N+j]*sol[j];
	}	
      }
      for(int i=0;i<N;i++){
	Ap[i]=0;
	for(int j=0;j<N;j++){
	  Ap[i]+=conj(A[j*N+i])*tmp[j];
	}	
      }
      double nv = 0;
      for(int i=0;i<N;i++){
	Ap[i]=Ap[i]-ssrc[i];
	nv+=real(conj(Ap[i])*Ap[i]);
      }

      linop_d->BossMessage("gmres true residual %le %g\n",nv,sqrt(nv/ssq));
      
    if(cp<rsq) {

      linop_d->BossMessage("gmres via CG converged : %d iterations , residual %g\n",k,
			   sqrt(nv/ssq));
      
      return 0;
    }
  }
  return 0;
}

template<class cFloat> 
int BfmHDCG<cFloat>::fPcg(BfmHDCG<double> * guess_solver,Fermion_t psi, Fermion_t src,Fermion_t resid)
{
  double f;
  double rtzp,rtz,a,d,b;
  double rptzp;

  int me = linop_d->thread_barrier();

  if ( linop_d->SPIcomms() && !me )  linop_d->comm_init();

  PcgReport();

  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterEnter();
  }


  Fermion_t x   = psi;
  std::vector<Fermion_t> p;
  std::vector<Fermion_t> mmp;
  std::vector<double> pAp;
  Fermion_t z   = linop_d->threadedAllocFermion(); 

  Fermion_t tmp = linop_d->threadedAllocFermion(); 
  Fermion_t mp  = linop_d->threadedAllocFermion(); 
  Fermion_t r   = linop_d->threadedAllocFermion(); 
  Fermion_t mu   = linop_d->threadedAllocFermion(); 

  /////////////////////////////
  // Set up history vectors
  /////////////////////////////
  int mmax = 5;
  p.resize(mmax);
  mmp.resize(mmax);
  pAp.resize(mmax);
  
  for (int i=0;i<mmax;i++)  p[i] = linop_d->threadedAllocFermion(); 
  for (int i=0;i<mmax;i++)mmp[i] = linop_d->threadedAllocFermion(); 
  linop_d->thread_barrier();

  //Initial residual computation & set up
  double guess = linop_d->norm(psi);
  double tn;

  //////////////////////////
  // x0 = Vstart -- possibly modify guess
  //////////////////////////
  PcgVstart(x,src,r,mp,mmp[0],tmp);

  // r0 = b -A x0
  linop_d->Mprec(x,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,mmp[0],tmp,DaggerYes);
  if ( PcgType == PcgAdef2fSingleShift ) { 
    linop_d->axpy(mmp[0],x,mmp[0],PcgSingleShift);
  }

  linop_d->axpy (r, mmp[0], src,-1.0);    // Recomputes r=src-Ax0

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  linop_d->ThreadBossMessage("HDCG: subspace residual for x0 = %le\n",nv);

  //////////////////////////////////
  // Compute z = M1 x
  //////////////////////////////////
  PcgM1(r,z,tmp,mp,SmootherMirs);
  rtzp =real(linop_d->inner(r,z));

  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  // Def2: p = 1 - Q Az = Pright z 
  // Other algos M2 is trivial
  ///////////////////////////////////////
  PcgM2(z,p[0]);

  double ssq =  linop_d->norm(src);
  double rsq =  linop_d->residual* linop_d->residual*ssq;

  linop_d->ThreadBossLog("HDCG: k=0 residual %le rsq %le\n",rtzp,rsq);

  uint64_t t_start=GetTimeBase();
  uint64_t t1;

  Fermion_t pp;
  static  uint64_t Tinner;
  Tinner=0;
  LittleDopSolverResidual=LittleDopSolverResidualInner;
  linop_d->thread_barrier();
  for (int k=0;k<=linop_d->max_iter;k++){
    
    int peri_k  = k % mmax;
    int peri_kp = (k+1) % mmax;

    rtz=rtzp;
    d= PcgM3(p[peri_k],mp,mmp[peri_k],tmp);
    a = rtz/d;
    
    // Memorise this
    pAp[peri_k] = d;

    t1 = GetTimeBase();
    linop_d->axpy(x,p[peri_k],x,a);
    double rn = linop_d->axpy_norm(r,mmp[peri_k],r,-a);
    if ( !me )     PcgLinalg+=GetTimeBase()-t1;

    // Compute z = M x
    if ( ParamSmoother == SmootherMirsPoly ) {
      if ( (k < 5)  ||  ((k%10) == 0) ) { 
	PcgM1(r,z,tmp,mp,SmootherMirs);
      } else { 
	PcgM1(r,z,tmp,mp,SmootherMirsPoly);
      }
    } else {  // Freeze out the polynomials
      PcgM1(r,z,tmp,mp,ParamSmoother);
    }

    t1 = GetTimeBase();
    rtzp =real(linop_d->inner(r,z));
    if ( !me )     PcgLinalg+=GetTimeBase()-t1;
    if ( !me )     Tinner+=GetTimeBase()-t1;

    PcgM2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

    t1 = GetTimeBase();
    linop_d->axpy(p[peri_kp],p[peri_k],mu,0);
    // Standard search direction  p -> z + b p    ; b = 
    b = (rtzp)/rtz;

    int northog;

    // k=zero  <=> peri_kp=1;        northog = 1
    // k=1     <=> peri_kp=2;        northog = 2
    // ...               ...                  ...
    // k=mmax-2<=> peri_kp=mmax-1;   northog = mmax-1
    // k=mmax-1<=> peri_kp=0;        northog = 1

    //    northog     = (peri_kp==0)?1:peri_kp; // This is the fCG(mmax) algorithm
    northog     = (k>mmax-1)?(mmax-1):k;        // This is the fCG-Tr(mmax-1) algorithm
    
    //    linop_d->ThreadBossMessage("HDCG::fPcg iteration %d : orthogonalising to last %d vectors\n",k,northog);
    for(int back=0; back < northog; back++){
      int peri_back = (k-back)%mmax;
      uint64_t tx = GetTimeBase();
      double pbApk= real(linop_d->inner(mmp[peri_back],p[peri_kp]));
      if ( !me )     Tinner+=GetTimeBase()-tx;
      double beta = -pbApk/pAp[peri_back];
      linop_d->axpy(p[peri_kp],p[peri_back],p[peri_kp],beta);
    }

    if ( !me )     PcgLinalg+=GetTimeBase()-t1;

    double rrn=sqrt(rn/ssq);
    double rtn=sqrt(rtz/ssq);
    double rtnp=sqrt(rtzp/ssq);
    linop_d->ThreadBossMessage("HDCG: fPcg k= %d residual = %le %le %le\n",k,rrn,rtn,rtnp);

    // Stopping condition
    if ( rn <= rsq ) { 

      t1=GetTimeBase();

	linop_d->ThreadBossLog("HDCG: fPcg converged in %d iterations %f s\n", k,1.0e-6*(t1-t_start)/MHz());
	linop_d->ThreadBossMessage("HDCG: Mprec  %f s\n",1.0e-6*(PcgMprec)/MHz());
	linop_d->ThreadBossMessage("HDCG: MprecO %f s\n",1.0e-6*(PcgMprecO)/MHz());
	linop_d->ThreadBossMessage("HDCG: MssInv %f s\n",1.0e-6*(PcgMssInv)/MHz());
	linop_d->ThreadBossMessage("HDCG: Proj   %f s\n",1.0e-6*(PcgProj)/MHz());
	linop_d->ThreadBossMessage("HDCG: Prom   %f s\n",1.0e-6*(PcgProm)/MHz());
	linop_d->ThreadBossMessage("HDCG: Linalg %f s\n",1.0e-6*(PcgLinalg)/MHz());
	linop_d->ThreadBossMessage("HDCG: Inner  %f s\n",1.0e-6*(Tinner)/MHz());
	linop_d->ThreadBossMessage("HDCG: Prec   %f s\n",1.0e-6*(PcgPrecChange)/MHz());
	
	linop_d->Mprec(x,mp,tmp,0);
	linop_d->Mprec(mp,mmp[0],tmp,1); 
	if ( PcgType == PcgAdef2fSingleShift ) { 
	  linop_d->axpy(mmp[0],x,mmp[0],PcgSingleShift); 
	}
	linop_d->axpy(tmp,src,mmp[0],-1.0);
	
	double psinorm = sqrt(linop_d->norm(x));
	double srcnorm = sqrt(linop_d->norm(src));
	double tmpnorm = sqrt(linop_d->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	linop_d->ThreadBossLog("HDCG: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	linop_d->ThreadBossMessage("HDCG: target residual was %le \n",linop_d->residual);
	
	
	linop_d->threadedFreeFermion(tmp);
	linop_d->threadedFreeFermion(z);
	linop_d->threadedFreeFermion(mu);
	linop_d->threadedFreeFermion(mp);
	for(int i=0;i<mmax;i++){
	  linop_d->threadedFreeFermion(p[i]);
	  linop_d->threadedFreeFermion(mmp[i]);
	}
	linop_d->threadedFreeFermion(r);

	if ( linop_d->isBoss() && (!me) ) { 
	  linop_d->InverterExit();
	}
	return k;
    }

  }
  linop_d->ThreadBossLog("HDCG: not converged \n");
  linop_d->threadedFreeFermion(tmp);
  for(int i=0;i<mmax;i++){
    linop_d->threadedFreeFermion(p[i]);
    linop_d->threadedFreeFermion(mmp[i]);
  }
  linop_d->threadedFreeFermion(mp);
  linop_d->threadedFreeFermion(z);
  linop_d->threadedFreeFermion(mu);
  linop_d->threadedFreeFermion(r);


  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterExit();
  }

  return -1;
}



template<class cFloat> 
void BfmHDCG<cFloat>::Pleft (Fermion_t in,Fermion_t out)
{
  // P_L  = [ 1  -Mbs Mss^-1] 
  //        [ 0   0         ] 
  Fermion_t in_sbar = linop_d->threadedAllocFermion();
  Fermion_t tmp2= linop_d->threadedAllocFermion();
  Fermion_t Mtmp= linop_d->threadedAllocFermion();

  int me = linop_d->thread_barrier();

  uint64_t t1=GetTimeBase();
  ProjectToSubspace(in,PleftProj);     
  PromoteFromSubspace(PleftProj,out);  
  linop_d->axpy(in_sbar,out,in,-1.0);      // in_sbar = in - in_s

  uint64_t t2=GetTimeBase();
  ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} in_s
  uint64_t t3=GetTimeBase();
  PromoteFromSubspace(PleftMss_proj,out);
  uint64_t t4=GetTimeBase();

  linop_d->Mprec(out,tmp2,Mtmp,DaggerNo);  // M Mss^{-1} in_s
  linop_d->Mprec(tmp2,out,Mtmp,DaggerYes);    
  uint64_t t5=GetTimeBase();

  ProjectToSubspace(out,PleftProj);      // Msbar s Mss^{-1}
  PromoteFromSubspace(PleftProj,tmp2);
  linop_d->axpy(out,tmp2,out,-1.0);

  linop_d->axpy(out,out,in_sbar,-1.0);     // in_sbar - Msbars Mss^{-1} in_s
  uint64_t t6=GetTimeBase();

  linop_d->ThreadBossPerformance("Pleft cyc ProjProm/AInv/Prom/MdagM/ProjProm\n\t%10ld\n\t%10ld\n\t%10ld\n\t%10ld\n\t%10ld\n",t2-t1,t3-t2,t4-t3,t5-t4,t6-t5);

  linop_d->threadedFreeFermion(in_sbar);
  linop_d->threadedFreeFermion(Mtmp);
  linop_d->threadedFreeFermion(tmp2);

}
template<class cFloat> 
void BfmHDCG<cFloat>::Pright(Fermion_t in,Fermion_t out)
{
  // P_R  = [ 1              0 ] 
  //        [ -Mss^-1 Msb    0 ] 
  Fermion_t in_sbar = linop_d->threadedAllocFermion();
  Fermion_t tmp = linop_d->threadedAllocFermion();
  Fermion_t Mtmp= linop_d->threadedAllocFermion();

  ProjectToSubspace(in,PleftProj);     
  PromoteFromSubspace(PleftProj,out);  
  linop_d->axpy(in_sbar,out,in,-1.0);       // in_sbar = in - in_s 

  linop_d->Mprec(in_sbar,tmp,Mtmp,DaggerNo);// M in_sbar
  linop_d->Mprec(tmp,out,Mtmp,DaggerYes);  
  ProjectToSubspace(out,PleftProj);           // Mssbar in_sbar  (project)

  ApplyInverse     (PleftProj,PleftMss_proj); // Mss^{-1} Mssbar 
  PromoteFromSubspace(PleftMss_proj,out);     // 

  linop_d->axpy(out,out,in_sbar,-1.0);     // in_sbar - Mss^{-1} Mssbar in_sbar

  linop_d->threadedFreeFermion(in_sbar);
  linop_d->threadedFreeFermion(Mtmp);
  linop_d->threadedFreeFermion(tmp);
}


template void BfmHDCG<float>::PcgMlinop<float>(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_internal<float> *lop,double shift,std::vector<double>&p);
template void BfmHDCG<float>::PcgMlinop<double>(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_internal<double> *lop,double shift,std::vector<double>&p);
template void BfmHDCG<double>::PcgMlinop<float>(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_internal<float> *lop,double shift,std::vector<double>&p);
template void BfmHDCG<double>::PcgMlinop<double>(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_internal<double> *lop,double shift,std::vector<double>&p);

template void BfmHDCG<float>::PolyMdagMprec<float>(bfm_internal<float> *lop,Fermion_t in,Fermion_t out,Chebyshev & approx );
template void BfmHDCG<float>::PolyMdagMprec<double>(bfm_internal<double> *lop,Fermion_t in,Fermion_t out,Chebyshev & approx );
template void BfmHDCG<double>::PolyMdagMprec<float>(bfm_internal<float> *lop,Fermion_t in,Fermion_t out,Chebyshev & approx );
template void BfmHDCG<double>::PolyMdagMprec<double>(bfm_internal<double> *lop,Fermion_t in,Fermion_t out,Chebyshev & approx );

template class BfmHDCG<double>;
template class BfmHDCG<float>;
