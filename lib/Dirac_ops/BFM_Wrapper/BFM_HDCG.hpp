/*!
 * @file BFM_HDCG.hpp
 * @brief Declares classes for P. Boyle HDCG inverter
 * Time-stamp: <2014-08-06 17:34:35 neo>
 */
#ifndef _BFM_MULTIGRID_H_
#define _BFM_MULTIGRID_H_

#include <bfm.h>
#include <vector>
#include <complex>

// I haven't given a name to these types yet. Perhaps do it this way
enum  { LittleDopSolverMCR,LittleDopSolverCG,LittleDopSolverDeflCG, LittleDopSolverADef2 };
enum  PcgAlgorithm { PcgPrec, PcgAD, PcgDef1, PcgDef2, PcgADef1, PcgAdef2, PcgBNN, PcgMssDef, PcgAdef2f,  PcgAdef2fSingleShift, PcgV11f };
enum  HDCGsmoother { SmootherMirs, SmootherNone, SmootherChebyshev, SmootherMirsPoly };
enum  LdopM1ControlType { LdopM1Chebyshev, LdopM1Mirs, LdopM1MirsPoly, LdopM1MirsPolyRecord, LdopM1DiagInv  };

void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out);
void mycgemv_f(int N,float  * __restrict A, float * __restrict in,float c, float * __restrict out);
void mycgemv(int N,float  * __restrict A, double * __restrict in,double c, double * __restrict out);

extern "C" { 
  void qpx_zgemv  (int N,double *  A, double *  in,double c, double *  out);
  void qpx_cgemv  (int N,float  *  A, double *  in,double c, double *  out);
  void qpx_cgemv_f(int N,float  *  A, float *  in ,double c, float *  out);
  void qpx_axpy(int N,double *  R, double *  X, double * Y,double c);
  void   qpx_zaxpy(int N, double *  R,  double *  X,  double * Y,double *c);
  void   qpx_axpby(int N, double *  R,  double *  X,  double * Y,double a,double b);
  void qpx_inner(double *dot,int N,double *A, double *B);
  double qpx_inner_real(int N,double *A, double *B);
  double qpx_norm      (int N,double *  A);

  void   qpx_axpy_f(int N,float *  R, float *  X, float * Y,double c);
  void   qpx_zaxpy_f(int N, float *  R,  float *  X,  float * Y,double *c);
  void   qpx_axpby_f(int N, float *  R,  float *  X,  float* Y,double a,double b);
  void   qpx_inner_f(double *dot,int N,float *A, float *B);
  double qpx_inner_real_f(int N,float *A, float *B);
  double qpx_norm_f      (int N,float *  A);

  void qpx_gather(int Nvec,
		  int *sbuf_gather,
		  int sbuf_gather_size,
		  double *sendbuf,
                  double *vec
		  );
  void qpx_gather_f(int Nvec,
		    int *sbuf_gather,
		    int sbuf_gather_size,
		    float *sendbuf,
		    float *vec
		    );
  void LapackHermitianInvert(int Dim,double *mat);
  void LapackHermitianDiagonalise(int Dim,double *mat,double *evecs,double *evals);
  void LapackHermitianInvert_f(int Dim,float *mat);
  void LapackHermitianDiagonalise_f(int Dim,float *mat,float *evecs,float *evals);
}

class Chebyshev {
public:
  double *Coeffs;
  int order;
  double hi;
  double lo;
  void Init(double lo,double hi,int order, double (* func)(double) );
  void End (void ) {  delete [] Coeffs; }
  double Evaluate(double x);

};

class BfmHDCGParams {
public:
  int    NumberSubspace;
  int    Ls;
  int    Block[5];
  int    SubspaceRationalRefine;
  double SubspaceRationalRefineLo; // Mass, Ls??
  double SubspaceRationalRefineResidual;

  int    SubspaceRationalLs;
  double SubspaceRationalLo;
  double SubspaceRationalMass;
  double SubspaceRationalResidual;

  int    SubspaceSurfaceDepth;

  double LittleDopSolverResidualInner;
  double LittleDopSolverResidualVstart;
  double LittleDopSolverResidualSubspace;
  int    LittleDopSubspaceRational;
  int    LittleDopSolverIterMax;
  int    LittleDopSolver;
  int    LdopDeflVecs;

  // 
  // Controls of smoother
  double PreconditionerKrylovResidual;
  int    PreconditionerKrylovIterMax;
  double PreconditionerKrylovShift;   // Duplicate annoying

  double PcgSingleShift; // The shift 

  // New options
  int SloppyComms;
  LdopM1ControlType LdopM1control;
  double LdopM1Hi;
  double LdopM1Lo;
  double LdopM1resid;
  int    LdopM1iter;


  PcgAlgorithm  PcgType;
  int           Flexible;
  HDCGsmoother  ParamSmoother;

  //  int           OuterSolver;  // FlexCG, CG, GMRES etc...
  //  int           Can drop the 

  BfmHDCGParams() { 
    //    LittleDopSolver = LittleDopSolverDeflCG; 
    LittleDopSolver = LittleDopSolverADef2; 
    ParamSmoother=SmootherMirs;
    LittleDopSubspaceRational=0;
  }
};

class BfmHDCGStatic {
public:
  static int StaticInitialised;
  static std::vector< std::complex<double> * > sendbufs_static; 
  static std::vector< std::complex<double> * > recvbufs_static;
};

template<class cFloat> class BfmHDCG : public BfmHDCGParams, BfmHDCGStatic { 
public:

  int AnalyseSpectrum;

  void SetParams ( BfmHDCGParams & parms ) { 
    *((BfmHDCGParams *) this) = parms;
  }
  void end(void) {};

  // Construct/destroy
  BfmHDCG(int _N5,
	  int _Nvec,
	  int _BlockSize[5],
	  int _QuadrantSize[4],
	  bfm_internal<double> * _linop_d,
	  bfm_internal<float>  * _linop_f) ;

  // Copy
  template<class oFloat> void CloneSubspace(BfmHDCG<oFloat> &ref);
  void FreeDoubleSubspace (void);
  void FreeSingleSubspace (void);

  bfm_internal<double> *linop_d;
  bfm_internal<float>  *linop_f;

  ~BfmHDCG() {
    HaloEnd();
  };

  int MHz(void) { return linop_d->MHz() ; }

  /////////////// Internal information //////////////////
  int BlockSize[5];

  int NblockS;
  int Nvec;
  int Nball;
  int N5;

  int LocalVol;
  int local_nb[5];
  int LocalNblock4;
  int LocalNblock;
  int LocalNsubspace;

  int glatt[5];
  int slatt[5];
  int ncoor[5];

  int global_nb[5];
  int Nblock4;
  int Nblock;
  int Nsubspace;
  int Nquadrant;
  int QuadrantSize[4];
  // Sparse matrix implementation [block][mu][Nvec x Nvec]
  // Within the stencil A becomes dense. We use a subblock to identify "quadrants" within 
  // the "block". quadrant_size = block_size is the original behaviour
  int s_min;
  int s_max;
  int stencil_lo[5];
  int stencil_hi[5];
  int stencil_size[5];
  std::vector<int> StencilNonZero;
  std::vector<int> StencilDepth;
  std::vector<int> StencilReord;
  std::vector<int> StencilDepthCount;
  std::vector<std::vector<int> > block_id;
  std::vector<std::vector<int> > local_block_id;
  std::vector<std::vector<int> > quadrant_id;

  // Sizing
  int SubspaceDimension(void){return Nsubspace;};
  int SubspaceLocalDimension(void){return LocalNsubspace;};
  int NeighbourBlockId(int myblock,int _mu, int &me,int &from,int &xblock,int &to,int rev);
  int NeighbourGlobalBlockId(int localBlock,int mu);
  int GetLocalBlockIdx(int delta[5]);
  void GetDelta(int delta[5],int mu,int rev);
  int ReverseDirection(int mu);
  void GlobalToLocalBlock(int gb[5],int &proc,int lb[5]);

  // Eventually compress the quadrant subdivision to decrease overhead
  int SubspaceIndex(int _block, int _quadrant,int _vec) {    return _block*Nvec*Nquadrant+ _vec*Nquadrant +_quadrant ; }
  int SubspaceIndex(int _block, int _vec) {    return _block*Nvec+ _vec ; }
  void OrthogonaliseSubspace(void);
  template<class lFloat> void ComputeLittleMatrixColored (bfm_internal<lFloat> *op,Fermion_t *subspace);
  void SinglePrecSubspace(void);

  double LittleDopSolverResidual;

  template<class Float> void RelaxSubspace(bfm_internal<Float> *rop);
  template<class Float> void RationalSubspace(bfm_internal<Float> *rop,Fermion_t src,Fermion_t sol);

  ///////////////////////////////////////////////////////
  // Cook up elegant multiprecision implementation.
  ///////////////////////////////////////////////////////
  Fermion_t *subspace_r;  // Non-orthogonalised subspace for refinement
  Fermion_t *subspace_d;  // double subspace, block orthog
  Fermion_t *subspace_f;  // float  subspace, block orthog

  ///////////////////////////////////////////////////////
  // Record and replay the CG polynomial
  ///////////////////////////////////////////////////////
  int CGpolynomialCounter;
  std::vector<double> CGpolynomial;
  std::vector<double> SparePolynomial;

  std::vector<  std::complex<cFloat> > Asparse;
  std::vector<  std::complex<float> >  AsparseSingle;
  std::vector<  std::complex<cFloat> > AsparseDiagInv;

  std::vector<std::complex<double> >  inner_reduce ;

  std::vector<std::complex<cFloat> >  KrylovNest_p ;
  std::vector<std::complex<cFloat> >  KrylovNest_r ;
  std::vector<std::complex<cFloat> >  KrylovNest_Ap ;

  std::vector<std::complex<cFloat> >  Krylov_p ;
  std::vector<std::complex<cFloat> >  Krylov_r ;
  std::vector<std::complex<cFloat> >  Krylov_mu ;
  std::vector<std::complex<cFloat> >  Krylov_Ap;
  std::vector<std::complex<cFloat> >  Krylov_Ar;
  std::vector<std::complex<cFloat> >  Krylov_Amu ;
  std::vector<std::complex<cFloat> >  Krylov_tmp ;

  std::vector<std::complex<cFloat> >  Krylov_Atmp;

  std::vector<std::complex<cFloat> > ProjectInnerProduct;
  std::vector<std::complex<cFloat> > PromoteCoeff;
  std::vector<std::complex<cFloat> > PleftProj;
  std::vector<std::complex<cFloat> > PleftMss_proj;

  template<class datum> std::vector<datum> *ThreadedAllocVector(int nvec);
  template<class datum> void ThreadedFreeVector(std::vector<datum> *);

  std::vector<std::complex<cFloat> > DeflKrylovProj;
  std::vector<std::complex<cFloat> > DeflKrylovMss;
  std::vector<std::complex<double> > DeflKrylovGsum;

  std::vector< std::vector<std::complex<cFloat> > > ComputeProj;
  std::vector<std::complex<double> > phases;

  //////////////////////////////////////////////
  ///////////  Ldop comms ///////////
  //////////////////////////////////////////////
  void HaloExchange(std::vector<std::complex<cFloat> >&my_data, int depth);
  void HaloInit(void);
  void HaloEnd(void);
  std::complex<cFloat> *HaloGetStencilData(int _idx,int _ballidx,std::vector<std::complex<cFloat> > &my_data);

  std::vector< int > tproc; // Processor we must send to and from
  std::vector< int > fproc;

  ///////////////////////////////////////////////////////////////
  // Halo depth 4
  ///////////////////////////////////////////////////////////////
  // How much data we must send and receive
  // buffers we must send and receive
  std::vector< int > sendbuf_bytes;
  std::vector< int > recvbuf_bytes;

  // How big buffers to each neighbour for each depth halo
  std::vector< std::vector<int > > sendbuf_depth_bytes;
  std::vector<int > sbuf_gather_depth_count;

  std::vector< std::complex<cFloat> * > sendbufs; 
  std::vector< std::complex<cFloat> * > recvbufs;






  std::vector<QMP_msgmem_t> xmit_msgmem;  // QMP data structures for the sending
  std::vector<QMP_msgmem_t> recv_msgmem;
  std::vector<QMP_msghandle_t> xmit_msghandle;
  std::vector<QMP_msghandle_t> recv_msghandle;









  std::vector< int > sbuf_id; // Tables mapping local block + mu combination to buffer and offset
  std::vector< int > sbuf_idx;
  std::vector< std::vector<int> > sbuf_gather;

  std::vector<int>  sbuf_gather_from_idx;
  std::vector<int>  sbuf_gather_buff_idx;
  std::vector<int>  sbuf_gather_buff_id;

  std::vector< int > rbuf_id;
  std::vector< int > rbuf_idx;
  
  ///////////////////////////////////////////////////////////////
  // Halo depth 1
  ///////////////////////////////////////////////////////////////


  ////// To from coarse grid //////
  void ProjectToSubspace(Fermion_t, std::vector<std::complex<cFloat> > &proj);
  void ProjectToSubspace_f(Fermion_t, std::vector<std::complex<cFloat> > &proj);
  void PromoteFromSubspace(std::vector<std::complex<cFloat> > &v, Fermion_t prom);
  void PromoteFromSubspace_f(std::vector<std::complex<cFloat> > &v, Fermion_t prom);

  /////// Coarse grid operations //////////////
  void scale(std::vector<std::complex<cFloat> > &r,
	     std::vector<std::complex<cFloat> > &x,
	     double a);
  void copy(std::vector<std::complex<cFloat> > &r,std::vector<std::complex<cFloat> > &x) { scale(r,x,1.0); }

  void zaxpy(std::vector<std::complex<cFloat> > &r,
	     std::vector<std::complex<cFloat> > &x,
	     std::vector<std::complex<cFloat> > &y,std::complex<double> a);

  void axpy(std::vector<std::complex<cFloat> > &r,
	    std::vector<std::complex<cFloat> > &x,
	    std::vector<std::complex<cFloat> > &y,double a);

  void axpby(std::vector<std::complex<cFloat> > &r,
	     std::vector<std::complex<cFloat> > &x,
	     std::vector<std::complex<cFloat> > &y,double a,double b);

  std::complex<double> innerProduct(std::vector<std::complex<cFloat> > &v1,
				    std::vector<std::complex<cFloat> > &v2) ;

  double innerProductReal(std::vector<std::complex<cFloat> > &v1,
			  std::vector<std::complex<cFloat> > &v2);

  void zeroOut(std::vector<std::complex<cFloat> > &v) ;

  double norm_vec(std::vector<std::complex<cFloat> > &v) ;


  // Preconditioner 1/[MdagM + shift]
  template<class Float>   void PcgMlinop(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_internal<Float> *lop,double shift,std::vector<double>&poly);
  //  template<class Float>   void PcgMDDlinop(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_qdp<Float> *lop_bulk,bfm_qdp<Float> *lop_wall);

  // Polynomials of the ldop
  void PolyMdagMprec(std::vector<std::complex<cFloat> > &in,std::vector<std::complex<cFloat> > &out,Chebyshev &Approx);
  void PolyMdagMprec(std::vector<std::complex<cFloat> > &in,std::vector<std::complex<cFloat> > &out,vector<double> &coeffs);

  // Polynomials of the fine dop
  template<class Float>   void PolyMdagMprec(bfm_internal<Float> *lop,Fermion_t in,Fermion_t out,Chebyshev &Approx);
  template<class Float>   void PolyMdagMprec(bfm_internal<Float> *lop,Fermion_t in,Fermion_t out,vector<double> & coeffs);
  template<class Float>   void SpectralDecomposition(Fermion_t psi,bfm_internal<Float> *lop,const char *tag);

  void   PcgM(Fermion_t in,Fermion_t out,Fermion_t tmp); 
  void   PcgM_f(Fermion_t in,Fermion_t out,Fermion_t tmp);  // single prec variant
  void   PcgM1(Fermion_t in, Fermion_t out,Fermion_t tmp,Fermion_t mp,int Smoother=SmootherMirs);
  void   PcgM2(Fermion_t in, Fermion_t out);
  double PcgM3(Fermion_t p, Fermion_t mp,Fermion_t mmp, Fermion_t tmp);
  void   PcgVstart(Fermion_t in, Fermion_t src, Fermion_t r, 
	  	   Fermion_t mp, Fermion_t mmp, Fermion_t tmp);
  void   PcgVout  (Fermion_t in, Fermion_t out,Fermion_t src,Fermion_t tmp);
  void Pleft (Fermion_t in,Fermion_t out);
  void Pright(Fermion_t in,Fermion_t out);


  // Timing in cycles
  uint64_t PcgMssInv;
  uint64_t PcgMprec;
  uint64_t PcgMprecO;
  uint64_t PcgProj;
  uint64_t PcgProm;
  uint64_t PcgLinalg;
  uint64_t PcgPrecChange;

  void   PcgReport (void);
  
  int    Pcg (BfmHDCG<double > * guess_solver,Fermion_t psi, Fermion_t src, Fermion_t resid);
  int    fPcg(BfmHDCG<double > * guess_solver,Fermion_t psi, Fermion_t src, Fermion_t resid);

  // Try a gmres preconditioned
  int    gmres(bfm_fermion psi,bfm_fermion src);
  int    gmres_solve(int N, std::vector<std::complex<double> > &sol,std::vector<std::complex<double> > &src, std::vector<std::complex<double> > &A);


  //  template<class Float> void SolverMatrix (bfm_qdp<Float> *lop,Fermion_t in,Fermion_t out);

  void ApplyThread(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> > & , int single=0,int depth=4);

  void ApplyDiagInv(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> > &);
  void ApplyThreadOpt(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> > & , int single=0, int depth=4);
  void ApplyInverseMultiShift(std::vector<std::complex<cFloat> >&src,
			      std::vector<std::vector<std::complex<cFloat> > >&psi,
			      std::vector<std::vector<std::complex<cFloat> > >&ps,
			      int nshift,
			      double mass[],
			      double alpha[],
			      double residual[] );

  void ApplyInverse(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);
  void ApplyInverseCG(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&, double shift=0,double resid=0,int maxit=0);
  void ApplyInverseCGopt(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&, 
			 std::vector<double> & poly,
			 double shift=0,double resid=0,int maxit=0,int halo=1);
  void ApplyInverseMCR(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);
  void ApplyInverseDeflCG(std::vector<std::complex<cFloat> > &src,std::vector<std::complex<cFloat> >&, double shift=0);
  void ApplyInverseADef2(std::vector<std::complex<cFloat> > &src,std::vector<std::complex<cFloat> >&, double shift=0);
  void LdopM1(std::vector<std::complex<cFloat> > &in,
	      std::vector<std::complex<cFloat> > &out, Chebyshev & Approx,LdopM1ControlType control);

  ///////////////////////////////////////////////
  // Deflating the little dop too!
  ///////////////////////////////////////////////
  int LdopDeflationBasisSize ;
  std::vector<std::vector<std::complex<cFloat> > > LdopDeflationBasis ;
  std::vector<std::vector<std::complex<cFloat> > > LdopDeflationAv;
  std::vector<double> LdopDeflationEvals;
  std::vector<double> LdopDeflationInvEvals;
  std::vector<std::complex<cFloat> >   LdopDeflationMatrix;
  std::vector<std::complex<cFloat> >   LdopDeflationInverseMatrix;
  int  LdopDeflationIsDiagonal;
  void LdopDeflationBasisInit(int Nbasis);
  void LdopDeflationBasisTrivial(void);
  void LdopDeflationBasisDiagonalise(int Nbasis);
  void LdopDeflationProject      (std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);
  void LdopDeflationPromote      (std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);
  void LdopDeflationMatPromote   (std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);
  void LdopDeflationMatrixInverseMult(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);
  void LdopDeflationMatrixMult       (std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> >&);

};


/////////////////////////////////
// Routines that call out to SPI 
/////////////////////////////////
template<class cFloat> void HaloInitBuffers(BfmHDCG<cFloat> *BMG);
template<class cFloat> void HaloExchangeCommStart(BfmHDCG<cFloat> *BMG,int depth);
template<class cFloat> void HaloExchangeCommComplete(BfmHDCG<cFloat> *BMG,int depth);


#endif
