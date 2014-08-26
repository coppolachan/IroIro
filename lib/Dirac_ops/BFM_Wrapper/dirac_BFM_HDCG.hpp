/*!
 * @file dirac_BFM_HDCG_wrapper.hpp
 * @brief Declares the wrapper class for P. Boyle HDCG inverter
 * Time-stamp: <2014-08-07 15:31:29 neo>
 */
#ifndef DIRAC_BFM_HDCG_WRAPPER_
#define DIRAC_BFM_HDCG_WRAPPER_

#include "dirac_BFM_wrapper.hpp"
#include "BFM_HDCG.hpp"

/* 
 *! @class Dirac_BFM_HDCG_Wrapper
 * @brief Peter Boyle's Hierarchically deflated CG (HDCG)
 * Converted from BFM (Chroma++ style)
 * Paper: arXiv 1402.2585
 * Subclass from Dirac_BFM_Wrapper
 */
class Dirac_BFM_HDCG_Wrapper: public Dirac_BFM_Wrapper {
private:
  typedef float rFloat;

  BFM_HDCG_Extend<double> *ldop_d;
  BfmHDCGParams HDCGParms;
  bfm_internal<rFloat> linop_r;

  //Hide default constructor
  Dirac_BFM_HDCG_Wrapper();
  Dirac_BFM_HDCG_Wrapper(Dirac_BFM_HDCG_Wrapper&);

  
public:
  Dirac_BFM_HDCG_Wrapper(XML::node node,
			 const Field* F, 
			 DiracWilsonLike_EvenOdd* DWEO,
			 const Type_5d_DWF = Regular5D)
    : Dirac_BFM_Wrapper(node, F, DWEO) {}
    

  void HDCG_init(BfmHDCGParams & Parms_, bfmActionParams &BAP);
  
  void HDCG_subspace_rotate(double *phases_by_pi,double *phases_dir, int sgn);
  void HDCG_set_mass(double mass);
  
  // Low modes subspace 
  void HDCG_subspace_init();
  void HDCG_subspace_refine();
  void HDCG_subspace_free();
  void HDCG_subspace_compute(int sloppy);
  
  void solve_HDCG(Fermion_t solution[2], Fermion_t source[2], double residual, int maxit);
  void solve_HDCG(FermionField &sol, FermionField &src,double residual, int maxit);
    
};


#endif
