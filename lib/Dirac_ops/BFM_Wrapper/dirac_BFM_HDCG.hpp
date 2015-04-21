/*!
 * @file dirac_BFM_HDCG_wrapper.hpp
 * @brief Declares the wrapper class for P. Boyle HDCG inverter
 * Time-stamp: <2015-04-17 23:13:03 cossu>
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

  BfmHDCGParams HDCGParams;
  BFM_HDCG_Extend<double> *ldop_d;  // internal double precision operator

  bfm_internal<rFloat> linop_r;     // internal operator
  BFM_Storage<rFloat>  BFM_interface_r;

  //Hide default constructor
  Dirac_BFM_HDCG_Wrapper();
  Dirac_BFM_HDCG_Wrapper(Dirac_BFM_HDCG_Wrapper&);
  
public:
  Dirac_BFM_HDCG_Wrapper(XML::node node,
			 const Field* F, 
			 DiracWilsonLike_EvenOdd* DWEO,
			 const Type_5d_DWF = Regular5D)
    : Dirac_BFM_Wrapper(node, F, DWEO), BFM_interface_r(linop_r) {}
    
  void open_comms();
  void close_comms();


  void HDCG_init(XML::node);// keep separated from the constructor
  
  void HDCG_subspace_rotate(double *phases_by_pi,double *phases_dir, int sgn);
  void HDCG_set_mass(double mass);
  
  // Low modes subspace 
  void HDCG_subspace_init();
  void HDCG_subspace_refine();
  void HDCG_subspace_free();
  void HDCG_subspace_compute(int sloppy);
  
  // Actual solve (one the subspace is ready)
  void solve_HDCG(Fermion_t solution[2], Fermion_t source[2], double residual, int maxit, int cb);
  void solve_HDCG(FermionField &sol, FermionField &src);

  bool do_refine(){return bool(HDCGParams.SubspaceRationalRefine);};
    
};


#endif
