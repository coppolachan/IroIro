/*!
 * @file rationalSolver_DWF_BGQ.cpp
 * @brief Wrapper around MultiShiftSolver_CG class
 * Calculates rational functions approximations
 */

#include "Solver/Architecture_Optimized/rationalSolver_DWF_BGQ.hpp"
#include "Tools/Architecture_Optimized/utils_BGQ.hpp"

#ifdef HAVE_LIBBFM
RationalSolver_DWF_Optimized::RationalSolver_DWF_Optimized(Dirac_BFM_Wrapper* BFMopr,
							   const XML::node node)
  :BFM_opr_(BFMopr),
   Params(MultiShiftSolver_CG_Params(node)),
   Residuals(0),
   Poles(0),
   ConstTerm(0),
   InvResiduals(0),
   InvPoles(0),
   InvConstTerm(0),
   is_BFM(true),
   Nsites(BFM_opr_->getN5()){
  BFM_opr_->set_SolverParams(node); 
};
#endif

RationalSolver_DWF_Optimized::RationalSolver_DWF_Optimized(const Dirac_DomainWall_EvenOdd* DWFopr,
							   const XML::node node,
							   RationalApprox& RA)
  :opr_(DWFopr),
   Params(MultiShiftSolver_CG_Params(node)),
   Residuals(RA.Residuals()),
   Poles(RA.Poles()),
   ConstTerm(RA.Const()),
   InvResiduals(RA.InvResiduals()),
   InvPoles(RA.InvPoles()),
   InvConstTerm(RA.InvConst()),
   Nsites(opr_->getN5()){};

// To be used in Action constructions
RationalSolver_DWF_Optimized::RationalSolver_DWF_Optimized(const Dirac_DomainWall_EvenOdd* DWFopr,
							   const XML::node node)
  :opr_(DWFopr),
   Params(MultiShiftSolver_CG_Params(node)),
   Residuals(0),
   Poles(0),
   ConstTerm(0),
   InvResiduals(0),
   InvPoles(0),
   InvConstTerm(0),
   Nsites(opr_->getN5()){};

void RationalSolver_DWF_Optimized::internal_solve(vector_Field& v_sol, 
						  const Field& source,
						  const vector_double& Shifts,
						  SolverOutput& out) const {
#ifdef HAVE_LIBBFM
  if (is_BFM){
    std::vector < FermionField > BFM_ms_solution(Shifts.size());
    FermionField source_FF(source);
    vector_double mresiduals(size_t(Shifts.size()), sqrt(Params.GoalPrecision));
    BFM_opr_->solve_CGNE_multishift(BFM_ms_solution,source_FF, Shifts, mresiduals);
    
    for (int i=0; i< v_sol.size(); ++i) {
      v_sol[i] = BFM_ms_solution[i].data;
    }
  } else 
#endif
    {
      opr_->solve_ms_eo(v_sol, source, out, Shifts, Params.MaxIter, Params.GoalPrecision);
    }

}

SolverOutput RationalSolver_DWF_Optimized::solve(Field& sol, 
						 const Field& source) const {
  SolverOutput out;
  if (Residuals.size() == 0)  {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  vector_Field shifted_sol;
  sol.resize(source.size());
  shifted_sol.resize(Residuals.size());
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
  
  internal_solve(shifted_sol, source, Poles, out);
  
  // Reconstruct solution (M^dag M)^(a/b)
  
  Spinor *sol_ptr = (Spinor*)sol.getaddr(0);
  Spinor *source_ptr = (Spinor*)(const_cast< Field& >(source).getaddr(0));
 
  int Nvol = CommonPrms::Nvol()/2; //even-odd

#pragma omp parallel 
  {
    int nid = omp_get_num_threads();
    int is = omp_get_thread_num()*Nvol/nid;
    int ns = Nvol/nid; 
    for (int s = 0; s < Nsites; s++)
      BGWilsonLA_MultScalar(sol_ptr+s*Nvol+is, 
			    source_ptr+s*Nvol+is, 
			    ConstTerm, ns);
    
    
    for (int i = 0; i < Poles.size(); ++i){
      Spinor *sh_sol_ptr = (Spinor*)(shifted_sol[i]).getaddr(0);
      for (int s = 0; s < Nsites; s++){
	BGWilsonLA_MultAddScalar(sol_ptr+s*Nvol+is, 
				 sh_sol_ptr+s*Nvol+is, 
				 Residuals[i], ns);
      }
    }
  }
  
  return out;
}

SolverOutput RationalSolver_DWF_Optimized::solve_inv(Field& sol, 
						     const Field& source) const {
  SolverOutput out;
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
 
  vector_Field shifted_sol;
  shifted_sol.resize(InvResiduals.size());
  sol.resize(source.size());

  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
  

  //opr_->solve_ms_eo(shifted_sol, source, out, InvPoles, Params.MaxIter, Params.GoalPrecision);
  internal_solve(shifted_sol, source, InvPoles, out);

  // Reconstruct solution (M^dag M)^(-a/b)
  /*
    Field temp;
    temp.resize(source.size());
    sol = source;
    sol *= InvConstTerm;
    
    for (int i = 0; i < InvPoles.size(); ++i){
    temp = shifted_sol[i];
    temp *= InvResiduals[i];
    sol += temp;
  }
  */
  Spinor *sol_ptr = (Spinor*)sol.getaddr(0);
  Spinor *source_ptr = (Spinor*)(const_cast< Field& >(source).getaddr(0));
 
  int Nvol = CommonPrms::Nvol()/2; //even-odd

#pragma omp parallel 
  {
    int nid = omp_get_num_threads();
    int is = omp_get_thread_num()*Nvol/nid;
    int ns = Nvol/nid; 
    for (int s = 0; s < Nsites; s++)
      BGWilsonLA_MultScalar(sol_ptr+s*Nvol+is, 
			    source_ptr+s*Nvol+is, 
			    InvConstTerm, ns);
    
    
    for (int i = 0; i < InvPoles.size(); ++i){
      Spinor *sh_sol_ptr = (Spinor*)(shifted_sol[i]).getaddr(0);
      for (int s = 0; s < Nsites; s++){
	BGWilsonLA_MultAddScalar(sol_ptr+s*Nvol+is, 
				 sh_sol_ptr+s*Nvol+is, 
				 InvResiduals[i], ns);
      }
    }
  }
  return out;
}


// Needed in force calculation
// It solves the inverse equation
SolverOutput RationalSolver_DWF_Optimized::solve_noReconstruct(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  SolverOutput out;
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  shifted_sol.resize(InvResiduals.size());
 
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }

  //opr_->solve_ms_eo(shifted_sol, source, out, InvPoles, Params.MaxIter, Params.GoalPrecision);
  internal_solve(shifted_sol, source, InvPoles, out);

  return out;
}


// Needed in force calculation
// It solves the direct equation (nevertheless the name says "inv" because in the action 
// such a term is associated to the inverse operator
SolverOutput RationalSolver_DWF_Optimized::solve_noReconstruct_inv(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  SolverOutput out;
  if (Residuals.size() == 0) {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  shifted_sol.resize(Residuals.size());
 
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
 
  //opr_->solve_ms_eo(shifted_sol, source, out, Poles, Params.MaxIter, Params.GoalPrecision);
  internal_solve(shifted_sol, source, Poles, out);

  return out;
}
