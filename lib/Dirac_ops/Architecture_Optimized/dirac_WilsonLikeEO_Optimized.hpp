/*
 * @file dirac_WilsonLikeEO_Optimized.hpp
 *
 * @brief Definition of Solver_CG_DWF_Optimized class
 *
 * Time-stamp: <2013-07-04 17:20:20 cossu> 
 */
#include "Dirac_ops/dirac_WilsonLike.hpp"

class DiracWilsonLike_Optimized: public DiracWilsonLike {
public:
  virtual void solve(Field&, 
		     const Field&, 
		     SolverOutput&, 
		     int, 
		     double);
  
  virtual void solve_multishift(prop_t& ,
				const Field& , 
				SolverOutput&, 
				const std::vector<double>&, 
				int,
				double);

};

class DiracDWF_Optimized_BGQ: public DiracWilsonLike_Optimized {
  Dirac_optimalDomainWall_EvenOdd *DWF;
public:
 
  void solve(Field&, 
	     const Field&, 
	     SolverOutput&, 
	     int, 
	     double);
  
  void solve_multishift(prop_t& ,
			const Field& , 
			SolverOutput&, 
			const std::vector<double>&, 
			int,
			double);
  
};

class Dirac_Optimized_BFM: public DiracWilsonLike_Optimized {
public:
  void solve(Field&, 
	     const Field&, 
	     SolverOutput&, 
	     int, 
	     double);
  
  void solve_multishift(prop_t& ,
			const Field& , 
			SolverOutput&, 
			const std::vector<double>&, 
			int,
			double);
  
};


