/*!--------------------------------------------------------------------------
 * @file domainWallSolver_BGQ.cpp
 *
 * @brief Definition of BGQ optimized solvers for Dirac_DomainWall (5d operator)
 *Time-stamp: <2013-12-10 18:30:57 noaki>
 *-------------------------------------------------------------------------*/
#include "Tools/utils_BGQ.hpp"
#include "Solver/solver.hpp"
#include "Fields/field_expressions.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"

void Dirac_DomainWall::solve_eo_5d(Field& w5,const Field& b, 
				   SolverOutput& Out, 
				   int maxIter, 
				   double targetPrec) const{
#if VERBOSITY>=SOLV_ITER_VERB_LEVEL
  CCIO::header("CG_BGQ solver start");
#endif

#ifdef ENABLE_THREADING
  BGQThread_Init(); //initializing BGQ fast threading routines

  Out.Msg = "CG solver BGQ";
  Out.Iterations = -1;
  double kernel_timing = 0.0; 

  register int Nvol = CommonPrms::instance()->Nvol()/2;
  register int N5 = prms_.N5_;
  double* u_ptr = const_cast<Field * >(Dw_->getGaugeField_ptr())->getaddr(0);

  Field temp(b.size()), temp2(b.size()), s(b.size()); //eventually eliminated
  Field p(b.size());
  int v_size = Nvol*N5; // << assumes 24 elements

  double pap, rrp, cr, rr, snorm;

  TIMING_START;
  Field x = b;//initial condition
  Field r = b;//initial residual

  Spinor* x_ptr = (Spinor*)x.getaddr(0);
  Spinor* r_ptr = (Spinor*)r.getaddr(0);
  Spinor* s_ptr = (Spinor*)s.getaddr(0);
  Spinor* w5_ptr = (Spinor*)w5.getaddr(0);
  Spinor* b_ptr = (Spinor*)const_cast<Field&>(b).getaddr(0);
  Spinor* temp_ptr = (Spinor*)temp.getaddr(0);
  Spinor* temp2_ptr = (Spinor*)temp2.getaddr(0);
  Spinor* p_ptr = (Spinor*)p.getaddr(0);

  BGWilson_DW_Init(N5,prms_.mq_,prms_.M0_,
		   (double*)&prms_.dp_[0],(double*)&prms_.dm_[0],
		   (double*)&prms_.bs_[0],(double*)&prms_.cs_[0],
		   (double*)&prms_.es_[0],(double*)&prms_.fs_[0]);

  double kappa = 0.5/(4.0+prms_.M0_);
 #pragma omp parallel 
   {
     double pap, rrp, cr, rr, snorm;
     double tSum,t, temp;
     double timer;

     int nid = omp_get_num_threads();
     int tid = omp_get_thread_num();
     int ns = Nvol/nid;
     int is = tid*ns;

     BGWilson_DW_Mult_hop(temp_ptr,(void*)(u_ptr),x_ptr, kappa,BGWILSON_DIRAC);
     BGWilson_DW_Mult_hop_dag(temp2_ptr,(void*)(u_ptr),temp_ptr, kappa,BGWILSON_DIRAC);

     for(int s5=0; s5<N5; ++s5)
       BGWilsonLA_Sub(r_ptr+s5*Nvol+is, temp2_ptr+s5*Nvol+is, ns);
     
     for(int s5=0; s5<N5; ++s5)
       BGWilsonLA_Equate(p_ptr+s5*Nvol+is, r_ptr+s5*Nvol+is, ns);

     tSum = 0.0;
     for(int s5=0; s5<N5; ++s5){
       BGWilsonLA_Norm(&t,r_ptr+s5*Nvol+is,ns);
       tSum += t;
     }

     tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
     if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
     rr = BGQThread_ScatterDouble(tSum,0,tid,nid);

     tSum = 0.0;
     for(int s5=0; s5<N5; ++s5){
       BGWilsonLA_Norm(&t,b_ptr+s5*Nvol+is,ns);
       tSum += t;
     }

     tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
     if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
     snorm = BGQThread_ScatterDouble(tSum,0,tid,nid);
     //snorm = sqrt(snorm);
     snorm = 1.0/snorm;

 #if VERBOSITY>=SOLV_ITER_VERB_LEVEL
 #pragma omp single
     {
       CCIO::cout<<" Snorm = "<< snorm << std::endl;
       CCIO::cout<<" Init  = "<< rr*snorm<< std::endl;
     }
 #endif

     for(int it=0; it<maxIter; ++it){
       double tPAP;

       if(tid == 0) FINE_TIMING_START(timer);

       BGWilson_DW_Mult_hop(temp_ptr,(void*)(u_ptr),    p_ptr, kappa,BGWILSON_DIRAC);
       BGWilson_DW_Mult_hop_dag(s_ptr,(void*)(u_ptr),temp_ptr, kappa,BGWILSON_DIRAC);

       ///////////////////////////////////////////////////
       tSum = 0.0;
       for(int s5=0; s5<N5; ++s5){
	 BGWilsonLA_DotProd(&t,p_ptr+s5*Nvol+is,s_ptr+s5*Nvol+is,ns);
	 tSum += t;
       }
       tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
       if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
       pap = BGQThread_ScatterDouble(tSum,0,tid,nid);
       rrp = rr;
       cr = rrp/pap;// (r,r)/(p,Ap)

       //  x += cr*p; // x = x + cr * p
       for(int s5=0; s5<N5; ++s5)
	 BGWilsonLA_MultAddScalar(x_ptr+s5*Nvol+is,p_ptr+s5*Nvol+is,cr,ns);
       //  r -= cr*s; // r_k = r_k - cr * Ap
       for(int s5=0; s5<N5; ++s5)
	 BGWilsonLA_MultAddScalar(r_ptr+s5*Nvol+is,s_ptr+s5*Nvol+is,-cr,ns);
       //  rr = r*r; // rr = (r_k,r_k)
       tSum = 0.0;
       for(int s5=0; s5<N5; ++s5){
	 BGWilsonLA_Norm(&t,r_ptr+s5*Nvol+is,ns);
	 tSum += t;
       }

       tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
       if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
       rr = BGQThread_ScatterDouble(tSum,0,tid,nid);

       for(int s5=0; s5<N5; ++s5)
	 BGWilsonLA_MultScalar_Add(p_ptr+s5*Nvol+is,r_ptr+s5*Nvol+is,rr/rrp,ns);

 #if VERBOSITY>=SOLV_ITER_VERB_LEVEL
       if (tid==0) {
	 CCIO::cout<< std::setw(5)<< "["<<it<<"] "
		   << std::setw(20) << rr*snorm<<"\n";
       } 
 #endif
       BGQThread_Barrier(0,nid);

       if(rr*snorm < targetPrec){
	 if (tid==0) Out.Iterations = it;
	 break;
       }
       if (tid==0){
       FINE_TIMING_END(timer);
	 _Message(DEBUG_VERB_LEVEL,"  CG iteration time: "<<timer<<" seconds\n");
       }
     }//end of iterations loop      

 #pragma omp single
     {
       if(Out.Iterations == -1) {
	 CCIO::cout<<" Not converged. Current residual: "<< rr*snorm << "\n";
	 abort();
       }
     }
     BGQThread_Barrier(0,nid);//after loop

     BGWilson_DW_Mult_hop(temp_ptr,(void*)(u_ptr),x_ptr,kappa,BGWILSON_DIRAC);
     BGWilson_DW_Mult_hop_dag(p_ptr,(void*)(u_ptr),temp_ptr,kappa,BGWILSON_DIRAC);

     //  p -= b;
     for(int s5=0; s5<N5; ++s5)
       BGWilsonLA_Sub(p_ptr+s5*Nvol+is, b_ptr+s5*Nvol+is, ns);

     tSum = 0.0;
     for(int s5=0; s5<N5; ++s5){
       BGWilsonLA_Norm(&t,p_ptr+s5*Nvol+is,ns);
       tSum += t;
     }

     tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
     if(tid == 0){
       tSum = Communicator::instance()->reduce_sum(tSum);
       Out.diff = sqrt(tSum*snorm);
    }
    rr = BGQThread_ScatterDouble(tSum,0,tid,nid);

    //w5 = x;  
    for(int s5=0; s5<N5; ++s5)
      BGWilsonLA_Equate(w5_ptr+s5*Nvol+is, x_ptr+s5*Nvol+is, ns);
  }

  TIMING_END(Out.timing);
  //CCIO::cout << "Kernel section timing: "<< kernel_timing << "\n";
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//// Multishift optimized solver
///////////////////////////////////////////////////////////////////////
void Dirac_DomainWall::solve_ms_init(std::vector<Field>& x,
				     std::vector<Field>& p,
				     Field& r,
				     Field& s,
				     double& rr,
				     std::vector<double>& zeta1,
				     std::vector<double>& zeta2,
				     std::vector<double>& csh2,
				     double& alphap,
				     double& betap) const{
  int Nshift = p.size();
  _Message(SOLV_ITER_VERB_LEVEL, "    MultiShiftSolver_CG Inizialitation\n");
  
  for(int i=0; i<Nshift; ++i){
    p[i] = s;
    x[i] = 0.0;
  }
  r = s;
  rr = r*r;  alphap = 0.0;
  betap  = 1.0;
  _Message(SOLV_ITER_VERB_LEVEL, "    | Initial residual |^2 = "<<rr<<"\n");
}

void Dirac_DomainWall::solve_ms_eo_5d(std::vector<Field>& xq, 
				      const Field& b,
				      SolverOutput& Out, 
				      const std::vector<double>& sigma, 
				      int maxIter, 
				      double targetPrec) const{ 
  using namespace FieldExpression;
 
  BGWilson_DW_Init(prms_.N5_,prms_.mq_,prms_.M0_,
		   (double*)&prms_.dp_[0],(double*)&prms_.dm_[0],
		   (double*)&prms_.bs_[0],(double*)&prms_.cs_[0],
		   (double*)&prms_.es_[0],(double*)&prms_.fs_[0]);

  _Message(SOLV_ITER_VERB_LEVEL, "Multi-shift solver Conjugate Gradient start - BGQ optimized\n");
  Out.Msg = "Multishift CG solver BGQ";
  Out.Iterations = -1;

  TIMING_START;

  int Nshift = sigma.size();
  size_t fsize = b.size();

  double snorm = 1.0/(b.norm()*b.norm());
  
  Field s = b;
  Field r = b;
  std::vector<Field> p(Nshift);
  std::vector<Field> x(Nshift);

  std::vector<double> zeta1(Nshift,1.0);
  std::vector<double> zeta2(Nshift,1.0);
  std::vector<double> csh2(Nshift);
  std::vector<double> pp(Nshift);

  double rr;
  double alphap, betap, rrp;
  int Nshift2 = Nshift;

  Field temp(b.size());

  //////////////////
  BGQThread_Init(); //initializing BGQ fast threading routines
  //note: only sigma[0] is needed...
  
  register int Nvol = CommonPrms::instance()->Nvol()/2;
  register int N5 = prms_.N5_;
  double* u_ptr = const_cast<Field * >(Dw_->getGaugeField_ptr())->getaddr(0);
  Spinor* s_ptr = (Spinor*)s.getaddr(0);
  Spinor* b_ptr = (Spinor*)const_cast<Field&>(b).getaddr(0);
  Spinor* r_ptr = (Spinor*)r.getaddr(0);
  Spinor* p_ptr;
  Spinor* x_ptr;
  Spinor* temp_ptr = (Spinor*)temp.getaddr(0);
  ////////////////////////////////////////

  // Initial messages
  _Message(SOLV_ITER_VERB_LEVEL, "    -------------------\n");
  _Message(SOLV_ITER_VERB_LEVEL, "    Number of shifts = "<< Nshift<<"\n");
  _Message(SOLV_ITER_VERB_LEVEL, "    Values of shifts:\n");
  for(int i = 0; i<Nshift; ++i){
    _Message(SOLV_ITER_VERB_LEVEL, "      #["<<i<<"] = "<< sigma[i]<<"\n");
  }
  _Message(SOLV_ITER_VERB_LEVEL, "    -------------------\n");

  // Initial condition
  for(int i=0; i<Nshift; ++i){
    p[i].resize(fsize);
    x[i].resize(fsize);
    csh2[i] = sigma[i] -sigma[0];
  }
  
  solve_ms_init(x,p,r,s,rr,zeta1,zeta2,csh2,alphap,betap);
  _Message(SOLV_ITER_VERB_LEVEL, "    | Init | = "<<rr*snorm<<"\n");
  
  double kappa = 0.5/(4.0+prms_.M0_);

#pragma omp parallel 
  {
    double timer;
    double t, tSum, rrp;
    double alpha, alphah, beta, pap;
    double zeta, zetas, alphas, betas, zr;
    double residual, diff1;

    int nid = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int ns = Nvol/nid;    
    int is = tid*ns;
  
    double* pp = (double*)BGQThread_Malloc(sizeof(double)*sigma.size(), nid);

    for(int it = 0; it < maxIter; it++){
      if (tid == 0) FINE_TIMING_START(timer);

      ////////////////////////////
      p_ptr = (Spinor*)(p[0].getaddr(0));
      x_ptr = (Spinor*)(x[0].getaddr(0));
        
      BGWilson_DW_Mult_hop(temp_ptr,(void*)(u_ptr),p_ptr,kappa,BGWILSON_DIRAC);
      BGWilson_DW_Mult_hop_dag(s_ptr,(void*)(u_ptr),temp_ptr,kappa,BGWILSON_DIRAC);
      
      //s = opr_->mult(p[0]);
      for(int s5=0; s5<N5; ++s5)
	BGWilsonLA_MultAddScalar(s_ptr+s5*Nvol+is, p_ptr+s5*Nvol+is, sigma[0], ns);
      //s += sigma[0]*p[0];   
      
      ///////////////////////////
      //pap = s*p[0]; c++ code
      tSum = 0.0;
      for(int s5=0; s5<N5; ++s5){
	BGWilsonLA_DotProd(&t, s_ptr+s5*Nvol+is, p_ptr+s5*Nvol+is, ns);
	tSum += t;
      }
      tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
      if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
      
      pap = BGQThread_ScatterDouble(tSum,0,tid,nid);
      ////////////////////////////  
      rrp = rr;
      beta = -rrp/pap;
      
      ///////////////////////////
      //  x[0] -= beta*p[0]; c++ code
      for(int s5=0; s5<N5; ++s5)
	BGWilsonLA_MultAddScalar(x_ptr+s5*Nvol+is, p_ptr+s5*Nvol+is, -beta, ns);
      
      ///////////////////////////
      //  r += beta*s; c++ code
      //  rr = r*r; c++ code
      tSum = 0.0;
      for(int s5=0;s5<N5;++s5){
	BGWilsonLA_MultAddScalar_Norm(r_ptr+s5*Nvol+is, &t, s_ptr+s5*Nvol+is, beta, ns);  
	tSum += t;
      }

      tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
      if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
      rr = BGQThread_ScatterDouble(tSum,0,tid,nid);
      
      alpha = rr/rrp;
      
      ///////////////////////////
      // p[0] *= alpha; p[0] += r;
      for(int s5=0;s5<N5;++s5)
	BGWilsonLA_MultScalar_Add(p_ptr + s5*Nvol+is,r_ptr + s5*Nvol+is,alpha,ns);
      
      pp[0] = rr; 
      
      alphah = 1.0 + alphap*beta/betap;
      
      for(int ish = 1; ish<Nshift2; ++ish){
	p_ptr = (Spinor*)(p[ish].getaddr(0));
	x_ptr = (Spinor*)(x[ish].getaddr(0));
	
	zeta =(alphah-csh2[ish]*beta)/zeta1[ish]+(1.0-alphah)/zeta2[ish];
	zeta = 1.0/zeta;
	zr = zeta/zeta1[ish];
	betas  = beta  *  zr;
	alphas = alpha * zr*zr;
	
	///////////////////////////
	//x[ish] -= betas * p[ish];
	for(int s5=0;s5<N5;++s5)
	  BGWilsonLA_MultAddScalar(x_ptr + s5*Nvol+is, p_ptr + s5*Nvol+is, -betas, ns);
	
	///////////////////////////
	// p[ish] *= alphas;
	// p[ish] += zeta * r;
	// pp[ish] = p[ish] * p[ish];
	for(int s5=0;s5<N5;++s5){
	  BGWilsonLA_MultScalar(p_ptr + s5*Nvol+is, p_ptr + s5*Nvol+is, alphas, ns);
	}
	tSum = 0.0;
	for(int s5=0;s5<N5;++s5){
	  BGWilsonLA_MultAddScalar_Norm(p_ptr+s5*Nvol+is, &t, r_ptr+s5*Nvol+is, zeta, ns);
	  tSum += t;
	}
	
	tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
	if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
       	pp[ish] = BGQThread_ScatterDouble(tSum,0,tid,nid);   
	pp[ish] *= snorm;
	
#pragma omp single 
	{
	  // it is critical to update these variables in the right sequence
	  zeta2[ish] = zeta1[ish];
	  zeta1[ish] = zeta;
	}
      }
      //end of loop
      BGQThread_Barrier(0,nid);    
 
      for(int ish = Nshift2-1; ish>=0; --ish){
	if(pp[ish]> targetPrec){
	  Nshift2 = ish+1;
	  break;
	}
      }
      alphap = alpha;
      betap  = beta;

      ///////////////////////////
      residual = rr*snorm; 
      if(tid==0) {
	_Message(SOLV_ITER_VERB_LEVEL, "   "<<std::setw(5)<<"["<<it<<"]  "
		 <<std::setw(20)<<residual<<"     Left: "<<Nshift2<<" \n");
      }
      if(residual < targetPrec){
	if (tid==0) Out.Iterations = it;
	break;
      }
      if(tid==0) {
	FINE_TIMING_END(timer);
	_Message(DEBUG_VERB_LEVEL,"  Multishift iteration time: "<<timer<<" seconds\n");
      }	
    }//for iter
    BGQThread_Barrier(0,nid);    

    if (tid==0){
      if(Out.Iterations == -1) {
	CCIO::cout << "Not converged.\n";
	exit(1);
      }
      _Message(SOLV_ITER_VERB_LEVEL, "  --- Summary of true residuals\n");
    }
    Out.diff = -1.0;

    double kappa = 0.5/(4.0+prms_.M0_);
    for(int i=0; i<Nshift; ++i){
      x_ptr  = (Spinor*)x[i].getaddr(0);
      BGWilson_DW_Mult_hop(temp_ptr,(void*)(u_ptr),x_ptr,kappa,BGWILSON_DIRAC);
      BGWilson_DW_Mult_hop_dag(s_ptr,(void*)(u_ptr),temp_ptr,kappa,BGWILSON_DIRAC);
      //s = opr_->mult(x[i]);
      //s += sigma[i]*x[i];
      for(int s5=0; s5<N5; ++s5)
	BGWilsonLA_MultAddScalar(s_ptr+s5*Nvol+is, x_ptr+s5*Nvol+is, sigma[i], ns);

      // s -= b;
      for(int s5=0; s5<N5; ++s5)
	BGWilsonLA_Sub(s_ptr+s5*Nvol+is, b_ptr+s5*Nvol+is,ns);

      // diff1 = s * s;
      tSum = 0.0;
      for(int s5=0; s5<N5; ++s5){
	BGWilsonLA_Norm(&t,s_ptr+s5*Nvol+is,ns);
	tSum += t;
      }
      tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
      if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
      diff1 = BGQThread_ScatterDouble(tSum,0,tid,nid);
      diff1 *= snorm;

      if(tid==0) _Message(SOLV_ITER_VERB_LEVEL,"       ["<<i<<"]  "<<diff1<<"\n");
      if(diff1>Out.diff) Out.diff = diff1;
    }
    if(tid==0)
      _Message(SOLV_ITER_VERB_LEVEL," Maximum residual  = "<<Out.diff<<"\n");
    
    for(int i=0; i<Nshift; ++i){
      //xq[i] = x[i];
      x_ptr  = (Spinor*)x[i].getaddr(0);
      Spinor* xqi_ptr = (Spinor*)(xq[i].getaddr(0));
      for(int s5=0; s5<N5; ++s5)
	BGWilsonLA_Equate(xqi_ptr+s5*Nvol+is, x_ptr+s5*Nvol+is, ns);
    }
    BGQThread_Free(pp,tid);
  }
  Out.diff = sqrt(Out.diff);
  TIMING_END(Out.timing);
}

