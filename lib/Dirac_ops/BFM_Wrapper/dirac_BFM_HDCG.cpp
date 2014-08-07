/*!
 * @file dirac_BFM_HDCG_wrapper.cpp
 * @brief Defines the wrapper class methods for P. Boyle HDCG inverter
 * Time-stamp: <2014-08-07 16:44:13 neo>
 */
#include <math.h>
#include "dirac_BFM_HDCG.hpp"

void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_rotate(double *phases_by_pi,double *phases_dir, int sgn){
  
  /* TwistedBC does the following
     for(int i=0; i < Nd-1; i++) {
     int direction = phases_dir[i]; 
     Real arg = phases_by_pi[i]*onepi / Real( nrow[ direction ] );
     Real re = cos(arg);
     Real im = sin(arg);
     Complex phase = cmplx( re, im );
     u[ direction ] *= phase;
     }
     * Suppose 
     *    D psi = 0
     * Want to modify subspace vector psi so that
     *    D^tw psi^tw = 0
     * If 
     *    psi^tw(x) = e^{-i x arg} psi(x)
     * Then 
     *    U_\mu^tw psi^tw(x+mu) = e^{i x arg} U_\mu psi(x+mu) 
     */
  linop.BossLog("HDCG Rotating subspace for twisted BCs\n");
  

  double onepi = M_PI;
  int LocalVol = CommonPrms::instance()->Nvol();
  int Ls=linop.Ls;  

  vector_double theta, eitheta_re, eitheta_im;// scalar fields
  theta.resize(LocalVol); 
  eitheta_re.resize(LocalVol); 
  eitheta_im.resize(LocalVol); 
  for(int i=0; i < ND_-1; i++) { 
    int mu = phases_dir[i];
    
    for(int site = 0; site < LocalVol; site++) {
	int coord[4];
	coord[0] = SiteIndex::instance()->c_x(site);
	coord[1] = SiteIndex::instance()->c_y(site);
	coord[2] = SiteIndex::instance()->c_z(site);
	coord[3] = SiteIndex::instance()->c_t(site);
	
	theta[site] = theta[site] + (SiteIndex::instance()->*SiteIndex::global_idx[mu])(coord[mu]) * phases_by_pi[i]* onepi / CommonPrms::instance()->global_size(mu) ;

	eitheta_re[site] = cos(theta[site]);
	eitheta_im[site] = -sgn*sin(theta[site]);
	
	
      }
	  }
    
    FermionField copy(LocalVol*Ls); // 5d object
  
    for(int v=0;v<ldop_d->Nvec;v++){
      BFM_interface.FermionExport_to_BFM(copy, ldop_d->subspace_d[v], 1); //last parameter Even/Odd
      //linop.exportFermion(copy,ldop_d->subspace_d[v],1);
      for(int s=0;s<Ls;s++){
	for(int site = 0; site < LocalVol; site++) {
	    for (int n_in = 0; n_in < copy.Nin()/2; n_in ++){
	      int idx = copy.format.index(n_in, site+LocalVol*s);
	      int idxp = copy.format.index(n_in+1, site+LocalVol*s);
	      copy.data.set(idx, eitheta_re[site]*copy[idx]);
	      copy.data.set(idxp, eitheta_im[site]*copy[idxp]);
	    }
	  }
	      }
	      BFM_interface.FermionImport_from_BFM(copy, ldop_d->subspace_d[v], 1); //last parameter Even/Odd     
	    //linop.importFermion(copy,ldop_d->subspace_d[v],1);
	    
	
      }
  
    
  }


void Dirac_BFM_HDCG_Wrapper::HDCG_set_mass(double mass){
    linop.BossLog("HDCG setting mass %le\n",mass);
    linop.GeneralisedFiveDimEnd();
    linop.mass = mass;
    linop.GeneralisedFiveDimInit();

    linop_single.GeneralisedFiveDimEnd();
    linop_single.mass = mass;
    linop_single.GeneralisedFiveDimInit();

    linop_r.GeneralisedFiveDimEnd();
    linop_r.mass = mass;
    linop_r.GeneralisedFiveDimInit();

}


void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_init(){
  //    int sloppy_relax=Parms.SubspaceRationalRefine;
  int sloppy_relax=1;
  
  uint32_t seed = ldop_d->ncoor[0]*17 + ldop_d->ncoor[1]*251;
  + ldop_d->ncoor[2]*1023  + ldop_d->ncoor[3]*8191;
  srand48(seed);
  
  if ( sloppy_relax ) { 
    ldop_d->RelaxSubspace<float>(&linop_r);
  } else {
    ldop_d->RelaxSubspace<double>(&linop);
  }
  ldop_d->SinglePrecSubspace();
}

void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_refine(){
  
  int Nvec = HDCGParms.NumberSubspace;
  
  Fermion_t src = linop.allocFermion();
  Fermion_t sol = linop.allocFermion();
  Fermion_t tmp = linop.allocFermion();
  
  int threads = linop.threads;
  BfmHDCG<float>  *ldop_f;
  ldop_f = new BfmHDCG<float>(HDCGParms.Ls,HDCGParms.NumberSubspace,HDCGParms.Block,HDCGParms.Block,&linop,&linop_single);
  ldop_f->CloneSubspace<double>(*ldop_d);
  ldop_f->SetParams(HDCGParms);

  ldop_f->PcgSingleShift = HDCGParms.SubspaceRationalRefineLo;
  linop.residual = HDCGParms.SubspaceRationalRefineResidual;
  ldop_f->PcgType        = PcgAdef2fSingleShift;
  
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      for(int v=0;v<Nvec;v++){
	
	double nn=linop.norm(ldop_d->subspace_r[v]);
	linop.ThreadBossLog("HDCG Refining vector[%d] of norm %le\n",v,nn);
	
	linop.copy(src,ldop_d->subspace_r[v]);
	linop.set_zero(sol);
	ldop_f->Pcg(ldop_d,sol,src,tmp);
	linop.copy(ldop_d->subspace_r[v],sol);
	
	nn=linop.norm(ldop_d->subspace_r[v]);
	//	  linop.ThreadBossLog(" Refined vector[%d] new norm %le\n",v,nn);
      }
      for(int v=0;v<Nvec;v++){
	linop.copy(ldop_d->subspace_d[v],ldop_d->subspace_r[v]);
      }
    }
  }
  ldop_f->end();
  delete ldop_f;
  
  linop.freeFermion(tmp);
  linop.freeFermion(src);
  linop.freeFermion(sol);
  
  ldop_d->FreeSingleSubspace();
  ldop_d->OrthogonaliseSubspace();
  ldop_d->SinglePrecSubspace();
  
}

void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_free(void){
  ldop_d->FreeDoubleSubspace();
  ldop_d->FreeSingleSubspace();
}


void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_compute(int sloppy){
  
  int threads = linop.threads;
  int Nvec = HDCGParms.NumberSubspace;
  
  if ( sloppy ) { 
    ldop_d->ComputeLittleMatrixColored<float>(&linop_single,ldop_d->subspace_f);
    ldop_d->LdopDeflationBasisSize=0;
  } else { 
    //ldop_d->ComputeLittleMatrixColored<float>(&dop_f,ldop_d->subspace_f);
    ldop_d->ComputeLittleMatrixColored<double>(&linop,ldop_d->subspace_d);
    // NB due to BG/Q QPX have hard coded DeflationBasisSize must be multiple of 4.
    // This may become even more on future machines.
#if 1
    if( HDCGParms.SubspaceSurfaceDepth < HDCGParms.Ls/2) {
      linop.BossMessage("Zeroing s=[%d,%d]\n",HDCGParms.SubspaceSurfaceDepth,HDCGParms.Ls-HDCGParms.SubspaceSurfaceDepth-1);
#pragma omp parallel 
      {
#pragma omp for 
	for(int i=0;i<threads;i++) {
	  for(int v=0;v<Nvec;v++){
	    double nno = linop.norm(ldop_d->subspace_d[v]);
	    for(int s=HDCGParms.SubspaceSurfaceDepth; s<HDCGParms.Ls-HDCGParms.SubspaceSurfaceDepth;s++){
	      linop.axpby_ssp(ldop_d->subspace_d[v],0.0,ldop_d->subspace_d[v],0.0,ldop_d->subspace_d[v],s,s);
	    }
	    double nnr = linop.norm(ldop_d->subspace_d[v]);
	      linop.ThreadBossLog("Restricted vec[%d] %le -> %le \n",v,nno,nnr);
	  }
	}
      }
    }
#endif
    
    uint32_t seed = ldop_d->ncoor[0]*11 + ldop_d->ncoor[1]*161
      + ldop_d->ncoor[2]*357 + ldop_d->ncoor[3]*1203;
    srand48(seed);
    
    int LdopDeflVecs = HDCGParms.LdopDeflVecs;
    ldop_d->LdopDeflationBasisSize=0;
    //      ldop_d->LdopDeflationBasisTrivial();
    int min =4;
      if ( ldop_d->LdopDeflationBasisSize> min )  min = ldop_d->LdopDeflationBasisSize;
      for (int v = min ; v<=LdopDeflVecs; v+=4){
	ldop_d->LdopDeflationBasisInit(v);
      }
      ldop_d->LdopDeflationBasisDiagonalise(LdopDeflVecs);
      
  }
}

void Dirac_BFM_HDCG_Wrapper::solve_HDCG(Fermion_t solution[2], Fermion_t source[2], double residual, int maxit){
  
  linop.residual = residual;
  linop.max_iter = maxit;
  
  Fermion_t src = linop.allocFermion();// also this may be useless - use psi_h in the parent class
  //Fermion_t tmp = linop.allocFermion(); // already in the definition of the parent class
  Fermion_t Mtmp= linop.allocFermion();
  Fermion_t resid= linop.allocFermion();

  // Even/Odd eventually must be changed... Peter has a different convention


  BfmHDCG<float>  *ldop_f;
  ldop_f = new BfmHDCG<float>(HDCGParms.Ls,HDCGParms.NumberSubspace,HDCGParms.Block,HDCGParms.Block,&linop,&linop_single);
  ldop_f->CloneSubspace<double>(*ldop_d);
  ldop_f->SetParams(HDCGParms);
  ldop_f->PcgType=PcgAdef2f;

  int threads = linop.threads;
#pragma omp parallel
  {
#pragma omp for
    for(int i=0;i<threads;i++) {

      int me = linop.thread_barrier();

      // checkboardsource
      linop.MooeeInv(source[Even],tmp,DaggerNo,Even);
      linop.Meo     (tmp,src,Odd,DaggerNo);
      linop.axpy    (tmp,src,source[Odd],-1.0);
      linop.Mprec(tmp,src,Mtmp,DaggerYes);

      ldop_f->AnalyseSpectrum=0;
      if ( HDCGParms.Flexible == 1 ) { 
	linop_single.ThreadBossLog("Solving with single precision Ldop, fPcg\n");
	ldop_f->fPcg(ldop_d,solution[Odd],src,tmp);
      } else if ( HDCGParms.Flexible == 2 ) {  
	linop.ThreadBossLog("Solving with double precision Ldop, fPcg\n");
	ldop_d->fPcg(ldop_d,solution[Odd],src,tmp);
      } else { 
	linop_single.ThreadBossLog("Solving with single precision Ldop, Pcg\n");
	ldop_f->Pcg(ldop_d,solution[Odd],src,tmp);
      }

      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      linop.Meo(solution[Odd],tmp,Even,DaggerNo);
      linop.axpy(src,tmp,source[Even],-1.0);
      linop.MooeeInv(src,solution[Even],DaggerNo,Even);

      // Even defect
      linop.Meo(solution[Odd],tmp,Even,DaggerNo);
      linop.Mooee(solution[Even],Mtmp,DaggerNo,Even);
      linop.axpy(Mtmp,Mtmp,tmp,1.0);
      linop.axpy(Mtmp,Mtmp,source[Even],-1.0);

      double rf=linop.norm(Mtmp);
      // Odd defect
      linop.Meo(solution[Even],tmp,Odd,DaggerNo);
      linop.Mooee(solution[Odd],Mtmp,DaggerNo,Odd);
      linop.axpy(Mtmp,Mtmp,tmp,1.0);
      linop.axpy(Mtmp,Mtmp,source[Odd],-1.0);
      rf+=linop.norm(Mtmp);

      double f = linop.norm(source[Odd]);
      f       += linop.norm(source[Even]);

      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: unprec sol true residual is %le\n",sqrt(rf/f));

    }
  }

  linop.freeFermion(src);
  //linop.freeFermion(tmp);
  linop.freeFermion(Mtmp);
  linop.freeFermion(resid);

  ldop_f->end();
  delete ldop_f;

}    


void Dirac_BFM_HDCG_Wrapper::solve_HDCG(FermionField &sol, FermionField &src,double residual, int maxit){
  // Wrapper of the previous basic method for the inverter
}
