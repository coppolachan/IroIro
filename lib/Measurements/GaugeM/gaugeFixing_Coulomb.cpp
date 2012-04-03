/*!
        @file    $Id:: gaugeFixing_Coulomb.cpp #$

        @brief

        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru) 
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2012-03-15 17:01:15 #$

        @version $LastChangedRevision: 137 $
*/

#include "gaugeFixing_Coulomb.h"

using namespace std;

#include "common_prms.h"
#include "communicator.h"
using CommunicatorItems::pprintf;

#include "shiftField_eo.h"
#include "randomNumbers.h"
#include "sun_mat.h"
using namespace SUNmat_utils;


//====================================================================
void GaugeFixing_Coulomb::set_prms(int Niter, int Nnaive, int Nmeas,
                                int Nreset, double Enorm, double wp){

  d_Niter  = Niter;
  d_Nnaive = Nnaive;
  d_Nmeas  = Nmeas;
  d_Nreset = Nreset;
  d_Enorm  = Enorm;
  d_wp     = wp;

  pprintf("Coulomb gauge fixing:\n");
  pprintf("  Niter  = %d\n",d_Niter);
  pprintf("  Nnaive = %d\n",d_Nnaive);
  pprintf("  Nmeas  = %d\n",d_Nmeas);
  pprintf("  Nreset = %d\n",d_Nreset);
  pprintf("  Enorm  = %12.4e\n",d_Enorm);
  pprintf("  wp     = %8.4f\n",d_wp);

}

//====================================================================
void GaugeFixing_Coulomb::set_prms(){
  // This function is just for test, and remained to give an idea
  // about typical values of paramters.

  d_Niter  = 5000;
  d_Nnaive = 50;
  d_Nmeas  = 10;
  d_Nreset = 1000;
  d_Enorm  = 1.0e-24;
  d_wp     = 1.6;

  pprintf("Coulomb gauge fixing:\n");
  pprintf("  Niter  = %d\n",d_Niter);
  pprintf("  Nnaive = %d\n",d_Nnaive);
  pprintf("  Nmeas  = %d\n",d_Nmeas);
  pprintf("  Nreset = %d\n",d_Nreset);
  pprintf("  Enorm  = %12.4e\n",d_Enorm);
  pprintf("  wp     = %8.4f\n",d_wp);

}

//====================================================================
void GaugeFixing_Coulomb::fix(Field_G& Ufix, const Field_G& Uorg){

  int Nvol = Uorg.nvol();
  int Nex  = Uorg.nex();
  int Lt   = CommonPrms::Lt();

  int Nvol2 = Nvol/2;
  Field_G Ue(Nvol2,Nex), Uo(Nvol2,Nex);

  Field_G Ge(Nvol2,1), Go(Nvol2,1);

  d_index.convertField(Ue,Uorg,0);
  d_index.convertField(Uo,Uorg,1);

  int Nconv = -1;

  valarray<double> sg(Lt);
  valarray<double> Fval(Lt);

  // gauge fixing iteration
  for(int iter = 0; iter < d_Niter; ++iter){

    if((iter % d_Nmeas) == 0){

      calc_SG(sg,Fval,Ue,Uo);

      double sg_max = 0.0;
      pprintf("  iter = %6d:\n",iter);
      for(int t = 0; t < Lt; ++t){
         pprintf("    t = %4d  sg = %16.8e  Fval = %16.8e\n",
                  t,sg[t],Fval[t]);
         if(sg[t] > sg_max) sg_max = sg[t];
      }

      if(sg_max < d_Enorm){
        Nconv = iter;
        pprintf("converged at iter = %d\n",Nconv);
        break;
      }

    }

    double wp2 = d_wp;
    if((iter % d_Nreset) < d_Nnaive) wp2 = 1.0;
    gfix_step(sg, Ue, Uo, wp2);

    if((iter % d_Nreset) == 0 && iter > 0){
      // random gauge transformation
      pprintf("  random gauge transformation performed.\n");
      set_randomGtrf(sg,Ge);
      gauge_trf(Ue, Uo, Ge, 0);
      set_randomGtrf(sg,Go);
      gauge_trf(Ue, Uo, Go, 1);
    }

  }

  d_index.reverseField(Ufix,Ue,0);
  d_index.reverseField(Ufix,Uo,1);

}
//====================================================================
void GaugeFixing_Coulomb::gfix_step(valarray<double>& sg,
                               Field_G& Ue, Field_G& Uo, double wp){

  int Nc = CommonPrms::Nc();
  int Nvol2 = Ue.nvol();
  int Nex   = Ue.nex();

  Field_G Weo(Nvol2,1), Geo(Nvol2,1);
  SUNmat ut(Nc), uwp(Nc);
  uwp.unit();
  uwp *= (1.0-wp);

  for(int ieo = 0; ieo < 2; ++ieo){

    calc_W(Weo, Ue, Uo, ieo);
    maxTr(Geo,Weo);
    Geo *= wp;

    //double wp_sbt = 1.0-wp;
    for(int site = 0; site < Nvol2; ++site){
      ut = Geo.mat(site,0);
      ut += uwp;
      ut.reunit();
      Geo.set_mat(site,0,ut);
    }

    gauge_trf(Ue, Uo, Geo, ieo);

  }

}
//====================================================================
void GaugeFixing_Coulomb::set_randomGtrf(valarray<double>& sg,
                                                     Field_G& Geo){

  int Lt = CommonPrms::Lt();
  int Nt = CommonPrms::Nt();
  int Nz = CommonPrms::Nz();
  int Ny = CommonPrms::Ny();
  int Nx2 = CommonPrms::Nx()/2;
  int Nc = CommonPrms::Nc();
  int ipet = Communicator::ipe(3);

  int Nvol = Geo.nvol();
  assert(Geo.nex() == 1);
  assert(sg.size() == Lt);

  SUNmat gt(Nc);

  for(int t = 0; t < Nt; ++t){
   int tg = t + ipet*Nt;

   if(sg[tg] > d_Enorm){
     for(int z = 0; z < Nz; ++z){
      for(int y = 0; y < Ny; ++y){
       for(int x = 0; x < Nx2; ++x){
         int site = d_index.siteh(x,y,z,t);
         gt.set_random(d_rnd);
         Geo.set_mat(site,0, gt);
       }
      }
     }
   }else{
     gt.unit();
     for(int z = 0; z < Nz; ++z){
      for(int y = 0; y < Ny; ++y){
       for(int x = 0; x < Nx2; ++x){
         int site = d_index.siteh(x,y,z,t);
         Geo.set_mat(site,0, gt);
       }
      }
     }
   }

  }

}
//====================================================================
void GaugeFixing_Coulomb::gauge_trf(Field_G& Ue, Field_G& Uo,
                                            Field_G& Geo, int Ieo){
  //  Ieo = 0: gauge transformation on even sites.
  //  Ieo = 1:                      on odd sites.

  int Nvol2 = Geo.nvol();
  int Nex   = Geo.nex();
  int Ndim  = CommonPrms::Ndim();

  ShiftField_eo shift;

  Field_G Ut(Nvol2,1), Gt(Nvol2,1);
  Field_G Ut2(Nvol2,1);

  if(Ieo == 0){
    for(int mu = 0; mu < Ndim; ++mu){
      Ut.mult_Field_Gnn(0,Geo,0,Ue,mu);
      Ue.setpart_ex(mu,Ut,0);

      shift.backward_h(Gt,Geo,mu,1);
      Ut.mult_Field_Gnd(0,Uo,mu,Gt,0);
      Uo.setpart_ex(mu,Ut,0);
    }
  }else{
    for(int mu = 0; mu < Ndim; ++mu){
      Ut.mult_Field_Gnn(0,Geo,0,Uo,mu);
      Uo.setpart_ex(mu,Ut,0);

      shift.backward_h(Gt,Geo,mu,0);
      Ut.mult_Field_Gnd(0,Ue,mu,Gt,0);
      Ue.setpart_ex(mu,Ut,0);
    }
  }

}
//====================================================================
void GaugeFixing_Coulomb::calc_SG(valarray<double>& sg,
                                  valarray<double>& Fval,
                                  Field_G& Ue, Field_G& Uo){

  int Nc = CommonPrms::Nc();
  int Lt = CommonPrms::Lt();
  int Nt = CommonPrms::Nt();
  int Nz = CommonPrms::Nz();
  int Ny = CommonPrms::Ny();
  int Nx2 = CommonPrms::Nx()/2;
  int Nvol2 = Nx2*Ny*Nz*Nt;
  int Ndim  = CommonPrms::Ndim();
  int NPE  = CommonPrms::NPE();

  assert(Ue.nex()  == Ndim);
  assert(Ue.nvol() == Nvol2);
  assert(Uo.nex()  == Ndim);
  assert(Uo.nvol() == Nvol2);

  valarray<double> sg_local(Nt);
  valarray<double> Fval_local(Nt);
  Field_G DLT(Nvol2,1);
  SUNmat ut(Nc);

  sg   = 0.0;
  Fval = 0.0;
  sg_local   = 0.0;
  Fval_local = 0.0;

  for(int ieo = 0; ieo < 2; ++ieo){
    calc_DLT(DLT, Ue, Uo, ieo);

    for(int t = 0; t < Nt; ++t){
     for(int z = 0; z < Nz; ++z){
      for(int y = 0; y < Ny; ++y){
       for(int x = 0; x < Nx2; ++x){
         int site = d_index.siteh(x,y,z,t);
         ut = DLT.mat(site,0);
         sg_local[t] += ut.norm2();
       }
      }
     }
    }
  }

  sum_global_t(sg,sg_local);
  for(int t = 0; t < Lt; ++t){
    sg[t] = sg[t]/(Ndim*Nc*(2*Nvol2*NPE)/Lt);
  }

  for(int mu = 0; mu < Ndim-1; ++mu){
    for(int t = 0; t < Nt; ++t){
     for(int z = 0; z < Nz; ++z){
      for(int y = 0; y < Ny; ++y){
       for(int x = 0; x < Nx2; ++x){
         int site = d_index.siteh(x,y,z,t);
         ut = Ue.mat(site,mu);
         Fval_local[t] += ReTr(ut);
         ut = Uo.mat(site,mu);
         Fval_local[t] += ReTr(ut);
       }
      }
     }
    }
  }

  sum_global_t(Fval,Fval_local);
  for(int t = 0; t < Lt; ++t){
    Fval[t] = Fval[t]/(Ndim*(2*Nvol2*NPE)/Lt);
  }

}
//====================================================================
void GaugeFixing_Coulomb::sum_global_t(valarray<double>& val_global,
                                       valarray<double>& val_local){

  int Lt = CommonPrms::Lt();
  int Nt = CommonPrms::Nt();
  int ipet = Communicator::ipe(3);

  assert(val_global.size() == Lt);
  assert(val_local.size()  == Nt);

  for(int t = 0; t< Lt; ++t){
    val_global[t] = 0.0;
  }

  for(int tl = 0; tl< Nt; ++tl){
    val_global[tl + ipet*Nt] = val_local[tl];
  }

  for(int t = 0; t< Lt; ++t){
    double val = val_global[t];
    val_global[t] = Communicator::reduce_sum(val);
  }

}
//====================================================================
void GaugeFixing_Coulomb::calc_DLT(Field_G& DLT,
                                  Field_G& Ue, Field_G& Uo, int Ieo){

  int Nvol2 = Ue.nvol();
  int Nc = CommonPrms::Nc();
  int Ndim = CommonPrms::Ndim();

  ShiftField_eo shift;

  Field_G Ut1(Nvol2,1), Ut2(Nvol2,1);
  SUNmat u_tmp(Nc);

  DLT = 0.0;

  if(Ieo == 0){ // on even sites

    for(int mu = 0; mu < Ndim-1; ++mu){
      DLT.addpart_ex(0,Ue,mu,-1.0);
      Ut1.setpart_ex(0,Uo,mu);
      shift.forward_h(Ut2,Ut1,mu,0);
      DLT.addpart_ex(0,Ut2,0);
    }

  }else{        // on odd sites

    for(int mu = 0; mu < Ndim-1; ++mu){
      DLT.addpart_ex(0,Uo,mu,-1.0);
      Ut1.setpart_ex(0,Ue,mu);
      shift.forward_h(Ut2,Ut1,mu,1);
      DLT.addpart_ex(0,Ut2,0);
    }

  }

  for(int site = 0; site < Nvol2; ++site){
    u_tmp = DLT.mat(site,0);
    u_tmp.at();
    u_tmp *=2.0;
    DLT.set_mat(site,0,u_tmp);
  }

}

//====================================================================
void GaugeFixing_Coulomb::calc_W(Field_G& Weo,
                                Field_G& Ue, Field_G& Uo, int Ieo){

  int Nvol2 = Ue.nvol();
  int Nc = CommonPrms::Nc();
  int Ndim = CommonPrms::Ndim();

  assert(Weo.nex() == 1);

  ShiftField_eo shift;

  Field_G Ut1(Nvol2,1), Ut2(Nvol2,1);
  SUNmat u_tmp(Nc);

  Weo = 0.0;

  if(Ieo == 0){       // on even sites

    for(int mu = 0; mu < Ndim-1; ++mu){
      Weo.addpart_ex(0,Ue,mu);
      Ut1.setpart_ex(0,Uo,mu);
      shift.forward_h(Ut2,Ut1,mu,0);
      for(int site = 0; site < Nvol2; ++site){
        u_tmp = Ut2.mat_dag(site,0);
        Weo.add_mat(site,0,u_tmp);
      }
    }

  }else if(Ieo == 1){ // on odd sites

    for(int mu = 0; mu < Ndim-1; ++mu){
      Weo.addpart_ex(0,Uo,mu);
      Ut1.setpart_ex(0,Ue,mu);
      shift.forward_h(Ut2,Ut1,mu,1);
      for(int site = 0; site < Nvol2; ++site){
        u_tmp = Ut2.mat_dag(site,0);
        Weo.add_mat(site,0,u_tmp);
      }
    }

  }else{
    pprintf("  Wrong ieo: gaugeFixing_Coulomb::calcW().\n");
    abort();
  }

}
//====================================================================
void GaugeFixing_Coulomb::maxTr(Field_G& G0, Field_G& W){
  // Present implementation only applys to SU(3) case.

  int Nc = CommonPrms::Nc();
  int Nvol2 = G0.nvol();

  int Nmt = 1;

  SUNmat unity(Nc);
  unity.unit();

  for(int site = 0; site < Nvol2; ++site){
    G0.set_mat(site,0,unity);
  }

  for(int imt = 0; imt < Nmt; ++imt){
    maxTr1(G0, W);
    maxTr2(G0, W);
    maxTr3(G0, W);
  }

}
//====================================================================
void GaugeFixing_Coulomb::maxTr1(Field_G& G, Field_G& W){

  int Nc = CommonPrms::Nc();
  int Nvol2 = W.nvol();

  SUNmat gt(Nc), wt(Nc), gt2(Nc), wt2(Nc);

  for(int site = 0; site < Nvol2; ++site){

    wt = W.mat(site,0);

    gt.set(2, 0.0, 0.0);
    gt.set(5, 0.0, 0.0);
    gt.set(6, 0.0, 0.0);
    gt.set(7, 0.0, 0.0);
    gt.set(8, 1.0, 0.0);

    double fn1 =  (wt.r(0) + wt.r(4))*(wt.r(0) + wt.r(4))
                + (wt.i(0) - wt.i(4))*(wt.i(0) - wt.i(4));
    double fn2 =  (wt.r(1) - wt.r(3))*(wt.r(1) - wt.r(3))
                + (wt.i(1) + wt.i(3))*(wt.i(1) + wt.i(3));
    double fn = 1.0/sqrt(fn1+fn2);

    gt.set( 0, fn*( wt.r(0)+wt.r(4)), fn*(-wt.i(0)+wt.i(4)) );
    gt.set( 1, fn*(-wt.r(1)+wt.r(3)), fn*(-wt.i(1)-wt.i(3)) );
    gt.set( 3, fn*( wt.r(1)-wt.r(3)), fn*(-wt.i(1)-wt.i(3)) );
    gt.set( 4, fn*( wt.r(0)+wt.r(4)), fn*( wt.i(0)-wt.i(4)) );

    wt2 = gt * wt;
    W.set_mat(site,0,wt2);
    gt2 = G.mat(site,0);
    wt2 = gt * gt2;
    G.set_mat(site,0,wt2);

  }

}
//====================================================================
void GaugeFixing_Coulomb::maxTr2(Field_G& G, Field_G& W){

  int Nc = CommonPrms::Nc();
  int Nvol2 = W.nvol();

  SUNmat gt(Nc), wt(Nc), gt2(Nc), wt2(Nc);

  for(int site = 0; site < Nvol2; ++site){

    wt = W.mat(site,0);

    gt.set(1, 0.0, 0.0);
    gt.set(3, 0.0, 0.0);
    gt.set(4, 1.0, 0.0);
    gt.set(5, 0.0, 0.0);
    gt.set(7, 0.0, 0.0);

    double fn1 =  (wt.r(8) + wt.r(0))*(wt.r(8) + wt.r(0))
                + (wt.i(8) - wt.i(0))*(wt.i(8) - wt.i(0));
    double fn2 =  (wt.r(2) - wt.r(6))*(wt.r(2) - wt.r(6))
                + (wt.i(2) + wt.i(6))*(wt.i(2) + wt.i(6));
    double fn = 1.0/sqrt(fn1+fn2);

    gt.set( 0, fn*( wt.r(8)+wt.r(0)), fn*( wt.i(8)-wt.i(0)) );
    gt.set( 2, fn*( wt.r(6)-wt.r(2)), fn*(-wt.i(6)-wt.i(2)) );
    gt.set( 6, fn*(-wt.r(6)+wt.r(2)), fn*(-wt.i(6)-wt.i(2)) );
    gt.set( 8, fn*( wt.r(8)+wt.r(0)), fn*(-wt.i(8)+wt.i(0)) );

    wt2 = gt * wt;
    W.set_mat(site,0,wt2);
    gt2 = G.mat(site,0);
    wt2 = gt * gt2;
    G.set_mat(site,0,wt2);

  }

}
//====================================================================
void GaugeFixing_Coulomb::maxTr3(Field_G& G, Field_G& W){

  int Nc = CommonPrms::Nc();
  int Nvol2 = W.nvol();

  SUNmat gt(Nc), wt(Nc), gt2(Nc), wt2(Nc);

  for(int site = 0; site < Nvol2; ++site){

    wt = W.mat(site,0);

    gt.set(0, 1.0, 0.0);
    gt.set(1, 0.0, 0.0);
    gt.set(2, 0.0, 0.0);
    gt.set(3, 0.0, 0.0);
    gt.set(6, 0.0, 0.0);

    double fn1 =  (wt.r(4) + wt.r(8))*(wt.r(4) + wt.r(8))
                + (wt.i(4) - wt.i(8))*(wt.i(4) - wt.i(8));
    double fn2 =  (wt.r(7) - wt.r(5))*(wt.r(7) - wt.r(5))
                + (wt.i(7) + wt.i(5))*(wt.i(7) + wt.i(5));
    double fn = 1.0/sqrt(fn1+fn2);

    gt.set( 4, fn*( wt.r(4)+wt.r(8)), fn*(-wt.i(4)+wt.i(8)) );
    gt.set( 5, fn*(-wt.r(5)+wt.r(7)), fn*(-wt.i(5)-wt.i(7)) );
    gt.set( 7, fn*( wt.r(5)-wt.r(7)), fn*(-wt.i(5)-wt.i(7)) );
    gt.set( 8, fn*( wt.r(4)+wt.r(8)), fn*( wt.i(4)-wt.i(8)) );

    wt2 = gt * wt;
    W.set_mat(site,0,wt2);
    gt2 = G.mat(site,0);
    wt2 = gt * gt2;
    G.set_mat(site,0,wt2);

  }

}
//====================================================================
//============================================================END=====
