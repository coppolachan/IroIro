/*!@file  gaugeFixing_Landau.cpp 
   @brief
   @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru) 
   LastChangedBy: noaki
*/
#include "gaugeFixing_Landau.hpp"

using namespace std;

#include "Main/commonPrms.h"
#include "Communicator/communicator.h"
#include "Tools/randNum.hpp"
#include "Tools/sunMatUtils.hpp"

void GaugeFixing_Landau::fix(GaugeField& Uout, const GaugeField& Uin){

  GaugeField Ue = Uin.even_sec();
  GaugeField Uo = Uin.odd_sec();

  int Nvh = Uin.Nvol()/2;
  GaugeField1D  Ge(Nvh), Go(Nvh);

  int Nconv = -1;

  // gauge fixing iteration
  for(int iter = 0; iter < d_Niter; ++iter){

    if((iter % d_Nmeas) == 0){
      double sg, Fval;
      calc_SG(sg,Fval,Ue,Uo);
      pprintf("  iter = %6d  sg = %16.8e  Fval = %16.8e\n",
              iter,sg,Fval);
      if(sg < d_Enorm){
        Nconv = iter;
        pprintf("converged at iter = %d\n",Nconv);
        break;
      }
    }

    double wp2 = d_wp;
    if((iter % d_Nreset) < d_Nnaive) wp2 = 1.0;
    gfix_step(Ue, Uo, wp2);

    if((iter % d_Nreset) == 0 && iter > 0){
      // random gauge transformation
      pprintf("  random gauge transformation performed.\n");
      set_randomGtrf(Ge);
      gauge_trf(Ue, Uo, Ge, 0);
      set_randomGtrf(Go);
      gauge_trf(Ue, Uo, Go, 1);
    }

  }

  d_index.reverseField(Ufix,Ue,0);
  d_index.reverseField(Ufix,Uo,1);

}
//====================================================================
void GaugeFixing_Landau::gfix_step(Field_G& Ue, Field_G& Uo,
                                                        double wp){

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

    double wp_sbt = 1.0-wp;
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
void GaugeFixing_Landau::set_randomGtrf(Field_G& Geo){

  int Nvol = Geo.nvol();
  int Nex  = Geo.nex();

  int Nc = CommonPrms::Nc();
  SUNmat gt(Nc);

  for(int ex = 0; ex < Nex; ++ex){
   for(int site = 0; site < Nvol; ++site){
     //gt = Geo.mat(site,ex);
     gt.set_random(d_rnd);
     Geo.set_mat(site,ex, gt);
   }
  }

}
//====================================================================
void GaugeFixing_Landau::gauge_trf(Field_G& Ue, Field_G& Uo,
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
void GaugeFixing_Landau::calc_SG(double& sg, double& Fval,
                                         Field_G& Ue, Field_G& Uo){

  int Nc = CommonPrms::Nc();
  int NPE = CommonPrms::NPE();
  int Nvol2 = Ue.nvol();
  int Nex   = Ue.nex();

  Field_G DLT(Nvol2,1);
  SUNmat ut(Nc);

  sg   = 0.0;
  Fval = 0.0;

  for(int ieo = 0; ieo < 2; ++ieo){
    calc_DLT(DLT, Ue, Uo, ieo);
    double tsg = DLT.norm2();
    sg += tsg;
  }
  sg = sg/(Nex*Nc*2*Nvol2*NPE);

  for(int mu = 0; mu < Nex; ++mu){
   for(int site = 0; site < Nvol2; ++site){
     ut = Ue.mat(site,mu);
     Fval += ReTr(ut);
     ut = Uo.mat(site,mu);
     Fval += ReTr(ut);
    }
  }
  Fval = Communicator::reduce_sum(Fval);
  Fval = Fval/(Nex*2*Nvol2*NPE);

}
//====================================================================
void GaugeFixing_Landau::calc_DLT(Field_G& DLT,
                                  Field_G& Ue, Field_G& Uo, int Ieo){

  int Nvol2 = Ue.nvol();
  int Nc = CommonPrms::Nc();
  int Ndim = CommonPrms::Ndim();

  ShiftField_eo shift;

  Field_G Ut1(Nvol2,1), Ut2(Nvol2,1);
  SUNmat u_tmp(Nc);

  DLT = 0.0;

  if(Ieo == 0){ // on even sites

    for(int mu = 0; mu < Ndim; ++mu){
      DLT.addpart_ex(0,Ue,mu,-1.0);
      Ut1.setpart_ex(0,Uo,mu);
      shift.forward_h(Ut2,Ut1,mu,0);
      DLT.addpart_ex(0,Ut2,0);
    }

  }else{        // on odd sites

    for(int mu = 0; mu < Ndim; ++mu){
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
void GaugeFixing_Landau::calc_W(Field_G& Weo,
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

    for(int mu = 0; mu < Ndim; ++mu){
      Weo.addpart_ex(0,Ue,mu);
      Ut1.setpart_ex(0,Uo,mu);
      shift.forward_h(Ut2,Ut1,mu,0);
      for(int site = 0; site < Nvol2; ++site){
        u_tmp = Ut2.mat_dag(site,0);
        Weo.add_mat(site,0,u_tmp);
      }
    }

  }else if(Ieo == 1){ // on odd sites

    for(int mu = 0; mu < Ndim; ++mu){
      Weo.addpart_ex(0,Uo,mu);
      Ut1.setpart_ex(0,Ue,mu);
      shift.forward_h(Ut2,Ut1,mu,1);
      for(int site = 0; site < Nvol2; ++site){
        u_tmp = Ut2.mat_dag(site,0);
        Weo.add_mat(site,0,u_tmp);
      }
    }

  }else{
    pprintf("  Wrong ieo: gaugeFixing_Landau::calcW().\n");
    abort();
  }

}
//====================================================================
void GaugeFixing_Landau::maxTr(Field_G& G0, Field_G& W){
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
void GaugeFixing_Landau::maxTr1(Field_G& G, Field_G& W){

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
void GaugeFixing_Landau::maxTr2(Field_G& G, Field_G& W){

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
void GaugeFixing_Landau::maxTr3(Field_G& G, Field_G& W){

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
