//----------------------------------------------------------------------
// dirac_wilson.cpp
//----------------------------------------------------------------------
#include "dirac_wilson.h"
#include "Tools/sunMatUtils.hpp"

using namespace SUNvec_utils;
using namespace std;

void Dirac_Wilson::mult_xp(Field& fp, ShiftField* sfp) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result
  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),0));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data is contiguous

 
    if (!sfp->on_bdry(site)) {vtmp = sfp->get_bulk_addr(site);}
    else {vtmp = sfp->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] - vtmp[2*c+6*NC_+1 ];
      v1tmp[c][1] = vtmp[2*c+1      ] + vtmp[2*c+6*NC_   ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] - vtmp[2*c+4*NC_+1 ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+4*NC_   ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {

	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v1tmp[c1][0] 
		   - utmp[NC_*2*c+2*c1+1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v1tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v1tmp[c1][1]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v2tmp[c1][0] 
                   - utmp[NC_*2*c+2*c1+1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v2tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v2tmp[c1][1]);
      }

      res[2*c        ] += v1[c][0];
      res[2*c+1      ] += v1[c][1];
      res[2*c+2*NC_  ] += v2[c][0];
      res[2*c+2*NC_+1] += v2[c][1];
      res[2*c+4*NC_  ] += v2[c][1];
      res[2*c+4*NC_+1] -= v2[c][0];
      res[2*c+6*NC_  ] += v1[c][1];
      res[2*c+6*NC_+1] -= v1[c][0];

    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, ShiftField* sfp) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),1));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data is contiguous

    //gamma_0
    if (!sfp->on_bdry(site)) {vtmp = sfp->get_bulk_addr(site);}
    else {vtmp = sfp->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] + vtmp[2*c+6*NC_   ];
      v1tmp[c][1] = vtmp[2*c+1      ] + vtmp[2*c+6*NC_+1 ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] - vtmp[2*c+4*NC_   ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] - vtmp[2*c+4*NC_+1 ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v1tmp[c1][0] 
		   - utmp[NC_*2*c+2*c1+1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v1tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v1tmp[c1][1]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v2tmp[c1][0] 
                   - utmp[NC_*2*c+2*c1+1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v2tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v2tmp[c1][1]);
      }

      res[2*c        ] += v1[c][0];
      res[2*c+1      ] += v1[c][1];
      res[2*c+2*NC_  ] += v2[c][0];
      res[2*c+2*NC_+1] += v2[c][1];
      res[2*c+4*NC_  ] -= v2[c][0];
      res[2*c+4*NC_+1] -= v2[c][1];
      res[2*c+6*NC_  ] += v1[c][0];
      res[2*c+6*NC_+1] += v1[c][1];

    }

  }
}

void Dirac_Wilson::mult_zp(Field& fp, ShiftField* sfp) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),2));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data is contiguous

    if (!sfp->on_bdry(site)) {vtmp = sfp->get_bulk_addr(site);}
    else {vtmp = sfp->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] - vtmp[2*c+4*NC_+1 ];
      v1tmp[c][1] = vtmp[2*c+1      ] + vtmp[2*c+4*NC_   ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] + vtmp[2*c+6*NC_+1 ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] - vtmp[2*c+6*NC_   ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v1tmp[c1][0] 
		   - utmp[NC_*2*c+2*c1+1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v1tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v1tmp[c1][1]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v2tmp[c1][0] 
                   - utmp[NC_*2*c+2*c1+1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v2tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v2tmp[c1][1]);
      }

      res[2*c        ] += v1[c][0];
      res[2*c+1      ] += v1[c][1];
      res[2*c+2*NC_  ] += v2[c][0];
      res[2*c+2*NC_+1] += v2[c][1];
      res[2*c+4*NC_  ] += v1[c][1];
      res[2*c+4*NC_+1] -= v1[c][0];
      res[2*c+6*NC_  ] -= v2[c][1];
      res[2*c+6*NC_+1] += v2[c][0];

    }

  }
}

void Dirac_Wilson::mult_tp(Field& fp, ShiftField* sfp) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),3));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data is contiguous

    if (!sfp->on_bdry(site)) {vtmp = sfp->get_bulk_addr(site);}
    else {vtmp = sfp->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c+4*NC_  ]*2.0;
      v1tmp[c][1] = vtmp[2*c+4*NC_+1]*2.0;
      v2tmp[c][0] = vtmp[2*c+6*NC_  ]*2.0;
      v2tmp[c][1] = vtmp[2*c+6*NC_+1]*2.0;
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v1tmp[c1][0] 
		   - utmp[NC_*2*c+2*c1+1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v1tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v1tmp[c1][1]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v2tmp[c1][0] 
                   - utmp[NC_*2*c+2*c1+1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v2tmp[c1][0] 
                   + utmp[NC_*2*c+2*c1  ]*v2tmp[c1][1]);
      }

      res[2*c+4*NC_  ] += v1[c][0];
      res[2*c+4*NC_+1] += v1[c][1];
      res[2*c+6*NC_  ] += v2[c][0];
      res[2*c+6*NC_+1] += v2[c][1];

    }

  }
}

void Dirac_Wilson::mult_xm(valarray<double>& w, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site),0));
    res  = &w[ff_->index_r(0,0,site)];
    //assumes matrix and fermion data is contiguous
    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
    
     for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] + vtmp[2*c+6*NC_+1 ];
      v1tmp[c][1] = vtmp[2*c+1      ] - vtmp[2*c+6*NC_   ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] + vtmp[2*c+4*NC_+1 ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] - vtmp[2*c+4*NC_   ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
	            + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	v1[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	v2[c][0] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
                    + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	v2[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);

      }
    
      res[2*c        ] =  v1[c][0];
      res[2*c+1      ] =  v1[c][1];
      res[2*c+2*NC_  ] =  v2[c][0];
      res[2*c+2*NC_+1] =  v2[c][1];
      res[2*c+4*NC_  ] = -v2[c][1];
      res[2*c+4*NC_+1] =  v2[c][0];
      res[2*c+6*NC_  ] = -v1[c][1];
      res[2*c+6*NC_+1] =  v1[c][0];

    }
  }
}

void Dirac_Wilson::mult_ym(valarray<double>& w, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site),1));
    res  = &w[ff_->index_r(0,0,site)];
    //assumes matrix and fermion data is contiguous
    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
    
     for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] - vtmp[2*c+6*NC_   ];
      v1tmp[c][1] = vtmp[2*c+1      ] - vtmp[2*c+6*NC_+1 ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] + vtmp[2*c+4*NC_   ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+4*NC_+1 ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
	            + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	v1[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	v2[c][0] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
                    + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	v2[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);

      }
    
      res[2*c        ] =  v1[c][0];
      res[2*c+1      ] =  v1[c][1];
      res[2*c+2*NC_  ] =  v2[c][0];
      res[2*c+2*NC_+1] =  v2[c][1];
      res[2*c+4*NC_  ] =  v2[c][0];
      res[2*c+4*NC_+1] =  v2[c][1];
      res[2*c+6*NC_  ] = -v1[c][0];
      res[2*c+6*NC_+1] = -v1[c][1];

    }
  }
}


void Dirac_Wilson::mult_zm(valarray<double>& w, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;

  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site),2));
    res  = &w[ff_->index_r(0,0,site)];
    //assumes matrix and fermion data is contiguous
    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
    
     for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] + vtmp[2*c+4*NC_+1 ];
      v1tmp[c][1] = vtmp[2*c+1      ] - vtmp[2*c+4*NC_   ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] - vtmp[2*c+6*NC_+1 ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+6*NC_   ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
	            + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	v1[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	v2[c][0] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
                    + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	v2[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);

      }
    
      res[2*c        ] =  v1[c][0];
      res[2*c+1      ] =  v1[c][1];
      res[2*c+2*NC_  ] =  v2[c][0];
      res[2*c+2*NC_+1] =  v2[c][1];
      res[2*c+4*NC_  ] = -v1[c][1];
      res[2*c+4*NC_+1] =  v1[c][0];
      res[2*c+6*NC_  ] =  v2[c][1];
      res[2*c+6*NC_+1] = -v2[c][0];

    }
  }
}

void Dirac_Wilson::mult_tm(valarray<double>& w, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;

  double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site),3));
    res  = &w[ff_->index_r(0,0,site)];
    //assumes matrix and fermion data is contiguous
    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
    
     for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ]*2.0;
      v1tmp[c][1] = vtmp[2*c+1      ]*2.0;
      v2tmp[c][0] = vtmp[2*c+2*NC_  ]*2.0;
      v2tmp[c][1] = vtmp[2*c+2*NC_+1]*2.0;
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
	            + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	v1[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	v2[c][0] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
                    + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	v2[c][1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
                    - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);

      }
    
      res[2*c        ] =  v1[c][0];
      res[2*c+1      ] =  v1[c][1];
      res[2*c+2*NC_  ] =  v2[c][0];
      res[2*c+2*NC_+1] =  v2[c][1];
      res[2*c+4*NC_  ] =  0.0;
      res[2*c+4*NC_+1] =  0.0;
      res[2*c+6*NC_  ] =  0.0;
      res[2*c+6*NC_+1] =  0.0;

    }
  }
}

void (Dirac_Wilson::*Dirac_Wilson::mult_p[])
(Field&,ShiftField*) const = {&Dirac_Wilson::mult_xp,
			      &Dirac_Wilson::mult_yp,
			      &Dirac_Wilson::mult_zp,
			      &Dirac_Wilson::mult_tp,};

void (Dirac_Wilson::*Dirac_Wilson::mult_m[])
(valarray<double>&,const Field&) const = {&Dirac_Wilson::mult_xm,
					  &Dirac_Wilson::mult_ym,
					  &Dirac_Wilson::mult_zm,
					  &Dirac_Wilson::mult_tm,};

const Field Dirac_Wilson::gamma5(const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();
  Field w(fsize_);
  for(int site = 0; site<Nvol_; ++site){
    for (int c = 0; c <Nc; ++c) {
      w.set(ff_->index_r(c,0,site),f[ff_->index_r(c,2,site)]);
      w.set(ff_->index_i(c,0,site),f[ff_->index_i(c,2,site)]);
      w.set(ff_->index_r(c,1,site),f[ff_->index_r(c,3,site)]);
      w.set(ff_->index_i(c,1,site),f[ff_->index_i(c,3,site)]);
      w.set(ff_->index_r(c,2,site),f[ff_->index_r(c,0,site)]);
      w.set(ff_->index_i(c,2,site),f[ff_->index_i(c,0,site)]);
      w.set(ff_->index_r(c,3,site),f[ff_->index_r(c,1,site)]);
      w.set(ff_->index_i(c,3,site),f[ff_->index_i(c,1,site)]);
    }
    //    w.set(ff_->cslice(0,site),f[ff_->cslice(2,site)]);
    //    w.set(ff_->cslice(1,site),f[ff_->cslice(3,site)]);
    //    w.set(ff_->cslice(2,site),f[ff_->cslice(0,site)]);
    //    w.set(ff_->cslice(3,site),f[ff_->cslice(1,site)]);
  }
  return w;
}

const Field Dirac_Wilson::proj_p(const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();
  Field w(fsize_);
  for(int site = 0; site<Nvol_; ++site){
    for (int c = 0; c <Nc; ++c) {
      double fup_r = 0.5*(f[ff_->index_r(c,0,site)]+f[ff_->index_r(c,2,site)]);
      double fup_i = 0.5*(f[ff_->index_i(c,0,site)]+f[ff_->index_i(c,2,site)]);
      double fdn_r = 0.5*(f[ff_->index_r(c,1,site)]+f[ff_->index_r(c,3,site)]);
      double fdn_i = 0.5*(f[ff_->index_i(c,1,site)]+f[ff_->index_i(c,3,site)]);
      w.set(ff_->index_r(c,0,site),fup_r);
      w.set(ff_->index_i(c,0,site),fup_i);
      w.set(ff_->index_r(c,1,site),fdn_r);
      w.set(ff_->index_i(c,1,site),fdn_i);
      w.set(ff_->index_r(c,2,site),fup_r);
      w.set(ff_->index_i(c,2,site),fup_i);
      w.set(ff_->index_r(c,3,site),fdn_r);
      w.set(ff_->index_i(c,3,site),fdn_i);
    }
  }
  return w;
}

const Field Dirac_Wilson::proj_m(const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();
  Field w(fsize_);
  for(int site = 0; site<Nvol_; ++site){
    for (int c = 0; c <Nc; ++c) {
      double fup_r = 0.5*(f[ff_->index_r(c,0,site)]-f[ff_->index_r(c,2,site)]);
      double fup_i = 0.5*(f[ff_->index_i(c,0,site)]-f[ff_->index_i(c,2,site)]);
      double fdn_r = 0.5*(f[ff_->index_r(c,1,site)]-f[ff_->index_r(c,3,site)]);
      double fdn_i = 0.5*(f[ff_->index_i(c,1,site)]-f[ff_->index_i(c,3,site)]);
      w.set(ff_->index_r(c,0,site),fup_r);
      w.set(ff_->index_i(c,0,site),fup_i);
      w.set(ff_->index_r(c,1,site),fdn_r);
      w.set(ff_->index_i(c,1,site),fdn_i);
      w.set(ff_->index_r(c,2,site),-fup_r);
      w.set(ff_->index_i(c,2,site),-fup_i);
      w.set(ff_->index_r(c,3,site),-fdn_r);
      w.set(ff_->index_i(c,3,site),-fdn_i);
    }
  }
  return w;
}
/*
void Dirac_Wilson::mult_core(Field& w,Field& f) const{

  valarray<double> wt(fsize_);

  for(int d=0; d <Ndim_; ++d){
    sf_up_[d]->setf(f);
    (this->*mult_p[d])(w,sf_up_[d]);
    (this->*mult_m[d])(wt,f);
    sf_dn_[d]->setf(wt);
    w += sf_dn_[d]->getva();
  }
  w *= -kpp_;
}
*/

void Dirac_Wilson::mult_a0(Field& w,Field& f) const{
  valarray<double> wt(fsize_);

  for(int d=0; d <Ndim_; ++d){
    sf_up_[d]->setf(f);
    (this->*mult_p[d])(w,sf_up_[d]);
    (this->*mult_m[d])(wt,f);
    sf_dn_[d]->setf(wt);
    w += sf_dn_[d]->getva();
  }
  w *= -kpp_;
}

void Dirac_Wilson::mult_a1(Field& w,Field& f) const{
  mult_a0(w,f);
  w += f;
}

const Field Dirac_Wilson::mult(const Field& f) const{
  Field w(fsize_);
  (this->*mult_core)(w,const_cast<Field&>(f));
  return w;
}

const Field Dirac_Wilson::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

/*!
 *  @brief MD-force contribution: \f$\zeta^\dagger\frac{dH_W}{d\tau}\eta\f$
 */
void Dirac_Wilson::md_force_p(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmat_utils;

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();
  SUNmat f;

  for(int mu=0; mu<Ndim_; ++mu){
    Field xie(fsize_);
  
    sf_up_[mu]->setf(const_cast<Field&>(eta));
    (this->*mult_p[mu])(xie, sf_up_[mu]);

    for(int site=0; site<Nvol_; ++site){
      f = 0.0;;
      for(int a=0; a<Nc; ++a){
        for(int b=0; b<Nc; ++b){
          double fre = 0.0;
          double fim = 0.0;
          for(int s=0; s<Nd; ++s){

	    size_t ra =ff_->index_r(a,s,site);
	    size_t ia =ff_->index_i(a,s,site);

	    size_t rb =ff_->index_r(b,s,site);
	    size_t ib =ff_->index_i(b,s,site);

	    fre += zeta[rb]*xie[ra] +zeta[ib]*xie[ia];
	    fim += zeta[rb]*xie[ia] -zeta[ib]*xie[ra];
          }
          f.set(a,b,fre,fim);
        }
      }
      // fce.add(gf_->cslice(0,gauge_site_p(site),mu),anti_hermite(f));
      int gsite = (this->*gp)(site);
      fce.add(gf_->cslice(0,gsite,mu),anti_hermite(f));
    }
  }
}

void Dirac_Wilson::md_force_m(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmat_utils;

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();
  SUNmat f;
  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

  for(int mu=0; mu<Ndim_; ++mu){
    Field xz5(fsize_);

    sf_up_[mu]->setf(const_cast<Field&>(zt5));
    (this->*mult_p[mu])(xz5, sf_up_[mu]);

    for(int site=0; site<Nvol_; ++site){
      f=0.0;
      for(int a=0; a<Nc; ++a){
        for(int b=0; b<Nc; ++b){
          double fre = 0.0;
          double fim = 0.0;
          for(int s=0; s<Nd; ++s){

	    size_t ra =ff_->index_r(a,s,site);
	    size_t ia =ff_->index_i(a,s,site);

	    size_t rb =ff_->index_r(b,s,site);
	    size_t ib =ff_->index_i(b,s,site);

	    fre -= xz5[rb]*et5[ra] +xz5[ib]*et5[ia];
	    fim -= xz5[rb]*et5[ia] -xz5[ib]*et5[ra];
          }
          f.set(a,b,fre,fim);
        }
      }
      
      int gsite = (this->*gp)(site);
      fce.add(gf_->cslice(0,gsite,mu),anti_hermite(f));
    }
  }
}

const Field Dirac_Wilson::md_force(const Field& eta,const Field& zeta)const{

  Field fp(gf_->size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);

  fp *= -kpp_;
  return fp;
}

const vector<int> Dirac_Wilson::get_gsite() const {
  return SiteIndex::instance()->get_gsite();
}

namespace Dw{
  double read_mass(const XML::node& node){
    double mass;
    XML::read(node, "mass", mass);
    return mass;
  }
}
