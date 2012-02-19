//----------------------------------------------------------------------
// dirac_wilson.cpp
//----------------------------------------------------------------------
#include "dirac_wilson.h"
#include "Tools/sunMatUtils.hpp"

using namespace SUNvec_utils;
using namespace std;

#ifdef IMPROVED_WILSON

void Dirac_Wilson::mult_xp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num of elements of a half spinor */
  int Nbdry = bdry_plw_[0].size();
  int Nbulk = bulk_pup_[0].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; /*!< @brief information on the lower boundary */   
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_plw_[0][k]));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r0(c)]-v[i3(c)]; vbd[i0(c)+is] = v[i0(c)]+v[r3(c)];
      vbd[r1(c)+is] = v[r1(c)]-v[i2(c)]; vbd[i1(c)+is] = v[i1(c)]+v[r2(c)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry]; /*!< @brief information on the upper neighbor */   
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,0);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  for(int k=0; k<Nbdry; ++k){   /*!< @brief calc on the upper boundary */   
    const double* U = u_->getaddr(gf_->index(0,(this->*gp)(bdry_pup_[0][k]),0));
    int id = ff_.index(0,bdry_pup_[0][k]);    
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r0(c), v1r[c]);   fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);   fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c), v2i[c]);   fp.add(id+i2(c),-v2r[c]);
      fp.add(id+r3(c), v1i[c]);   fp.add(id+i3(c),-v1r[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  for(int k=0; k<Nbulk; ++k){   /*!< @brief calc on the bulk */   
    const double* v = f.getaddr(ff_.index(0,(this->*x_p)(bulk_pup_[0][k],0)));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] -v[i3(c)];  w1i[c] = v[i0(c)] +v[r3(c)];
      w2r[c] = v[r1(c)] -v[i2(c)];  w2i[c] = v[i1(c)] +v[r2(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bulk_pup_[0][k]),0));
    int id = ff_.index(0,bulk_pup_[0][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(id+r0(c), v1r[c]);  fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);  fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c), v2i[c]);  fp.add(id+i2(c),-v2r[c]);
      fp.add(id+r3(c), v1i[c]);  fp.add(id+i3(c),-v1r[c]);
    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_plw_[1].size();
  int Nbulk = bulk_pup_[1].size();

  /// boundary part ///
  double vbd[Nih*Nbdry];
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_plw_[1][k]));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r0(c)]+v[r3(c)]; vbd[i0(c)+is] = v[i0(c)]+v[i3(c)];
      vbd[r1(c)+is] = v[r1(c)]-v[r2(c)]; vbd[i1(c)+is] = v[i1(c)]-v[i2(c)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,1);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bdry_pup_[1][k]),1));
    int id = ff_.index(0,bdry_pup_[1][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r0(c), v1r[c]);   fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);   fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c),-v2r[c]);   fp.add(id+i2(c),-v2i[c]);
      fp.add(id+r3(c), v1r[c]);   fp.add(id+i3(c), v1i[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int k=0; k<Nbulk; ++k){
    const double* v = f.getaddr(ff_.index(0,(this->*x_p)(bulk_pup_[1][k],1)));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[r3(c)];  w1i[c] = v[i0(c)] +v[i3(c)];
      w2r[c] = v[r1(c)] -v[r2(c)];  w2i[c] = v[i1(c)] -v[i2(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bulk_pup_[1][k]),1));
    int id = ff_.index(0,bulk_pup_[1][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(id+r0(c), v1r[c]);  fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);  fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c),-v2r[c]);  fp.add(id+i2(c),-v2i[c]);
      fp.add(id+r3(c), v1r[c]);  fp.add(id+i3(c), v1i[c]);
    }
  }
}

void Dirac_Wilson::mult_zp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_plw_[2].size();
  int Nbulk = bulk_pup_[2].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_plw_[2][k]));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r0(c)]-v[i2(c)]; vbd[i0(c)+is] = v[i0(c)]+v[r2(c)];
      vbd[r1(c)+is] = v[r1(c)]+v[i3(c)]; vbd[i1(c)+is] = v[i1(c)]-v[r3(c)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,2);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bdry_pup_[2][k]),2));
    int id = ff_.index(0,bdry_pup_[2][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r0(c), v1r[c]);   fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);   fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c), v1i[c]);   fp.add(id+i2(c),-v1r[c]);
      fp.add(id+r3(c),-v2i[c]);   fp.add(id+i3(c), v2r[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int k=0; k<Nbulk; ++k){
    const double* v = f.getaddr(ff_.index(0,(this->*x_p)(bulk_pup_[2][k],2)));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] -v[i2(c)];  w1i[c] = v[i0(c)] +v[r2(c)];
      w2r[c] = v[r1(c)] +v[i3(c)];  w2i[c] = v[i1(c)] -v[r3(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bulk_pup_[2][k]),2));
    int id = ff_.index(0,bulk_pup_[2][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(id+r0(c), v1r[c]);  fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);  fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c), v1i[c]);  fp.add(id+i2(c),-v1r[c]);
      fp.add(id+r3(c),-v2i[c]);  fp.add(id+i3(c), v2r[c]);
    }
  }
}

void Dirac_Wilson::mult_tp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_plw_[3].size();
  int Nbulk = bulk_pup_[3].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k) {
    const double* v = f.getaddr(ff_.index(0,bdry_plw_[3][k]));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r2(c)]*2.0;  vbd[i0(c)+is] = v[i2(c)]*2.0;
      vbd[r1(c)+is] = v[r3(c)]*2.0;  vbd[i1(c)+is] = v[i3(c)]*2.0;
    }
    is += Nih;
  }  
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,3);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bdry_pup_[3][k]),3));
    int id = ff_.index(0,bdry_pup_[3][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r2(c), v1r[c]);   fp.add(id+i2(c), v1i[c]);
      fp.add(id+r3(c), v2r[c]);   fp.add(id+i3(c), v2i[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int k=0; k<Nbulk; ++k) {
    const double* v = f.getaddr(ff_.index(0,(this->*x_p)(bulk_pup_[3][k],3)));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r2(c)]*2.0;  w1i[c] = v[i2(c)]*2.0;
      w2r[c] = v[r3(c)]*2.0;  w2i[c] = v[i3(c)]*2.0;
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gp)(bulk_pup_[3][k]),3));
    int id = ff_.index(0,bulk_pup_[3][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
	
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(id+r2(c), v1r[c]);  fp.add(id+i2(c), v1i[c]);
      fp.add(id+r3(c), v2r[c]);  fp.add(id+i3(c), v2i[c]);
    }
  }
}

void Dirac_Wilson::mult_xm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[0].size();
  int Nbulk = bulk_mlw_[0].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  /// boundary part ///
  double vbd[Nih*Nbdry]; /*!< @brief information on the upper boundary */
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_mup_[0][k]));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[i3(c)];  w1i[c] = v[i0(c)] -v[r3(c)];
      w2r[c] = v[r1(c)] +v[i2(c)];  w2i[c] = v[i1(c)] -v[r2(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(bdry_mup_[0][k]),0));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,0);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int id = ff_.index(0,bdry_mlw_[0][k]);
    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]);   fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]);   fm.add(id+i1(c), vbc[i1(c)+is]);
      fm.add(id+r2(c),-vbc[i1(c)+is]);   fm.add(id+i2(c), vbc[r1(c)+is]);
      fm.add(id+r3(c),-vbc[i0(c)+is]);   fm.add(id+i3(c), vbc[r0(c)+is]);
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; 
  for(int k=0; k<Nbulk; ++k){
    int xm = (this->*x_m)(bulk_mlw_[0][k],0);
    const double* v = f.getaddr(ff_.index(0,xm));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[i3(c)];  w1i[c] = v[i0(c)] -v[r3(c)];
      w2r[c] = v[r1(c)] +v[i2(c)];  w2i[c] = v[i1(c)] -v[r2(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(xm),0));
    int id = ff_.index(0,bulk_mlw_[0][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(id+r0(c), v1r[c]); 
      fm.add(id+i0(c), v1i[c]);
      fm.add(id+r1(c), v2r[c]);
      fm.add(id+i1(c), v2i[c]);
      fm.add(id+r2(c),-v2i[c]);
      fm.add(id+i2(c), v2r[c]);
      fm.add(id+r3(c),-v1i[c]); 
      fm.add(id+i3(c), v1r[c]);
    }
  }
}

void Dirac_Wilson::mult_ym(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[1].size();
  int Nbulk = bulk_mlw_[1].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  // boundary part  
  double vbd[Nih*Nbdry];
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_mup_[1][k]));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] -v[r3(c)];  w1i[c] = v[i0(c)] -v[i3(c)];
      w2r[c] = v[r1(c)] +v[r2(c)];  w2i[c] = v[i1(c)] +v[i2(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(bdry_mup_[1][k]),1));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;
	
      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy v1 from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,1);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int id = ff_.index(0,bdry_mlw_[1][k]);
    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]);   fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]);   fm.add(id+i1(c), vbc[i1(c)+is]);
      fm.add(id+r2(c), vbc[r1(c)+is]);   fm.add(id+i2(c), vbc[i1(c)+is]);
      fm.add(id+r3(c),-vbc[r0(c)+is]);   fm.add(id+i3(c),-vbc[i0(c)+is]);
    }
    is += Nih;
  }
  //bulk part
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];
  for(int k=0; k<Nbulk; ++k){
    int ym = (this->*x_m)(bulk_plw_[1][k],1);
    const double* v = f.getaddr(ff_.index(0,ym));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] -v[r3(c)];  w1i[c] = v[i0(c)] -v[i3(c)];
      w2r[c] = v[r1(c)] +v[r2(c)];  w2i[c] = v[i1(c)] +v[i2(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(ym),1));
    int id = ff_.index(0,bulk_mlw_[1][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
      fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
      fm.add(id+r2(c), v2r[c]);  fm.add(id+i2(c), v2i[c]);
      fm.add(id+r3(c),-v1r[c]);  fm.add(id+i3(c),-v1i[c]); 
    }
  }
}

void Dirac_Wilson::mult_zm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[2].size();
  int Nbulk = bulk_mlw_[2].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  // boundary part
  double vbd[Nih*Nbdry];
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_mup_[2][k]));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[i2(c)];  w1i[c] = v[i0(c)] -v[r2(c)];
      w2r[c] = v[r1(c)] -v[i3(c)];  w2i[c] = v[i1(c)] +v[r3(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(bdry_mup_[2][k]),2));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;    vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;    vbd[i1(c)+is] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];   //Copy v1 from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,2);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int id = ff_.index(0,bdry_mlw_[2][k]);
    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]);   fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]);   fm.add(id+i1(c), vbc[i1(c)+is]);
      fm.add(id+r2(c),-vbc[i0(c)+is]);   fm.add(id+i2(c), vbc[r0(c)+is]);
      fm.add(id+r3(c), vbc[i1(c)+is]);   fm.add(id+i3(c),-vbc[r1(c)+is]);
    }
    is += Nih;
  }
  //bulk part
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];
  for(int k=0; k<Nbulk; ++k){
    int zm = (this->*x_m)(bulk_mlw_[2][k],2);
    const double* v = f.getaddr(ff_.index(0,zm));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[i2(c)];  w1i[c] = v[i0(c)] -v[r2(c)];
      w2r[c] = v[r1(c)] -v[i3(c)];  w2i[c] = v[i1(c)] +v[r3(c)];
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(zm),2));
    int id = ff_.index(0,bulk_mlw_[2][k]);
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
      fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
      fm.add(id+r2(c),-v1i[c]);  fm.add(id+i2(c), v1r[c]);
      fm.add(id+r3(c), v2i[c]);  fm.add(id+i3(c),-v2r[c]); 
    }
  }
}

void Dirac_Wilson::mult_tm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[3].size();
  int Nbulk = bulk_mlw_[3].size();
  
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  /// boundary part ///
  double vbd[Nih*Nbdry];
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    const double* v = f.getaddr(ff_.index(0,bdry_mup_[3][k]));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)]*2.0;   w1i[c] = v[i0(c)]*2.0;
      w2r[c] = v[r1(c)]*2.0;   w2i[c] = v[i1(c)]*2.0;
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(bdry_mup_[3][k]),3));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	vbd[is+r0(c)] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[is+i0(c)] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[is+r1(c)] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[is+i1(c)] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,3);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int id = ff_.index(0,bdry_mlw_[3][k]);
    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]); fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]); fm.add(id+i1(c), vbc[i1(c)+is]);
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];
  for(int k=0; k<Nbulk; ++k){
    int tm = (this->*x_m)(bulk_mlw_[3][k],3);
    const double* v = f.getaddr(ff_.index(0,tm));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)]*2.0;   w1i[c] = v[i0(c)]*2.0;
      w2r[c] = v[r1(c)]*2.0;   w2i[c] = v[i1(c)]*2.0;
    }
    const double* U = u_->getaddr(gf_.index(0,(this->*gm)(tm),3));
    int id = ff_.index(0,bulk_mlw_[3][k]);       
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0; 
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
      fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
    }
  }
}

const Field Dirac_Wilson::gamma5(const Field& f) const{
  Field w(fsize_);
  
  for(int site=0; site<Nvol_; ++site){
    int id = ff_.index(0,site);
    const double* ft = f.getaddr(id);

    for(int c=0; c <NC_; ++c){
      w.set(id+r0(c), ft[r2(c)]);   w.set(id+i0(c), ft[i2(c)]);
      w.set(id+r1(c), ft[r3(c)]);   w.set(id+i1(c), ft[i3(c)]);
      w.set(id+r2(c), ft[r0(c)]);   w.set(id+i2(c), ft[i0(c)]);
      w.set(id+r3(c), ft[r1(c)]);   w.set(id+i3(c), ft[i1(c)]);
    }
  }
  return w;
}

const Field Dirac_Wilson::proj_p(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    int id = ff_.index(0,site);
    const double* ft = f.getaddr(id);

    for(int c=0; c<NC_; ++c){
      double fup_r = 0.5*(ft[r0(c)] +ft[r2(c)]);
      double fup_i = 0.5*(ft[i0(c)] +ft[i2(c)]);
      double fdn_r = 0.5*(ft[r1(c)] +ft[r3(c)]);
      double fdn_i = 0.5*(ft[i1(c)] +ft[i3(c)]);

      w.set(id+r0(c), fup_r);   w.set(id+i0(c), fup_i);
      w.set(id+r1(c), fdn_r);   w.set(id+i1(c), fdn_i);
      w.set(id+r2(c), fup_r);   w.set(id+i2(c), fup_i);
      w.set(id+r3(c), fdn_r);   w.set(id+i3(c), fdn_i);
    }
  }
  return w;
}

const Field Dirac_Wilson::proj_m(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    int id = ff_.index(0,site);
    const double* ft = f.getaddr(id);

    for (int c=0; c<NC_; ++c){
      double fup_r = 0.5*(ft[r0(c)] -ft[r2(c)]);
      double fup_i = 0.5*(ft[i0(c)] -ft[i2(c)]);
      double fdn_r = 0.5*(ft[r1(c)] -ft[r3(c)]);
      double fdn_i = 0.5*(ft[i1(c)] -ft[i3(c)]);

      w.set(id+r0(c), fup_r);   w.set(id+i0(c), fup_i);
      w.set(id+r1(c), fdn_r);   w.set(id+i1(c), fdn_i);
      w.set(id+r2(c),-fup_r);   w.set(id+i2(c),-fup_i);
      w.set(id+r3(c),-fdn_r);   w.set(id+i3(c),-fdn_i);
    }
  }
  return w;
}

#endif /*IMPROVED_WILSON*/

#ifndef IMPROVED_WILSON

void Dirac_Wilson::mult_xp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; 
  int Nbdry = bdry_plw_[0].size();
  int Nbulk = bulk_pup_[0].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_plw_[0][k];
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = f[ff_.index_r(c,0,site)] -f[ff_.index_i(c,3,site)];
      vbd[i0(c)+is] = f[ff_.index_i(c,0,site)] +f[ff_.index_r(c,3,site)];
      vbd[r1(c)+is] = f[ff_.index_r(c,1,site)] -f[ff_.index_i(c,2,site)];
      vbd[i1(c)+is] = f[ff_.index_i(c,1,site)] +f[ff_.index_r(c,2,site)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry]; 
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,0);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];   
  double U[2*NC_*NC_];
  is = 0;
  for(int k=0; k<Nbdry; ++k){  
    int site = bdry_pup_[0][k];    

    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),0)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),0)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(ff_.index_r(c,0,site), v1r[c]);
      fp.add(ff_.index_i(c,0,site), v1i[c]);
      fp.add(ff_.index_r(c,1,site), v2r[c]);
      fp.add(ff_.index_i(c,1,site), v2i[c]);
      fp.add(ff_.index_r(c,2,site), v2i[c]);
      fp.add(ff_.index_i(c,2,site),-v2r[c]);
      fp.add(ff_.index_r(c,3,site), v1i[c]);
      fp.add(ff_.index_i(c,3,site),-v1r[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  for(int k=0; k<Nbulk; ++k){

    int site = bulk_pup_[0][k];
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),0)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),0)];
      }
    }
    int xp = (this->*x_p)(site,0);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,xp)] -f[ff_.index_i(c,3,xp)];
      w1i[c] = f[ff_.index_i(c,0,xp)] +f[ff_.index_r(c,3,xp)];
      w2r[c] = f[ff_.index_r(c,1,xp)] -f[ff_.index_i(c,2,xp)];
      w2i[c] = f[ff_.index_i(c,1,xp)] +f[ff_.index_r(c,2,xp)];
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(ff_.index_r(c,0,site), v1r[c]);
      fp.add(ff_.index_i(c,0,site), v1i[c]);
      fp.add(ff_.index_r(c,1,site), v2r[c]);
      fp.add(ff_.index_i(c,1,site), v2i[c]);
      fp.add(ff_.index_r(c,2,site), v2i[c]);
      fp.add(ff_.index_i(c,2,site),-v2r[c]);
      fp.add(ff_.index_r(c,3,site), v1i[c]);
      fp.add(ff_.index_i(c,3,site),-v1r[c]);
    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; 
  int Nbdry = bdry_plw_[1].size();
  int Nbulk = bulk_pup_[1].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_plw_[1][k];
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = f[ff_.index_r(c,0,site)] +f[ff_.index_r(c,3,site)];
      vbd[i0(c)+is] = f[ff_.index_i(c,0,site)] +f[ff_.index_i(c,3,site)];
      vbd[r1(c)+is] = f[ff_.index_r(c,1,site)] -f[ff_.index_r(c,2,site)];
      vbd[i1(c)+is] = f[ff_.index_i(c,1,site)] -f[ff_.index_i(c,2,site)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry]; 
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,1);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];   
  double U[2*NC_*NC_];
  is = 0;
  for(int k=0; k<Nbdry; ++k){  
    int site = bdry_pup_[1][k];    

    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),1)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),1)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(ff_.index_r(c,0,site), v1r[c]);
      fp.add(ff_.index_i(c,0,site), v1i[c]);
      fp.add(ff_.index_r(c,1,site), v2r[c]);
      fp.add(ff_.index_i(c,1,site), v2i[c]);
      fp.add(ff_.index_r(c,2,site),-v2r[c]);
      fp.add(ff_.index_i(c,2,site),-v2i[c]);
      fp.add(ff_.index_r(c,3,site), v1r[c]);
      fp.add(ff_.index_i(c,3,site), v1i[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  for(int k=0; k<Nbulk; ++k){

    int site = bulk_pup_[1][k];
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),1)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),1)];
      }
    }
    int yp = (this->*x_p)(site,1);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,yp)] +f[ff_.index_r(c,3,yp)];
      w1i[c] = f[ff_.index_i(c,0,yp)] +f[ff_.index_i(c,3,yp)];
      w2r[c] = f[ff_.index_r(c,1,yp)] -f[ff_.index_r(c,2,yp)];
      w2i[c] = f[ff_.index_i(c,1,yp)] -f[ff_.index_i(c,2,yp)];
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(ff_.index_r(c,0,site), v1r[c]);
      fp.add(ff_.index_i(c,0,site), v1i[c]);
      fp.add(ff_.index_r(c,1,site), v2r[c]);
      fp.add(ff_.index_i(c,1,site), v2i[c]);
      fp.add(ff_.index_r(c,2,site),-v2r[c]);
      fp.add(ff_.index_i(c,2,site),-v2i[c]);
      fp.add(ff_.index_r(c,3,site), v1r[c]);
      fp.add(ff_.index_i(c,3,site), v1i[c]);
    }
  }
}

void Dirac_Wilson::mult_zp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; 
  int Nbdry = bdry_plw_[2].size();
  int Nbulk = bulk_pup_[2].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_plw_[2][k];
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = f[ff_.index_r(c,0,site)] -f[ff_.index_i(c,2,site)];
      vbd[i0(c)+is] = f[ff_.index_i(c,0,site)] +f[ff_.index_r(c,2,site)];
      vbd[r1(c)+is] = f[ff_.index_r(c,1,site)] +f[ff_.index_i(c,3,site)];
      vbd[i1(c)+is] = f[ff_.index_i(c,1,site)] -f[ff_.index_r(c,3,site)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry]; 
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,2);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];   
  double U[2*NC_*NC_];
  is = 0;
  for(int k=0; k<Nbdry; ++k){  
    int site = bdry_pup_[2][k];    

    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),2)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),2)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(ff_.index_r(c,0,site), v1r[c]);
      fp.add(ff_.index_i(c,0,site), v1i[c]);
      fp.add(ff_.index_r(c,1,site), v2r[c]);
      fp.add(ff_.index_i(c,1,site), v2i[c]);
      fp.add(ff_.index_r(c,2,site), v1i[c]);
      fp.add(ff_.index_i(c,2,site),-v1r[c]);
      fp.add(ff_.index_r(c,3,site),-v2i[c]);
      fp.add(ff_.index_i(c,3,site), v2r[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  for(int k=0; k<Nbulk; ++k){

    int site = bulk_pup_[2][k];
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),2)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),2)];
      }
    }
    int zp = (this->*x_p)(site,2);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,zp)] -f[ff_.index_i(c,2,zp)];
      w1i[c] = f[ff_.index_i(c,0,zp)] +f[ff_.index_r(c,2,zp)];
      w2r[c] = f[ff_.index_r(c,1,zp)] +f[ff_.index_i(c,3,zp)];
      w2i[c] = f[ff_.index_i(c,1,zp)] -f[ff_.index_r(c,3,zp)];
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(ff_.index_r(c,0,site), v1r[c]);
      fp.add(ff_.index_i(c,0,site), v1i[c]);
      fp.add(ff_.index_r(c,1,site), v2r[c]);
      fp.add(ff_.index_i(c,1,site), v2i[c]);
      fp.add(ff_.index_r(c,2,site), v1i[c]);
      fp.add(ff_.index_i(c,2,site),-v1r[c]);
      fp.add(ff_.index_r(c,3,site),-v2i[c]);
      fp.add(ff_.index_i(c,3,site), v2r[c]);
    }
  }
}

void Dirac_Wilson::mult_tp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; 
  int Nbdry = bdry_plw_[3].size();
  int Nbulk = bulk_pup_[3].size();

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_plw_[3][k];
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = f[ff_.index_r(c,2,site)]*2.0;
      vbd[i0(c)+is] = f[ff_.index_i(c,2,site)]*2.0;
      vbd[r1(c)+is] = f[ff_.index_r(c,3,site)]*2.0;
      vbd[i1(c)+is] = f[ff_.index_i(c,3,site)]*2.0;
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry]; 
  Communicator::instance()->transfer_fw(vbc,vbd,Nih*Nbdry,2);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];   
  double U[2*NC_*NC_];
  is = 0;
  for(int k=0; k<Nbdry; ++k){  
    int site = bdry_pup_[3][k];    

    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),3)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),3)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(ff_.index_r(c,2,site), v1r[c]);
      fp.add(ff_.index_i(c,2,site), v1i[c]);
      fp.add(ff_.index_r(c,3,site), v2r[c]);
      fp.add(ff_.index_i(c,3,site), v2i[c]);
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  for(int k=0; k<Nbulk; ++k){

    int site = bulk_pup_[3][k];
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gp)(site),3)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gp)(site),3)];
      }
    }
    int tp = (this->*x_p)(site,3);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,2,tp)]*2.0;
      w1i[c] = f[ff_.index_i(c,2,tp)]*2.0;
      w2r[c] = f[ff_.index_r(c,3,tp)]*2.0;
      w2i[c] = f[ff_.index_i(c,3,tp)]*2.0;
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0;  v1i[c] = 0.0;
      v2r[c] = 0.0;  v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
      }
      fp.add(ff_.index_r(c,2,site), v1r[c]);
      fp.add(ff_.index_i(c,2,site), v1i[c]);
      fp.add(ff_.index_r(c,3,site), v2r[c]);
      fp.add(ff_.index_i(c,3,site), v2i[c]);
    }
  }
}

void Dirac_Wilson::mult_xm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[0].size();
  int Nbulk = bulk_mlw_[0].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  double U[2*NC_*NC_];

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mup_[0][k];
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,site)] +f[ff_.index_i(c,3,site)];  
      w1i[c] = f[ff_.index_i(c,0,site)] -f[ff_.index_r(c,3,site)];
      w2r[c] = f[ff_.index_r(c,1,site)] +f[ff_.index_i(c,2,site)]; 
      w2i[c] = f[ff_.index_i(c,1,site)] -f[ff_.index_r(c,2,site)];
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(site),0)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(site),0)];
      }
    }
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,0);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mlw_[0][k];
    for(int c=0; c<NC_; ++c){  
      fm.add(ff_.index_r(c,0,site), vbc[r0(c)+is]);  
      fm.add(ff_.index_i(c,0,site), vbc[i0(c)+is]);
      fm.add(ff_.index_r(c,1,site), vbc[r1(c)+is]);   
      fm.add(ff_.index_i(c,1,site), vbc[i1(c)+is]);
      fm.add(ff_.index_r(c,2,site),-vbc[i1(c)+is]);
      fm.add(ff_.index_i(c,2,site), vbc[r1(c)+is]);
      fm.add(ff_.index_r(c,3,site),-vbc[i0(c)+is]); 
      fm.add(ff_.index_i(c,3,site), vbc[r0(c)+is]);
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; 
  for(int k=0; k<Nbulk; ++k){
    int site = bulk_mlw_[0][k];
    int xm = (this->*x_m)(site,0);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,xm)] +f[ff_.index_i(c,3,xm)]; 
      w1i[c] = f[ff_.index_i(c,0,xm)] -f[ff_.index_r(c,3,xm)];
      w2r[c] = f[ff_.index_r(c,1,xm)] +f[ff_.index_i(c,2,xm)];
      w2i[c] = f[ff_.index_i(c,1,xm)] -f[ff_.index_r(c,2,xm)];
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(xm),0)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(xm),0)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(ff_.index_r(c,0,site), v1r[c]);  
      fm.add(ff_.index_i(c,0,site), v1i[c]);
      fm.add(ff_.index_r(c,1,site), v2r[c]); 
      fm.add(ff_.index_i(c,1,site), v2i[c]);
      fm.add(ff_.index_r(c,2,site),-v2i[c]);  
      fm.add(ff_.index_i(c,2,site), v2r[c]);
      fm.add(ff_.index_r(c,3,site),-v1i[c]); 
      fm.add(ff_.index_i(c,3,site), v1r[c]);
    }
  }
}

void Dirac_Wilson::mult_ym(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[1].size();
  int Nbulk = bulk_mlw_[1].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  double U[2*NC_*NC_];

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mup_[1][k];
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,site)] -f[ff_.index_r(c,3,site)];  
      w1i[c] = f[ff_.index_i(c,0,site)] -f[ff_.index_i(c,3,site)];
      w2r[c] = f[ff_.index_r(c,1,site)] +f[ff_.index_r(c,2,site)]; 
      w2i[c] = f[ff_.index_i(c,1,site)] +f[ff_.index_i(c,2,site)];
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(site),1)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(site),1)];
      }
    }
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,1);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mlw_[1][k];
    for(int c=0; c<NC_; ++c){  
      fm.add(ff_.index_r(c,0,site), vbc[r0(c)+is]);  
      fm.add(ff_.index_i(c,0,site), vbc[i0(c)+is]);
      fm.add(ff_.index_r(c,1,site), vbc[r1(c)+is]);   
      fm.add(ff_.index_i(c,1,site), vbc[i1(c)+is]);
      fm.add(ff_.index_r(c,2,site), vbc[r1(c)+is]);
      fm.add(ff_.index_i(c,2,site), vbc[i1(c)+is]);
      fm.add(ff_.index_r(c,3,site),-vbc[r0(c)+is]); 
      fm.add(ff_.index_i(c,3,site),-vbc[i0(c)+is]);
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; 
  for(int k=0; k<Nbulk; ++k){
    int site = bulk_mlw_[1][k];
    int ym = (this->*x_m)(site,1);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,ym)] -f[ff_.index_r(c,3,ym)]; 
      w1i[c] = f[ff_.index_i(c,0,ym)] -f[ff_.index_i(c,3,ym)];
      w2r[c] = f[ff_.index_r(c,1,ym)] +f[ff_.index_r(c,2,ym)];
      w2i[c] = f[ff_.index_i(c,1,ym)] +f[ff_.index_i(c,2,ym)];
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(ym),1)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(ym),1)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(ff_.index_r(c,0,site), v1r[c]);  
      fm.add(ff_.index_i(c,0,site), v1i[c]);
      fm.add(ff_.index_r(c,1,site), v2r[c]); 
      fm.add(ff_.index_i(c,1,site), v2i[c]);
      fm.add(ff_.index_r(c,2,site), v2r[c]);  
      fm.add(ff_.index_i(c,2,site), v2i[c]);
      fm.add(ff_.index_r(c,3,site),-v1r[c]); 
      fm.add(ff_.index_i(c,3,site),-v1i[c]);
    }
  }
}

void Dirac_Wilson::mult_zm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[2].size();
  int Nbulk = bulk_mlw_[2].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  double U[2*NC_*NC_];

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mup_[2][k];
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,site)] +f[ff_.index_i(c,2,site)];  
      w1i[c] = f[ff_.index_i(c,0,site)] -f[ff_.index_r(c,2,site)];
      w2r[c] = f[ff_.index_r(c,1,site)] -f[ff_.index_i(c,3,site)]; 
      w2i[c] = f[ff_.index_i(c,1,site)] +f[ff_.index_r(c,3,site)];
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(site),2)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(site),2)];
      }
    }
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,2);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mlw_[2][k];
    for(int c=0; c<NC_; ++c){  
      fm.add(ff_.index_r(c,0,site), vbc[r0(c)+is]);  
      fm.add(ff_.index_i(c,0,site), vbc[i0(c)+is]);
      fm.add(ff_.index_r(c,1,site), vbc[r1(c)+is]);   
      fm.add(ff_.index_i(c,1,site), vbc[i1(c)+is]);
      fm.add(ff_.index_r(c,2,site),-vbc[i0(c)+is]);
      fm.add(ff_.index_i(c,2,site), vbc[r0(c)+is]);
      fm.add(ff_.index_r(c,3,site), vbc[i1(c)+is]); 
      fm.add(ff_.index_i(c,3,site),-vbc[r1(c)+is]);
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; 
  for(int k=0; k<Nbulk; ++k){
    int site = bulk_mlw_[2][k];
    int zm = (this->*x_m)(site,2);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,zm)] +f[ff_.index_i(c,2,zm)]; 
      w1i[c] = f[ff_.index_i(c,0,zm)] -f[ff_.index_r(c,2,zm)];
      w2r[c] = f[ff_.index_r(c,1,zm)] -f[ff_.index_i(c,3,zm)];
      w2i[c] = f[ff_.index_i(c,1,zm)] +f[ff_.index_r(c,3,zm)];
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(zm),2)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(zm),2)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(ff_.index_r(c,0,site), v1r[c]);  
      fm.add(ff_.index_i(c,0,site), v1i[c]);
      fm.add(ff_.index_r(c,1,site), v2r[c]); 
      fm.add(ff_.index_i(c,1,site), v2i[c]);
      fm.add(ff_.index_r(c,2,site),-v1i[c]);  
      fm.add(ff_.index_i(c,2,site), v1r[c]);
      fm.add(ff_.index_r(c,3,site), v2i[c]); 
      fm.add(ff_.index_i(c,3,site),-v2r[c]);
    }
  }
}

void Dirac_Wilson::mult_tm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  int Nbdry = bdry_mup_[3].size();
  int Nbulk = bulk_mlw_[3].size();

  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];
  double U[2*NC_*NC_];

  /// boundary part ///
  double vbd[Nih*Nbdry]; 
  int is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mup_[3][k];
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,site)]*2.0;
      w1i[c] = f[ff_.index_i(c,0,site)]*2.0;
      w2r[c] = f[ff_.index_r(c,1,site)]*2.0;
      w2i[c] = f[ff_.index_i(c,1,site)]*2.0;
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(site),3)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(site),3)];
      }
    }
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	vbd[i0(c)+is] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	vbd[r1(c)+is] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	vbd[i1(c)+is] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Nbdry,3);
  is = 0;
  for(int k=0; k<Nbdry; ++k){
    int site = bdry_mlw_[3][k];
    for(int c=0; c<NC_; ++c){  
      fm.add(ff_.index_r(c,0,site), vbc[r0(c)+is]);  
      fm.add(ff_.index_i(c,0,site), vbc[i0(c)+is]);
      fm.add(ff_.index_r(c,1,site), vbc[r1(c)+is]);   
      fm.add(ff_.index_i(c,1,site), vbc[i1(c)+is]);
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; 
  for(int k=0; k<Nbulk; ++k){
    int site = bulk_mlw_[3][k];
    int tm = (this->*x_m)(site,3);
    for(int c=0; c<NC_; ++c){
      w1r[c] = f[ff_.index_r(c,0,tm)]*2.0;
      w1i[c] = f[ff_.index_i(c,0,tm)]*2.0;
      w2r[c] = f[ff_.index_r(c,1,tm)]*2.0;
      w2i[c] = f[ff_.index_i(c,1,tm)]*2.0;
    }
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	U[re(c1,c2)] = (*u_)[gf_.index_r(c1,c2,(this->*gm)(tm),3)];
	U[im(c1,c2)] = (*u_)[gf_.index_i(c1,c2,(this->*gm)(tm),3)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
      }
      fm.add(ff_.index_r(c,0,site), v1r[c]);  
      fm.add(ff_.index_i(c,0,site), v1i[c]);
      fm.add(ff_.index_r(c,1,site), v2r[c]); 
      fm.add(ff_.index_i(c,1,site), v2i[c]);
    }
  }
}

const Field Dirac_Wilson::gamma5(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site), f[ff_.index_r(c,2,site)]);
      w.set(ff_.index_i(c,0,site), f[ff_.index_i(c,2,site)]);
      w.set(ff_.index_r(c,1,site), f[ff_.index_r(c,3,site)]);
      w.set(ff_.index_i(c,1,site), f[ff_.index_i(c,3,site)]);
      w.set(ff_.index_r(c,2,site), f[ff_.index_r(c,0,site)]);
      w.set(ff_.index_i(c,2,site), f[ff_.index_i(c,0,site)]);
      w.set(ff_.index_r(c,3,site), f[ff_.index_r(c,1,site)]);
      w.set(ff_.index_i(c,3,site), f[ff_.index_i(c,1,site)]);
    }
  }
  return w;
}

const Field Dirac_Wilson::proj_p(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      double fup_r = 0.5*(f[ff_.index_r(c,0,site)]+f[ff_.index_r(c,2,site)]);
      double fup_i = 0.5*(f[ff_.index_i(c,0,site)]+f[ff_.index_i(c,2,site)]);
      double fdn_r = 0.5*(f[ff_.index_r(c,1,site)]+f[ff_.index_r(c,3,site)]);
      double fdn_i = 0.5*(f[ff_.index_i(c,1,site)]+f[ff_.index_i(c,3,site)]);

      w.set(ff_.index_r(c,0,site),fup_r);
      w.set(ff_.index_i(c,0,site),fup_i);
      w.set(ff_.index_r(c,1,site),fdn_r);
      w.set(ff_.index_i(c,1,site),fdn_i);
      w.set(ff_.index_r(c,2,site),fup_r);
      w.set(ff_.index_i(c,2,site),fup_i);
      w.set(ff_.index_r(c,3,site),fdn_r);
      w.set(ff_.index_i(c,3,site),fdn_i);
    }
  }
  return w;
}

const Field Dirac_Wilson::proj_m(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      double fup_r = 0.5*(f[ff_.index_r(c,0,site)]-f[ff_.index_r(c,2,site)]);
      double fup_i = 0.5*(f[ff_.index_i(c,0,site)]-f[ff_.index_i(c,2,site)]);
      double fdn_r = 0.5*(f[ff_.index_r(c,1,site)]-f[ff_.index_r(c,3,site)]);
      double fdn_i = 0.5*(f[ff_.index_i(c,1,site)]-f[ff_.index_i(c,3,site)]);

      w.set(ff_.index_r(c,0,site),fup_r);
      w.set(ff_.index_i(c,0,site),fup_i);
      w.set(ff_.index_r(c,1,site),fdn_r);
      w.set(ff_.index_i(c,1,site),fdn_i);
      w.set(ff_.index_r(c,2,site),-fup_r);
      w.set(ff_.index_i(c,2,site),-fup_i);
      w.set(ff_.index_r(c,3,site),-fdn_r);
      w.set(ff_.index_i(c,3,site),-fdn_i);
    }
  }
  return w;
}

#endif  /*no IMPROVED_WILSON*/
/////////////////////////////////////////////////////////////////////////

void (Dirac_Wilson::*Dirac_Wilson::mult_p[])
(Field&,const Field&) const = {&Dirac_Wilson::mult_xp,
			       &Dirac_Wilson::mult_yp,
			       &Dirac_Wilson::mult_zp,
			       &Dirac_Wilson::mult_tp,};

void (Dirac_Wilson::*Dirac_Wilson::mult_m[])
(Field&,const Field&) const = {&Dirac_Wilson::mult_xm,
			       &Dirac_Wilson::mult_ym,
			       &Dirac_Wilson::mult_zm,
			       &Dirac_Wilson::mult_tm,};

void Dirac_Wilson::mult_offdiag(Field& w, const Field& f) const{
  for(int d=0; d <Ndim_; ++d){
    //  for(int d=0; d <1; ++d){
    (this->*mult_p[d])(w,f);
    (this->*mult_m[d])(w,f);
  }
  w *= -kpp_;
}
void Dirac_Wilson::mult_full(Field& w, const Field& f) const{
  mult_offdiag(w,f);
  w += f;
}

const Field Dirac_Wilson::mult(const Field& f) const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
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
  SUNmat f;

  for(int mu=0; mu<Ndim_; ++mu){
    Field xie(fsize_);

    (this->*mult_p[mu])(xie, eta);

    for(int site=0; site<Nvol_; ++site){
      f = 0.0;
      for(int a=0; a<NC_; ++a){
        for(int b=0; b<NC_; ++b){
          double fre = 0.0;
          double fim = 0.0;
          for(int s=0; s<Nd_; ++s){

	    size_t ra =ff_.index_r(a,s,site);
	    size_t ia =ff_.index_i(a,s,site);

	    size_t rb =ff_.index_r(b,s,site);
	    size_t ib =ff_.index_i(b,s,site);

	    fre += zeta[rb]*xie[ra] +zeta[ib]*xie[ia];
	    fim += zeta[rb]*xie[ia] -zeta[ib]*xie[ra];
          }
          f.set(a,b,fre,fim);
        }
      }
      int gsite = (this->*gp)(site);
      fce.add(gf_.cslice(0,gsite,mu),f.getva());
    }
  }
}

void Dirac_Wilson::md_force_m(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmat_utils;
  SUNmat f;
  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

  for(int mu=0; mu<Ndim_; ++mu){
    Field xz5(fsize_);
    (this->*mult_p[mu])(xz5, zt5);

    for(int site=0; site<Nvol_; ++site){
      f=0.0;
      for(int a=0; a<NC_; ++a){
        for(int b=0; b<NC_; ++b){
          double fre = 0.0;
          double fim = 0.0;
          for(int s=0; s<Nd_; ++s){

	    size_t ra =ff_.index_r(a,s,site);
	    size_t ia =ff_.index_i(a,s,site);

	    size_t rb =ff_.index_r(b,s,site);
	    size_t ib =ff_.index_i(b,s,site);

	    fre -= xz5[rb]*et5[ra] +xz5[ib]*et5[ia];
	    fim -= xz5[rb]*et5[ia] -xz5[ib]*et5[ra];
          }
          f.set(a,b,fre,fim);
        }
      }
      int gsite = (this->*gp)(site);
      fce.add(gf_.cslice(0,gsite,mu),f.getva());
    }
  }
}

const Field Dirac_Wilson::
md_force(const Field& eta,const Field& zeta)const{
  
  Field fp(gf_.size());
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
