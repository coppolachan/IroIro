/* Improved version of Dirac Kernel 
   Time-stamp: <2013-05-03 09:19:16 noaki>
*/

void Dirac_Wilson::mult_xp(Field& fp, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num of elements of a half spinor */

  /// boundary part ///
  int is = 0;
  int Xb = 0;
  int Nbdry = (this->*slice_isize)(Xb,XDIR);
  double vbd[Nih*Nbdry]; /*!< @brief data on the lower slice */   
  for(int k=0; k<Nbdry; ++k){
    const double* v 
      = const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(Xb,k,XDIR)));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r0(c)]-v[i3(c)]; vbd[i0(c)+is] = v[i0(c)]+v[r3(c)];
      vbd[r1(c)+is] = v[r1(c)]-v[i2(c)]; vbd[i1(c)+is] = v[i1(c)]+v[r2(c)];
    }
    is += Nih;
  }

  double vbc[Nih*Nbdry]; /*!< @brief data on the upper slice */   
  comm_->transfer_fw(vbc,vbd,Nih*Nbdry,XDIR);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  Xb = Nx_-1;
  Nbdry = (this->*slice_osize)(Xb,XDIR);

  for(int k=0; k<Nbdry; ++k){  /*!< @brief calc on the upper boundary */   
    int xc = (this->*slice_out)(Xb,k,XDIR);
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(xc),XDIR));
    double* res = fp.getaddr(ff_.index(0,xc));    

    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
      res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
      res[r2(c)] += v2i[c];      res[i2(c)] -= v2r[c];
      res[r3(c)] += v1i[c];      res[i3(c)] -= v1r[c];
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int x=0; x<Xb; ++x){
    int Nslice = (this->*slice_osize)(x,XDIR);
    for(int k=0; k<Nslice; ++k){   
      const double* v 
      	= const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(x+1,k,XDIR)));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)] -v[i3(c)];  w1i[c] = v[i0(c)] +v[r3(c)];
	w2r[c] = v[r1(c)] -v[i2(c)];  w2i[c] = v[i1(c)] +v[r2(c)];
      }
      int xc = (this->*slice_out)(x,k,XDIR);
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(xc),XDIR));
      double* res = fp.getaddr(ff_.index(0,xc));

      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	  v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	  v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	  v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
	res[r2(c)] += v2i[c];      res[i2(c)] -= v2r[c];
	res[r3(c)] += v1i[c];      res[i3(c)] -= v1r[c];
      }
    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num ob elements of a half spinor */
  /// boundary part ///
  int is = 0;
  int Yb = 0;
  int Nbdry = (this->*slice_isize)(Yb,YDIR);
  double vbd[Nih*Nbdry];
  for(int k=0; k<Nbdry; ++k){
    const double* v 
      = const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(Yb,k,YDIR)));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r0(c)]+v[r3(c)]; vbd[i0(c)+is] = v[i0(c)]+v[i3(c)];
      vbd[r1(c)+is] = v[r1(c)]-v[r2(c)]; vbd[i1(c)+is] = v[i1(c)]-v[i2(c)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry]; 
  comm_->transfer_fw(vbc,vbd,Nih*Nbdry,YDIR);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  Yb = Ny_-1;
  Nbdry = (this->*slice_osize)(Yb,YDIR);

  for(int k=0; k<Nbdry; ++k){

    int yc = (this->*slice_out)(Yb,k,YDIR);
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(yc),YDIR));
    double* res = fp.getaddr(ff_.index(0,yc));
    
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
      res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
      res[r2(c)] -= v2r[c];      res[i2(c)] -= v2i[c];
      res[r3(c)] += v1r[c];      res[i3(c)] += v1i[c];
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int y=0; y<Yb; ++y){
    int Nslice = (this->*slice_osize)(y,YDIR);
    for(int k=0; k<Nslice; ++k){
      const double* v 
	= const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(y+1,k,YDIR)));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)] +v[r3(c)];  w1i[c] = v[i0(c)] +v[i3(c)];
	w2r[c] = v[r1(c)] -v[r2(c)];  w2i[c] = v[i1(c)] -v[i2(c)];
      }
      int yc = (this->*slice_out)(y,k,YDIR);
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(yc),YDIR));
      double* res = fp.getaddr(ff_.index(0,yc));
      
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	  v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	  v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	  v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
	res[r2(c)] -= v2r[c];      res[i2(c)] -= v2i[c];
	res[r3(c)] += v1r[c];      res[i3(c)] += v1i[c];
      }
    }
  }
}

void Dirac_Wilson::mult_zp(Field& fp, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num ob elements of a half spinor */

  /// boundary part ///
  int is = 0;
  int Zb = 0;
  int Nbdry = (this->*slice_isize)(Zb,ZDIR);
  double vbd[Nih*Nbdry]; 
  for(int k=0; k<Nbdry; ++k){
    const double* v 
      = const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(Zb,k,ZDIR)));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r0(c)]-v[i2(c)]; vbd[i0(c)+is] = v[i0(c)]+v[r2(c)];
      vbd[r1(c)+is] = v[r1(c)]+v[i3(c)]; vbd[i1(c)+is] = v[i1(c)]-v[r3(c)];
    }
    is += Nih;
  }
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  comm_->transfer_fw(vbc,vbd,Nih*Nbdry,ZDIR);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  Zb = Nz_-1;
  Nbdry = (this->*slice_osize)(Zb,ZDIR);
  for(int k=0; k<Nbdry; ++k){

    int zc = (this->*slice_out)(Zb,k,ZDIR);
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(zc),ZDIR));
    double* res = fp.getaddr(ff_.index(0,zc));
      
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
      res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
      res[r2(c)] += v1i[c];      res[i2(c)] -= v1r[c];
      res[r3(c)] -= v2i[c];      res[i3(c)] += v2r[c];
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int z=0; z<Zb; ++z){
    int Nslice = (this->*slice_osize)(z,ZDIR);
    for(int k=0; k<Nslice; ++k){
      const double* v 
	= const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(z+1,k,ZDIR)));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)] -v[i2(c)];  w1i[c] = v[i0(c)] +v[r2(c)];
	w2r[c] = v[r1(c)] +v[i3(c)];  w2i[c] = v[i1(c)] -v[r3(c)];
      }
      int zc = (this->*slice_out)(z,k,ZDIR);
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(zc),ZDIR));
      double* res = fp.getaddr(ff_.index(0,zc));
      
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	  v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	  v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	  v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
	res[r2(c)] += v1i[c];      res[i2(c)] -= v1r[c];
	res[r3(c)] -= v2i[c];      res[i3(c)] += v2r[c];
      }
    }
  }
}

void Dirac_Wilson::mult_tp(Field& fp, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num ob elements of a half spinor */

  /// boundary part ///
  int is = 0;
  int Tb = 0;
  int Nbdry = (this->*slice_isize)(Tb,TDIR);
  double vbd[Nih*Nbdry]; 
  for(int k=0; k<Nbdry; ++k){
    const double* v 
      = const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(Tb,k,TDIR)));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = v[r2(c)]*2.0;  vbd[i0(c)+is] = v[i2(c)]*2.0;
      vbd[r1(c)+is] = v[r3(c)]*2.0;  vbd[i1(c)+is] = v[i3(c)]*2.0;
    }
    is += Nih;
  }  
  double vbc[Nih*Nbdry];  //Copy vbd from backward processor
  comm_->transfer_fw(vbc,vbd,Nih*Nbdry,TDIR);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     
  is = 0;
  Tb = Nt_-1;
  Nbdry = (this->*slice_osize)(Tb,TDIR);
  for(int k=0; k<Nbdry; ++k){

    int tc = (this->*slice_out)(Tb,k,TDIR);
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(tc),TDIR));
    double* res = fp.getaddr(ff_.index(0,tc));
    
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= U[re(c,c1)]*vbc[r0(c1)+is] -U[im(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= U[im(c,c1)]*vbc[r0(c1)+is] +U[re(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= U[re(c,c1)]*vbc[r1(c1)+is] -U[im(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= U[im(c,c1)]*vbc[r1(c1)+is] +U[re(c,c1)]*vbc[i1(c1)+is];
      }
      res[r2(c)] += v1r[c];      res[i2(c)] += v1i[c];
      res[r3(c)] += v2r[c];      res[i3(c)] += v2i[c];
    }
    is += Nih;
  }
  /// bulk part ///
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  for(int t=0; t<Tb; ++t){
    int Nslice = (this->*slice_osize)(t,TDIR);
    for(int k=0; k<Nslice; ++k){
      const double* v 
	= const_cast<Field&>(f).getaddr(ff_.index(0,(this->*slice_in)(t+1,k,TDIR)));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r2(c)]*2.0;  w1i[c] = v[i2(c)]*2.0;
	w2r[c] = v[r3(c)]*2.0;  w2i[c] = v[i3(c)]*2.0;
      }

      int tc = (this->*slice_out)(t,k,TDIR);
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gp)(tc),TDIR));
      double* res = fp.getaddr(ff_.index(0,tc));

      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c,c1)]*w1r[c1] -U[im(c,c1)]*w1i[c1];
	  v1i[c] += U[im(c,c1)]*w1r[c1] +U[re(c,c1)]*w1i[c1];
	  v2r[c] += U[re(c,c1)]*w2r[c1] -U[im(c,c1)]*w2i[c1];
	  v2i[c] += U[im(c,c1)]*w2r[c1] +U[re(c,c1)]*w2i[c1];
	}
	res[r2(c)] += v1r[c];      res[i2(c)] += v1i[c];
	res[r3(c)] += v2r[c];      res[i3(c)] += v2i[c];
      }
    }
  }
}

void Dirac_Wilson::mult_xm(Field& fm, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num ob elements of a half spinor */
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  /// boundary part ///
  int is = 0;
  int Xb = Nx_-1;
  int Nbdry = (this->*slice_isize)(Xb,XDIR);
  double vbd[Nih*Nbdry]; /*!< @brief data on the upper boundary */
  for(int k=0; k<Nbdry; ++k){
    int xc = (this->*slice_in)(Xb,k,XDIR);
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xc));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[i3(c)];  w1i[c] = v[i0(c)] -v[r3(c)];
      w2r[c] = v[r1(c)] +v[i2(c)];  w2i[c] = v[i1(c)] -v[r2(c)];
    }
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(xc),XDIR));
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
  comm_->transfer_bk(vbc,vbd,Nih*Nbdry,XDIR);
  is = 0;
  Xb = 0;
  Nbdry = (this->*slice_osize)(Xb,XDIR);
  
  for(int k=0; k<Nbdry; ++k){
    double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(Xb,k,XDIR)));

    for(int c=0; c<NC_; ++c){  
      res[r0(c)] += vbc[r0(c)+is];      res[i0(c)] += vbc[i0(c)+is];
      res[r1(c)] += vbc[r1(c)+is];      res[i1(c)] += vbc[i1(c)+is];
      res[r2(c)] -= vbc[i1(c)+is];      res[i2(c)] += vbc[r1(c)+is];
      res[r3(c)] -= vbc[i0(c)+is];      res[i3(c)] += vbc[r0(c)+is];
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; 

  for(int x=1; x<Nx_; ++x){
    int Nslice = (this->*slice_isize)(x,XDIR);
    for(int k=0; k<Nslice; ++k){
      int xm = (this->*slice_in)(x-1,k,XDIR);
      const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xm));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)] +v[i3(c)];  w1i[c] = v[i0(c)] -v[r3(c)];
	w2r[c] = v[r1(c)] +v[i2(c)];  w2i[c] = v[i1(c)] -v[r2(c)];
      }
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(xm),XDIR));
      double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(x,k,XDIR)));

      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0; v1i[c] = 0.0; 
	v2r[c] = 0.0; v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	  v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	  v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	  v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
	res[r2(c)] -= v2i[c];      res[i2(c)] += v2r[c];
	res[r3(c)] -= v1i[c];      res[i3(c)] += v1r[c];
      }
    }
  }
}

void Dirac_Wilson::mult_ym(Field& fm, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num of elements of a half spinor */
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  // boundary part  
  int is = 0;
  int Yb = Ny_-1;
  int Nbdry = (this->*slice_isize)(Yb,YDIR);
  double vbd[Nih*Nbdry];
  for(int k=0; k<Nbdry; ++k){
    int yc = (this->*slice_in)(Yb,k,YDIR);
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,yc));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] -v[r3(c)];  w1i[c] = v[i0(c)] -v[i3(c)];
      w2r[c] = v[r1(c)] +v[r2(c)];  w2i[c] = v[i1(c)] +v[i2(c)];
    }
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(yc),YDIR));
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
  comm_->transfer_bk(vbc,vbd,Nih*Nbdry,YDIR);
  is = 0;
  Yb = 0;
  Nbdry = (this->*slice_osize)(Yb,YDIR);

  for(int k=0; k<Nbdry; ++k){
    double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(Yb,k,YDIR)));
    for(int c=0; c<NC_; ++c){  
      res[r0(c)] += vbc[r0(c)+is];      res[i0(c)] += vbc[i0(c)+is];
      res[r1(c)] += vbc[r1(c)+is];      res[i1(c)] += vbc[i1(c)+is];
      res[r2(c)] += vbc[r1(c)+is];      res[i2(c)] += vbc[i1(c)+is];
      res[r3(c)] -= vbc[r0(c)+is];      res[i3(c)] -= vbc[i0(c)+is];
    }
    is += Nih;
  }
  //bulk part
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int y=1; y<Ny_; ++y){
    int Nslice = (this->*slice_isize)(y,YDIR);
    for(int k=0; k<Nslice; ++k){
      int ym = (this->*slice_in)(y-1,k,YDIR);
      const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,ym));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)] -v[r3(c)];  w1i[c] = v[i0(c)] -v[i3(c)];
	w2r[c] = v[r1(c)] +v[r2(c)];  w2i[c] = v[i1(c)] +v[i2(c)];
      }
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(ym),YDIR));
      double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(y,k,YDIR)));

      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0; v1i[c] = 0.0; 
	v2r[c] = 0.0; v2i[c] = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	  v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	  v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	  v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
	res[r2(c)] += v2r[c];      res[i2(c)] += v2i[c];
	res[r3(c)] -= v1r[c];      res[i3(c)] -= v1i[c];
      }
    }
  }
}
  
void Dirac_Wilson::mult_zm(Field& fm, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num ob elements of a half spinor */
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  // boundary part
  int is = 0;
  int Zb = Nz_-1;
  int Nbdry = (this->*slice_isize)(Zb,ZDIR);
  double vbd[Nih*Nbdry];
  for(int k=0; k<Nbdry; ++k){
    int zc = (this->*slice_in)(Zb,k,ZDIR);
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,zc));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)] +v[i2(c)];  w1i[c] = v[i0(c)] -v[r2(c)];
      w2r[c] = v[r1(c)] -v[i3(c)];  w2i[c] = v[i1(c)] +v[r3(c)];
    }
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(zc),ZDIR));
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
  double vbc[Nih*Nbdry];   //Copy vbd from backward processor
  comm_->transfer_bk(vbc,vbd,Nih*Nbdry,ZDIR);
  is = 0;
  Zb = 0;
  Nbdry = (this->*slice_osize)(Zb,ZDIR);

  for(int k=0; k<Nbdry; ++k){
    double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(Zb,k,ZDIR)));
    
    for(int c=0; c<NC_; ++c){  
      res[r0(c)] += vbc[r0(c)+is];      res[i0(c)] += vbc[i0(c)+is];
      res[r1(c)] += vbc[r1(c)+is];      res[i1(c)] += vbc[i1(c)+is];
      res[r2(c)] -= vbc[i0(c)+is];      res[i2(c)] += vbc[r0(c)+is];
      res[r3(c)] += vbc[i1(c)+is];      res[i3(c)] -= vbc[r1(c)+is];
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int z=1; z<Nz_; ++z){
    int Nslice = (this->*slice_isize)(z,ZDIR);
    for(int k=0; k<Nslice; ++k){
      int zm = (this->*slice_in)(z-1,k,ZDIR);
      const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,zm));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)] +v[i2(c)];  w1i[c] = v[i0(c)] -v[r2(c)];
	w2r[c] = v[r1(c)] -v[i3(c)];  w2i[c] = v[i1(c)] +v[r3(c)];
      }
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(zm),ZDIR));
      double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(z,k,ZDIR)));
      
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0; 
	v2r[c] = 0.0;  v2i[c] = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	  v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	  v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	  v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
	res[r2(c)] -= v1i[c];      res[i2(c)] += v1r[c];
	res[r3(c)] += v2i[c];      res[i3(c)] -= v2r[c];
      }
    }
  }
}

void Dirac_Wilson::mult_tm(Field& fm, const Field& f) const{
  int Nih = ND_*NC_; /*!< @brief num ob elements of a half spinor */
  
  double w1r[NC_], w1i[NC_], w2r[NC_], w2i[NC_];

  /// boundary part ///
  int is = 0;
  int Tb = Nt_-1;
  int Nbdry = (this->*slice_isize)(Tb,TDIR);
  double vbd[Nih*Nbdry];
  for(int k=0; k<Nbdry; ++k){
    int tc = (this->*slice_in)(Tb,k,TDIR);
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,tc));
    for(int c=0; c<NC_; ++c){
      w1r[c] = v[r0(c)]*2.0;   w1i[c] = v[i0(c)]*2.0;
      w2r[c] = v[r1(c)]*2.0;   w2i[c] = v[i1(c)]*2.0;
    }
    const double* U 
      = const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(tc),TDIR));
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
  comm_->transfer_bk(vbc,vbd,Nih*Nbdry,TDIR);
  is = 0;
  Tb = 0;
  Nbdry = (this->*slice_osize)(Tb,TDIR);

  for(int k=0; k<Nbdry; ++k){
    double* res = fm.getaddr(ff_.index(0,(this->*slice_out)(Tb,k,TDIR)));
    for(int c=0; c<NC_; ++c){  
      res[r0(c)] += vbc[r0(c)+is];      res[i0(c)] += vbc[i0(c)+is];
      res[r1(c)] += vbc[r1(c)+is];      res[i1(c)] += vbc[i1(c)+is];
    }
    is += Nih;
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];
  
  for(int t=1; t<Nt_; ++t){
    int Nslice = (this->*slice_isize)(t,TDIR);
    for(int k=0; k<Nslice; ++k){
      int tm = (this->*slice_in)(t-1,k,TDIR);
      const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,tm));
      for(int c=0; c<NC_; ++c){
	w1r[c] = v[r0(c)]*2.0;   w1i[c] = v[i0(c)]*2.0;
	w2r[c] = v[r1(c)]*2.0;   w2i[c] = v[i1(c)]*2.0;
      }
      const double* U 
	= const_cast<Field*>(u_)->getaddr(gf_.index(0,(this->*gm)(tm),TDIR));
      double* res = fm.getaddr( ff_.index(0,(this->*slice_out)(t,k,TDIR)));
      
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0; 
	v2r[c] = 0.0;  v2i[c] = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += U[re(c1,c)]*w1r[c1] +U[im(c1,c)]*w1i[c1];
	  v1i[c] -= U[im(c1,c)]*w1r[c1] -U[re(c1,c)]*w1i[c1];
	  v2r[c] += U[re(c1,c)]*w2r[c1] +U[im(c1,c)]*w2i[c1];
	  v2i[c] -= U[im(c1,c)]*w2r[c1] -U[re(c1,c)]*w2i[c1];
	}
	res[r0(c)] += v1r[c];      res[i0(c)] += v1i[c];
	res[r1(c)] += v2r[c];      res[i1(c)] += v2i[c];
      }
    }
  }
}

