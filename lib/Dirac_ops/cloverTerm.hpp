class CloverTerm{
private:
  double csw_;
  int Nvol_;

  const Field* const u_;
  const ffmt_t ff_;
  const gfmt_t gf_;
  GaugeField1D Bx_, By_, Bz_, Ex_, Ey_, Ez_;

  void isigma_12(FermionField&, const FermionField&)const;
  void isigma_13(FermionField&, const FermionField&)const;
  void isigma_14(FermionField&, const FermionField&)const;

  void isigma_21(FermionField&, const FermionField&)const;
  void isigma_23(FermionField&, const FermionField&)const;
  void isigma_24(FermionField&, const FermionField&)const;

  void isigma_31(FermionField&, const FermionField&)const;
  void isigma_32(FermionField&, const FermionField&)const;
  void isigma_34(FermionField&, const FermionField&)const;

  void isigma_41(FermionField&, const FermionField&)const;
  void isigma_42(FermionField&, const FermionField&)const;
  void isigma_43(FermionField&, const FermionField&)const;

  void set_FieldStrength(GaugeField1D&, const int, const int);
public:
  CloverTerm(const Field* u, double csw)
    :u_(u),
     csw_(csw),
     Nvol_(CommonPrms::instance()->Nvol()),
     ff_(Nvol_),gf_(Nvol_){
    set_fieldstrength(Bx_, 1, 2);
    set_fieldstrength(By_, 2, 0);
    set_fieldstrength(Bz_, 0, 1);
    set_fieldstrength(Ex_, 3, 0);
    set_fieldstrength(Ey_, 3, 1);
    set_fieldstrength(Ez_, 3, 2);
  }
  const Field mult(const Field& ){  /* Dirac_Clover::mult_sw */ }
};


/////////////////////////////////////////////////
/*
CloverTerm* clover_;

Dirac_Wilson_Brillouin(double mass,const Field* u,
                       ImpType imp=Standard,csw=0.0)
		       : ,,clover_(null){
 if(csw) clover_= new CloverTerm(u_,csw);
}

Dirac_Wilson_Brillouin(const XML::node& node,const Field* u,
		       ImpType imp=Standard){
  double csw;
  XML::read(node,"csw",csw);
  if(csw) clover_= new CloverTerm(u_,csw);
}

~Dirac_Wilson_Brillouin(){ if(clover_) delete clover_; }

const Field mult(const Field& f){
  Field w(fsize_);
  ..

  if(clover_) {w +=clover_->mult(f) }
  return w;
}

*/
