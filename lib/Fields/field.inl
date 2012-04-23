inline Field& Field::operator-(){
  field= -field;
  return *this;
}
inline Field& Field::operator=(const Field& rhs){
  field= rhs.field;
  return *this;
}
inline Field& Field::operator=(const double r){
  field= r;
  return *this;
}
inline Field& Field::operator=(const std::valarray<double>& rhs){
  field= rhs;
  return *this;
}

inline Field& Field::operator+=(const Field& rhs){
  field+= rhs.field;
  return *this;
}
inline Field& Field::operator+=(const std::valarray<double>& rhs){
  field+= rhs;
  return *this;
}

inline Field& Field::operator-=(const Field& rhs){
  field-= rhs.field;
  return *this;
}
inline Field& Field::operator-=(const std::valarray<double>& rhs){
  field-= rhs;
  return *this;
}

inline Field& Field::operator*=(const double rhs){
  field*= rhs;
  return *this;
}

inline Field& Field::operator/=(const double rhs){
  field/= rhs;
  return *this;
}

inline double Field::operator*(const Field& rhs) const{
  double a = (field*rhs.field).sum();
  double b = Communicator::instance()->reduce_sum(a);
  return b;
}
inline double Field::operator*(const std::valarray<double>& rhs) const{
  double a = (field*rhs).sum();
  double b = Communicator::instance()->reduce_sum(a);
  return b;
}

inline double Field::im_prod(const Field& rhs) const{
  std::slice re(0,field.size()/2,2);
  std::slice im(1,field.size()/2,2);
  
  std::valarray<double> lhs_im = field[im];
  std::valarray<double> lhs_re = field[re];

  double a = (lhs_re*rhs[im]).sum()-(lhs_im*rhs[re]).sum();
  double b = Communicator::instance()->reduce_sum(a);
  return b;
}
inline double Field::im_prod(const std::valarray<double>& rhs) const{
  std::slice re(0,field.size()/2,2);
  std::slice im(1,field.size()/2,2);

  std::valarray<double> lhs_im = field[im];
  std::valarray<double> lhs_re = field[re];

  double a = (lhs_re*rhs[im]).sum()-(lhs_im*rhs[re]).sum();
  double b =  Communicator::instance()->reduce_sum(a);
  return b;
}
