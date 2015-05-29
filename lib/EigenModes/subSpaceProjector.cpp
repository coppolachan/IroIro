/*!@file subSpaceProjector.cpp
 * @brief implementation of the subSpaceProjector methods
 */
#include "subSpaceProjector.hpp"

using namespace std;

namespace SubSpace{
   
  void project_real(Field& w,const Field& f,
		    const vector_Field& vsub,const vector<double>& q){

    size_t vsize = f.size();   assert(vsize%2 ==0);
    slice re(0,vsize/2,2);
    slice im(1,vsize/2,2);

    int Nq = vsub.size();

    vector<double> sr(Nq);
    vector<double> si(Nq);

    for(int i=0; i<Nq; ++i){
      sr[i] = vsub[i]*f;
      si[i] = vsub[i].im_prod(f);
    }
    w = 0.0;
    for(int i=0; i<Nq; ++i){
      w.add(re,q[i]*(sr[i]*vsub[i][re] -si[i]*vsub[i][im]));
      w.add(im,q[i]*(sr[i]*vsub[i][im] +si[i]*vsub[i][re]));
    }
  }

  void project_complex(Field& w,const Field& f,
		       const vector_Field& vsub,const vector<double>& q){

    size_t vsize = f.size();   assert(vsize%2 ==0);
    slice re(0,vsize/2,2);
    slice im(1,vsize/2,2);

    int Nq = vsub.size();;
    
    vector<double> sr(Nq);
    vector<double> si(Nq);

    for(int i=0; i<Nq; ++i){
      double cr = vsub[i]*f;
      double ci = vsub[i].im_prod(f);
      sr[i] = cr*q[2*i  ] -ci*q[2*i+1];
      si[i] = cr*q[2*i+1] +ci*q[2*i];
    }
    w = 0.0;
    for(int i=0; i<Nq; ++i){
      w.add(re,sr[i]*vsub[i][re] -si[i]*vsub[i][im]);
      w.add(im,sr[i]*vsub[i][im] +si[i]*vsub[i][re]);
    }
  }

  void project(Field& w,const Field& f,const vector_Field& vsub){

    size_t vsize = f.size();   assert(vsize%2 ==0);
    slice re(0,vsize/2,2);
    slice im(1,vsize/2,2);

    int Nq = vsub.size();

    vector<double> sr(Nq);
    vector<double> si(Nq);

    for(int i=0; i<Nq; ++i){
      sr[i] = vsub[i]*f;
      si[i] = vsub[i].im_prod(f);
    }
    w = 0.0;
    for(int i=0; i<Nq; ++i){
      w.add(re,sr[i]*vsub[i][re] -si[i]*vsub[i][im]);
      w.add(im,sr[i]*vsub[i][im] +si[i]*vsub[i][re]);
    }
  }

  void projectOut(Field& w,const Field& f,const vector_Field& vsub){
    project(w,f,vsub);
    w -= f;
    w *= -1.0;
  } 
}
