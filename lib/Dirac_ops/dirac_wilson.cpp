//----------------------------------------------------------------------
// dirac_wilson.cpp
//----------------------------------------------------------------------
#include "dirac_wilson.h"
#include "Tools/sunMatUtils.hpp"

using namespace SUNvec_utils;
using namespace std;

void Dirac_Wilson::mult_xp(Field& fp, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;                   //result
  double v[NC_*4];               //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];     

  //boundary part
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  
  int site_index, internal_idx, current_idx;
  int Xbdry;

  //boundary half spinors
  double vbdry[NC_*4*Ny*Nz*Nt], vbdry_copy[NC_*4*Ny*Nz*Nt];

  // boundary part
  Xbdry = 0;
  for(int bdry_site = 0; bdry_site < Ny*Nz*Nt; ++bdry_site) {
    //write func in SiteIndex
    site_index = Xbdry + Nx*(bdry_site);
    internal_idx = 4*NC_*bdry_site;

    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
    
    for (int c = 0; c < NC_; ++c) {
      vbdry[internal_idx + 2*c          ] = vtmp[2*c        ] - vtmp[2*c+6*NC_+1 ];
      vbdry[internal_idx + 2*c + 1      ] = vtmp[2*c+1      ] + vtmp[2*c+6*NC_   ];
      vbdry[internal_idx + 2*c + 2*NC_  ] = vtmp[2*c+2*NC_  ] - vtmp[2*c+4*NC_+1 ];
      vbdry[internal_idx + 2*c + 2*NC_+1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+4*NC_   ];
    }

  }

  //Copy v1 from backward processor
  Communicator::instance()->transfer_fw(vbdry_copy,vbdry, NC_*4*Ny*Nz*Nt,0);

  Xbdry = Nx-1;
  for(int bdry_site = 0; bdry_site < Ny*Nz*Nt; ++bdry_site) {
    //write func in SiteIndex
    site_index = Xbdry + Nx*(bdry_site);
    internal_idx = 4*NC_*bdry_site;
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),0));
    res = const_cast<Field*>(&fp)->getaddr(ff_->index_r(0,0,site_index));
    
    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
      
      for (int c1 = 0; c1 < NC_; ++c1) {
	
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*vbdry_copy[internal_idx + 2*c1      ] 
		   - utmp[NC_*2*c+2*c1+1]*vbdry_copy[internal_idx + 2*c1+1    ]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*vbdry_copy[internal_idx + 2*c1      ] 
		   + utmp[NC_*2*c+2*c1  ]*vbdry_copy[internal_idx + 2*c1 +1   ]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*vbdry_copy[internal_idx + 2*c1+NC_*2] 
		   - utmp[NC_*2*c+2*c1+1]*vbdry_copy[internal_idx + 2*c1+NC_*2+1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*vbdry_copy[internal_idx + 2*c1+NC_*2] 
		   + utmp[NC_*2*c+2*c1  ]*vbdry_copy[internal_idx + 2*c1+NC_*2+1]);
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
  
  //bulk part
  for(int bulk_site = 0; bulk_site < Ny*Nz*Nt; ++bulk_site) {
    for (int ix = 0; ix < Nx-1; ++ix){
      //write func in SiteIndex
      site_index  = ix+1 + Nx*(bulk_site);
      current_idx = ix   + Nx*(bulk_site);
      
      utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(current_idx),0));
      vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
      res  = fp.getaddr(ff_->index_r(0,0,current_idx));
      //assumes matrix and fermion data are contiguous
      
      for (int c = 0; c < NC_; ++c) {
	v[2*c        ] = vtmp[2*c        ] - vtmp[2*c+6*NC_+1 ];
	v[2*c+1      ] = vtmp[2*c+1      ] + vtmp[2*c+6*NC_   ];
	v[2*c+NC_*2  ] = vtmp[2*c+2*NC_  ] - vtmp[2*c+4*NC_+1 ];
	v[2*c+NC_*2+1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+4*NC_   ];
      }
      
      for (int c = 0; c < NC_; ++c) {
	v1[c][0] = 0.0; v1[c][1] = 0.0;
	v2[c][0] = 0.0; v2[c][1] = 0.0;
	
	for (int c1 = 0; c1 < NC_; ++c1) {
	  
	  v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1      ] 
		     - utmp[NC_*2*c+2*c1+1]*v[2*c1+1    ]);
	  v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1      ] 
		     + utmp[NC_*2*c+2*c1  ]*v[2*c1 +1   ]);
	  v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2] 
		     - utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2+1]);
	  v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2] 
		     + utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2+1]);
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
}

void Dirac_Wilson::mult_yp(Field& fp, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;            //auxiliary vectors
  double* res;                   //result
  double v[NC_*4];
  double v1[NC_][2], v2[NC_][2];

  sf_up_[1]->setf(f);

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),1));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data are contiguous

    if (!sf_up_[1]->on_bdry(site)) {vtmp = sf_up_[1]->get_bulk_addr(site);}
    else {vtmp = sf_up_[1]->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v[2*c        ] = vtmp[2*c        ] + vtmp[2*c+6*NC_   ];
      v[2*c+1      ] = vtmp[2*c+1      ] + vtmp[2*c+6*NC_+1 ];
      v[2*c+NC_*2  ] = vtmp[2*c+2*NC_  ] - vtmp[2*c+4*NC_   ];
      v[2*c+NC_*2+1] = vtmp[2*c+2*NC_+1] - vtmp[2*c+4*NC_+1 ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1      ] 
		   - utmp[NC_*2*c+2*c1+1]*v[2*c1+1    ]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1      ] 
                   + utmp[NC_*2*c+2*c1  ]*v[2*c1 +1   ]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2] 
                   - utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2+1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2] 
                   + utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2+1]);
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

void Dirac_Wilson::mult_zp(Field& fp, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;                   //result
  double v[NC_*4];
  double v1[NC_][2], v2[NC_][2];      

  sf_up_[2]->setf(f);

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),2));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data is contiguous

    if (!sf_up_[2]->on_bdry(site)) {vtmp = sf_up_[2]->get_bulk_addr(site);}
    else {vtmp = sf_up_[2]->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v[2*c        ] = vtmp[2*c        ] - vtmp[2*c+4*NC_+1 ];
      v[2*c+1      ] = vtmp[2*c+1      ] + vtmp[2*c+4*NC_   ];
      v[2*c+NC_*2  ] = vtmp[2*c+2*NC_  ] + vtmp[2*c+6*NC_+1 ];
      v[2*c+NC_*2+1] = vtmp[2*c+2*NC_+1] - vtmp[2*c+6*NC_   ];
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1      ] 
		   - utmp[NC_*2*c+2*c1+1]*v[2*c1+1    ]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1      ] 
                   + utmp[NC_*2*c+2*c1  ]*v[2*c1 +1   ]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2] 
                   - utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2+1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2] 
                   + utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2+1]);
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

void Dirac_Wilson::mult_tp(Field& fp, const Field& f) const{
  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;
  double v[NC_*4];
  //double v1tmp[NC_][2], v2tmp[NC_][2]; //auxiliary vectors
  double v1[NC_][2], v2[NC_][2];       //result

  sf_up_[3]->setf(f);

  for(int site = 0; site <Nvol_; ++site){
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(site),3));
    res  = fp.getaddr(ff_->index_r(0,0,site));
    //assumes matrix and fermion data is contiguous

    if (!sf_up_[3]->on_bdry(site)) {vtmp = sf_up_[3]->get_bulk_addr(site);}
    else {vtmp = sf_up_[3]->get_bdry_addr(site);}

    for (int c = 0; c < NC_; ++c) {
      v[2*c        ] = vtmp[2*c+4*NC_  ]*2.0;
      v[2*c+1      ] = vtmp[2*c+4*NC_+1]*2.0;
      v[2*c+NC_*2  ] = vtmp[2*c+6*NC_  ]*2.0;
      v[2*c+NC_*2+1] = vtmp[2*c+6*NC_+1]*2.0;
    }

    for (int c = 0; c < NC_; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
   
      for (int c1 = 0; c1 < NC_; ++c1) {
	v1[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1      ] 
		   - utmp[NC_*2*c+2*c1+1]*v[2*c1+1    ]);
	v1[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1      ] 
                   + utmp[NC_*2*c+2*c1  ]*v[2*c1 +1   ]);
	v2[c][0] += (utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2] 
                   - utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2+1]);
	v2[c][1] += (utmp[NC_*2*c+2*c1+1]*v[2*c1+NC_*2] 
                   + utmp[NC_*2*c+2*c1  ]*v[2*c1+NC_*2+1]);
      }

      res[2*c+4*NC_  ] += v1[c][0];
      res[2*c+4*NC_+1] += v1[c][1];
      res[2*c+6*NC_  ] += v2[c][0];
      res[2*c+6*NC_+1] += v2[c][1];

    }

  }
}

void Dirac_Wilson::mult_xm(Field& w, const Field& f) const{
  double* utmp;                        // auxiliary matrix
  const double* vtmp;
  double* res;                         // pointer to result
  double v1tmp[NC_][2], v2tmp[NC_][2]; // auxiliary vectors
  double v1[NC_][2], v2[NC_][2];     
  
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  
  int site_index, internal_idx, current_idx;
  int Xbdry;
  
  //boundary half spinors
  double vbdry[NC_*4*Ny*Nz*Nt], vbdry_copy[NC_*4*Ny*Nz*Nt];
  
  // boundary part
  Xbdry = Nx -1;
  for(int bdry_site = 0; bdry_site < Ny*Nz*Nt; ++bdry_site) {
    //write func in SiteIndex
    site_index = Xbdry + Nx*(bdry_site);
    internal_idx = 4*NC_*bdry_site;
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),0));
    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
  
    for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ] + vtmp[2*c+6*NC_+1 ];
      v1tmp[c][1] = vtmp[2*c+1      ] - vtmp[2*c+6*NC_   ];
      v2tmp[c][0] = vtmp[2*c+2*NC_  ] + vtmp[2*c+4*NC_+1 ];
      v2tmp[c][1] = vtmp[2*c+2*NC_+1] - vtmp[2*c+4*NC_   ];
     }
    for (int c = 0; c < NC_; ++c) {
      vbdry[internal_idx + 2*c          ] = 0.0; 
      vbdry[internal_idx + 2*c + 1      ] = 0.0;
      vbdry[internal_idx + 2*c + 2*NC_  ] = 0.0;
      vbdry[internal_idx + 2*c + 2*NC_+1] = 0.0;

      for (int c1 = 0; c1 < NC_; ++c1) {
	vbdry[internal_idx + 2*c          ] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
					       + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	vbdry[internal_idx + 2*c + 1      ] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
					       - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	vbdry[internal_idx + 2*c + 2*NC_  ] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
					       + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	vbdry[internal_idx + 2*c + 2*NC_+1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
					       - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);
      }
    }
  }
  //Copy v1 from backward processor
  Communicator::instance()->transfer_bk(vbdry_copy,vbdry, NC_*4*Ny*Nz*Nt,0);

  Xbdry = 0;
  for(int bdry_site = 0; bdry_site < Ny*Nz*Nt; ++bdry_site) {
    //write func in SiteIndex
    site_index = Xbdry + Nx*(bdry_site);
    res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,site_index));
    //res  = &w[ff_->index_r(0,0,site_index)];
    for (int c = 0; c < NC_; ++c) {  
      res[2*c        ] +=  vbdry_copy[4*NC_*bdry_site + 2*c          ];
      res[2*c+1      ] +=  vbdry_copy[4*NC_*bdry_site + 2*c + 1      ];
      res[2*c+2*NC_  ] +=  vbdry_copy[4*NC_*bdry_site + 2*c + 2*NC_  ];
      res[2*c+2*NC_+1] +=  vbdry_copy[4*NC_*bdry_site + 2*c + 2*NC_+1];
      res[2*c+4*NC_  ] += -vbdry_copy[4*NC_*bdry_site + 2*c + 2*NC_+1];
      res[2*c+4*NC_+1] +=  vbdry_copy[4*NC_*bdry_site + 2*c + 2*NC_  ];
      res[2*c+6*NC_  ] += -vbdry_copy[4*NC_*bdry_site + 2*c + 1      ];
      res[2*c+6*NC_+1] +=  vbdry_copy[4*NC_*bdry_site + 2*c          ];
    }
  }
   

  //bulk part

  for(int bulk_site = 0; bulk_site <Ny*Nz*Nt; ++bulk_site){
    for (int ix = 1; ix < Nx; ++ix){
      //write func in SiteIndex
      site_index  = ix-1 + Nx*(bulk_site);
      current_idx = ix   + Nx*(bulk_site);
      
      utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),0));
      vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
      res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,current_idx));
      //      res  = &w[ff_->index_r(0,0,current_idx)];
      
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
	
	res[2*c        ] +=  v1[c][0];
	res[2*c+1      ] +=  v1[c][1];
	res[2*c+2*NC_  ] +=  v2[c][0];
	res[2*c+2*NC_+1] +=  v2[c][1];
	res[2*c+4*NC_  ] += -v2[c][1];
	res[2*c+4*NC_+1] +=  v2[c][0];
	res[2*c+6*NC_  ] += -v1[c][1];
	res[2*c+6*NC_+1] +=  v1[c][0];
      }
    }
  }
}


void Dirac_Wilson::mult_ym(Field& w, const Field& f) const{
  double* utmp;                        // auxiliary matrix
  const double* vtmp;
  double* res;                         // pointer to result
  double v1tmp[NC_][2], v2tmp[NC_][2]; // auxiliary vectors
  double v1[NC_][2], v2[NC_][2];     
  
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  
  int site_index, internal_idx, current_idx;
  int Ybdry;
  
  //boundary half spinors
  double vbdry[NC_*4*Nx*Nz*Nt], vbdry_copy[NC_*4*Nx*Nz*Nt];
  
  // boundary part
  Ybdry = Ny -1;
  for(int bdry_site = 0; bdry_site < Nz*Nt; ++bdry_site) {
    for (int ix = 0; ix < Nx; ++ix) {
      //write func in SiteIndex
      site_index = ix + Nx*(Ybdry + Ny*(bdry_site));
      internal_idx = 4*NC_*(ix+ Nx*bdry_site);
      utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),1));
      vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
  
      for (int c = 0; c < NC_; ++c) {
	v1tmp[c][0] = vtmp[2*c        ] - vtmp[2*c+6*NC_   ];
	v1tmp[c][1] = vtmp[2*c+1      ] - vtmp[2*c+6*NC_+1 ];
	v2tmp[c][0] = vtmp[2*c+2*NC_  ] + vtmp[2*c+4*NC_   ];
	v2tmp[c][1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+4*NC_+1 ];
	
      }
      for (int c = 0; c < NC_; ++c) {
	vbdry[internal_idx + 2*c          ] = 0.0; 
	vbdry[internal_idx + 2*c + 1      ] = 0.0;
	vbdry[internal_idx + 2*c + 2*NC_  ] = 0.0;
	vbdry[internal_idx + 2*c + 2*NC_+1] = 0.0;
	
	for (int c1 = 0; c1 < NC_; ++c1) {
	  vbdry[internal_idx + 2*c          ] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
						   + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	  vbdry[internal_idx + 2*c + 1      ] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
						   - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	  vbdry[internal_idx + 2*c + 2*NC_  ] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
						   + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	  vbdry[internal_idx + 2*c + 2*NC_+1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
						   - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);
	}
      }
    }
  }
  //Copy v1 from backward processor
  Communicator::instance()->transfer_bk(vbdry_copy,vbdry, NC_*4*Nx*Nz*Nt,1);
  
  Ybdry = 0;
  for(int bdry_site = 0; bdry_site < Nz*Nt; ++bdry_site) {
    for (int ix = 0; ix < Nx; ++ix) {
      //write func in SiteIndex
      internal_idx = 4*NC_*(ix+ Nx*bdry_site);
      site_index   = ix+ Nx*(Ybdry + Ny*(bdry_site));
      res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,site_index));
      //res  = &w[ff_->index_r(0,0,site_index)];
      for (int c = 0; c < NC_; ++c) {  
	res[2*c        ] +=  vbdry_copy[internal_idx + 2*c          ];
	res[2*c+1      ] +=  vbdry_copy[internal_idx + 2*c + 1      ];
	res[2*c+2*NC_  ] +=  vbdry_copy[internal_idx + 2*c + 2*NC_  ];
	res[2*c+2*NC_+1] +=  vbdry_copy[internal_idx + 2*c + 2*NC_+1];
	res[2*c+4*NC_  ] +=  vbdry_copy[internal_idx + 2*c + 2*NC_  ];
	res[2*c+4*NC_+1] +=  vbdry_copy[internal_idx + 2*c + 2*NC_+1];
	res[2*c+6*NC_  ] += -vbdry_copy[internal_idx + 2*c          ];
	res[2*c+6*NC_+1] += -vbdry_copy[internal_idx + 2*c + 1      ];
      }
    }
  }
  
  
  //bulk part
  
  for(int bulk_site = 0; bulk_site <Nz*Nt; ++bulk_site){
    for (int iy = 1; iy < Ny; ++iy){
      for (int ix = 0; ix < Nx; ++ix){
	//write func in SiteIndex
	site_index  = ix + Nx*(iy-1 + Ny*(bulk_site));
	current_idx = ix + Nx*(iy   + Ny*(bulk_site));
	
	utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),1));
	vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
	res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,current_idx));
	//res  = &w[ff_->index_r(0,0,current_idx)];
	
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

	  res[2*c        ] +=  v1[c][0];
	  res[2*c+1      ] +=  v1[c][1];
	  res[2*c+2*NC_  ] +=  v2[c][0];
	  res[2*c+2*NC_+1] +=  v2[c][1];
	  res[2*c+4*NC_  ] +=  v2[c][0];
	  res[2*c+4*NC_+1] +=  v2[c][1];
	  res[2*c+6*NC_  ] += -v1[c][0];
	  res[2*c+6*NC_+1] += -v1[c][1];
	  
	}
      }
    }
  }
}
    

void Dirac_Wilson::mult_zm(Field& w, const Field& f) const{
  double* utmp;                        // auxiliary matrix
  const double* vtmp;
  double* res;                         // pointer to result
  double v1tmp[NC_][2], v2tmp[NC_][2]; // auxiliary vectors
  double v1[NC_][2], v2[NC_][2];     
  
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  
  int site_index, internal_idx, current_idx;
  int Zbdry;
  
  //boundary half spinors
  double vbdry[NC_*4*Nx*Ny*Nt], vbdry_copy[NC_*4*Nx*Ny*Nt];
  
  // boundary part
  Zbdry = Nz -1;
  for(int it = 0; it < Nt; ++it) {
    for (int ixy = 0; ixy < Nx*Ny; ++ixy) {
      //write func in SiteIndex
      site_index = ixy + Zbdry*Nx*Ny + Nx*Ny*Nz*it;
      internal_idx = 4*NC_*(ixy+ Nx*Ny*it);
      utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),2));
      vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
      
      for (int c = 0; c < NC_; ++c) {
	v1tmp[c][0] = vtmp[2*c        ] + vtmp[2*c+4*NC_+1 ];
	v1tmp[c][1] = vtmp[2*c+1      ] - vtmp[2*c+4*NC_   ];
	v2tmp[c][0] = vtmp[2*c+2*NC_  ] - vtmp[2*c+6*NC_+1 ];
	v2tmp[c][1] = vtmp[2*c+2*NC_+1] + vtmp[2*c+6*NC_   ];
      }
      for (int c = 0; c < NC_; ++c) {
	vbdry[internal_idx + 2*c          ] = 0.0; 
	vbdry[internal_idx + 2*c + 1      ] = 0.0;
	vbdry[internal_idx + 2*c + 2*NC_  ] = 0.0;
	vbdry[internal_idx + 2*c + 2*NC_+1] = 0.0;
	
	for (int c1 = 0; c1 < NC_; ++c1) {
	  vbdry[internal_idx + 2*c          ] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
						   + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	  vbdry[internal_idx + 2*c + 1      ] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
						   - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	  vbdry[internal_idx + 2*c + 2*NC_  ] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
						   + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	  vbdry[internal_idx + 2*c + 2*NC_+1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
						   - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);
	}
      }
    }
  }
  //Copy v1 from backward processor
  Communicator::instance()->transfer_bk(vbdry_copy,vbdry, NC_*4*Nx*Ny*Nt,2);
  
  Zbdry = 0;
  for(int it = 0; it < Nt; ++it) {
    for (int ixy = 0; ixy < Nx*Ny; ++ixy) {
      //write func in SiteIndex
      internal_idx = 4*NC_*(ixy+ Nx*Ny*it);
      site_index = ixy + Zbdry*Nx*Ny + Nx*Ny*Nz*it;
      res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,site_index));
      //res  = &w[ff_->index_r(0,0,site_index)];
      for (int c = 0; c < NC_; ++c) {  
	res[2*c        ] +=  vbdry_copy[internal_idx + 2*c          ];
	res[2*c+1      ] +=  vbdry_copy[internal_idx + 2*c + 1      ];
	res[2*c+2*NC_  ] +=  vbdry_copy[internal_idx + 2*c + 2*NC_  ];
	res[2*c+2*NC_+1] +=  vbdry_copy[internal_idx + 2*c + 2*NC_+1];
	res[2*c+4*NC_  ] += -vbdry_copy[internal_idx + 2*c + 1      ];
	res[2*c+4*NC_+1] +=  vbdry_copy[internal_idx + 2*c          ];
	res[2*c+6*NC_  ] +=  vbdry_copy[internal_idx + 2*c + 2*NC_+1];
	res[2*c+6*NC_+1] += -vbdry_copy[internal_idx + 2*c + 2*NC_  ];
      }
    }
  }
  
  
  //bulk part
  
  for(int it = 0; it <Nt; ++it){
    for (int iz = 1; iz < Nz; ++iz){
      for (int ixy = 0; ixy < Nx*Ny; ++ixy){
	//write func in SiteIndex
	site_index  = ixy + (iz-1)*Nx*Ny + Nx*Ny*Nz*it;
	current_idx = ixy +     iz*Nx*Ny + Nx*Ny*Nz*it;
	
	utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),2));
	vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
	res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,current_idx));
	//	res  = &w[ff_->index_r(0,0,current_idx)];
	
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
	  
	  res[2*c        ] +=  v1[c][0];
	  res[2*c+1      ] +=  v1[c][1];
	  res[2*c+2*NC_  ] +=  v2[c][0];
	  res[2*c+2*NC_+1] +=  v2[c][1];
	  res[2*c+4*NC_  ] += -v1[c][1];
	  res[2*c+4*NC_+1] +=  v1[c][0];
	  res[2*c+6*NC_  ] +=  v2[c][1];
	  res[2*c+6*NC_+1] += -v2[c][0];
	  
	  
	}
      }
    }
  }
}

void Dirac_Wilson::mult_tm(Field& w, const Field& f) const{
  double* utmp;                        // auxiliary matrix
  const double* vtmp;
  double* res;                         // pointer to result
  double v1tmp[NC_][2], v2tmp[NC_][2]; // auxiliary vectors
  double v1[NC_][2], v2[NC_][2];     
  
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  
  int site_index, internal_idx, current_idx;
  int Tbdry;
  
  //boundary half spinors
  double vbdry[NC_*4*Nx*Ny*Nz], vbdry_copy[NC_*4*Nx*Ny*Nz];
  
  // boundary part
  Tbdry = Nt -1;
  for (int ixyz = 0; ixyz < Nx*Ny*Nz; ++ixyz) {
    //write func in SiteIndex
    site_index = ixyz + Nx*Ny*Nz*Tbdry;
    internal_idx = 4*NC_*ixyz;
    utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),3));
    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
    
    for (int c = 0; c < NC_; ++c) {
      v1tmp[c][0] = vtmp[2*c        ]*2.0;
      v1tmp[c][1] = vtmp[2*c+1      ]*2.0;
      v2tmp[c][0] = vtmp[2*c+2*NC_  ]*2.0;
      v2tmp[c][1] = vtmp[2*c+2*NC_+1]*2.0;
    }
    for (int c = 0; c < NC_; ++c) {
      vbdry[internal_idx + 2*c          ] = 0.0; 
      vbdry[internal_idx + 2*c + 1      ] = 0.0;
      vbdry[internal_idx + 2*c + 2*NC_  ] = 0.0;
      vbdry[internal_idx + 2*c + 2*NC_+1] = 0.0;
      
      for (int c1 = 0; c1 < NC_; ++c1) {
	vbdry[internal_idx + 2*c          ] += ( utmp[NC_*2*c1+2*c  ] *v1tmp[c1][0] 
						 + utmp[NC_*2*c1+2*c+1] *v1tmp[c1][1]);
	vbdry[internal_idx + 2*c + 1      ] -= ( utmp[NC_*2*c1+2*c+1] *v1tmp[c1][0] 
						 - utmp[NC_*2*c1+2*c  ] *v1tmp[c1][1]);
	vbdry[internal_idx + 2*c + 2*NC_  ] += ( utmp[NC_*2*c1+2*c  ] *v2tmp[c1][0] 
						 + utmp[NC_*2*c1+2*c+1] *v2tmp[c1][1]);
	vbdry[internal_idx + 2*c + 2*NC_+1] -= ( utmp[NC_*2*c1+2*c+1] *v2tmp[c1][0] 
						 - utmp[NC_*2*c1+2*c  ] *v2tmp[c1][1]);
      }
    }
  }
  
  //Copy v1 from backward processor
  Communicator::instance()->transfer_bk(vbdry_copy,vbdry, NC_*4*Nx*Ny*Nz,3);
  
  Tbdry = 0;
  for (int ixyz = 0; ixyz < Nx*Ny*Nz; ++ixyz) {
    //write func in SiteIndex
    internal_idx = 4*NC_*ixyz;
    site_index = ixyz;//Tbdry=0
    res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,site_index));
    //res  = &w[ff_->index_r(0,0,site_index)];
    for (int c = 0; c < NC_; ++c) {  
      res[2*c        ] +=  vbdry_copy[internal_idx + 2*c          ];
      res[2*c+1      ] +=  vbdry_copy[internal_idx + 2*c + 1      ];
      res[2*c+2*NC_  ] +=  vbdry_copy[internal_idx + 2*c + 2*NC_  ];
      res[2*c+2*NC_+1] +=  vbdry_copy[internal_idx + 2*c + 2*NC_+1];
    }
  }
  
  //bulk part  
  for(int it = 1; it <Nt; ++it){
    for (int ixyz = 0; ixyz < Nx*Ny*Nz; ++ixyz){
      //write func in SiteIndex
      site_index  = ixyz + (it-1)*Nx*Ny*Nz;
      current_idx = ixyz +     it*Nx*Ny*Nz;
      
      utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site_index),3));
      vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site_index));
      res = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,current_idx)); 
      
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
	
	res[2*c        ] +=  v1[c][0];
	res[2*c+1      ] +=  v1[c][1];
	res[2*c+2*NC_  ] +=  v2[c][0];
	res[2*c+2*NC_+1] +=  v2[c][1];
	
      }
    }
  }
}

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

const Field Dirac_Wilson::gamma5(const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();
  Field w(fsize_);
  double* ftmp;
  double* wtmp;
  for(int site = 0; site<Nvol_; ++site){
    ftmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
    wtmp = const_cast<Field*>(&w)->getaddr(ff_->index_r(0,0,site));
    for (int c = 0; c <Nc; ++c) {
      wtmp[ 2*c           ] = ftmp[ 2*c + 4*NC_  ];
      wtmp[ 2*c + 1       ] = ftmp[ 2*c + 4*NC_+1];
      wtmp[ 2*c + 2*NC_   ] = ftmp[ 2*c + 6*NC_  ];
      wtmp[ 2*c + 2*NC_+1 ] = ftmp[ 2*c + 6*NC_+1];
      wtmp[ 2*c + 4*NC_   ] = ftmp[ 2*c          ];
      wtmp[ 2*c + 4*NC_+1 ] = ftmp[ 2*c + 1      ];
      wtmp[ 2*c + 6*NC_   ] = ftmp[ 2*c + 2*NC_  ];
      wtmp[ 2*c + 6*NC_+1 ] = ftmp[ 2*c + 2*NC_+1];
    }
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

void Dirac_Wilson::mult_a0(Field& w,const Field& f) const{
  for(int d=0; d <Ndim_; ++d){
    (this->*mult_p[d])(w,f);
    (this->*mult_m[d])(w,f);
  }
  w *= -kpp_;
}

void Dirac_Wilson::mult_a1(Field& w, const Field& f) const{
  mult_a0(w,f);
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

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();
  SUNmat f;

  for(int mu=0; mu<Ndim_; ++mu){
    Field xie(fsize_);
  
    (this->*mult_p[mu])(xie, eta);

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

    (this->*mult_p[mu])(xz5, zt5);

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
