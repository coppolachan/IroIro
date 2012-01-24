//----------------------------------------------------------------------
// dirac_wilson.cpp
//----------------------------------------------------------------------
#include "dirac_wilson.h"
#include "Tools/sunMatUtils.hpp"

using namespace SUNvec_utils;
using namespace std;




#ifdef IMPROVED_WILSON
void Dirac_Wilson::mult_xp(Field& fp, const Field& f) const{

  int Nd = CommonPrms::instance()->Nd();   /*!< @brief spinor dof */
  int Ndh = CommonPrms::instance()->Nd()/2;/*!< @brief half spionor dof */
  int Ndd = CommonPrms::instance()->Nd();  /*!< @brief Ndh*2 */
  int Nih = Ndd*NC_;          /*!< @brief internal dof of a half spinor */

  double* utmp;                  //auxiliary matrix
  const double* vtmp;
  double* res;                   //result
  double v[NC_*Nd];              //auxiliary vectors
  double v1[NC_][Ndh], v2[NC_][Ndh];     

  //boundary part
  int Nx = CommonPrms::instance()->Nx();
  int Vsl = SiteIndex::instance()->Vdir(0);  /*!< @brief Ny*Nz*Nt */
  int current_idx;
  int Xbdry;

  // boundary part
  double vbdry[Nih*Vsl]; 
  Xbdry = 0;
  for(int bsite=0; bsite<Vsl; ++bsite) {
    //write func in SiteIndex
    int site = Xbdry +Nx*bsite;
    int is = Nih*bsite;

    vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
    
    for(int c=0; c<NC_; ++c) {
      vbdry[is +2*c        ] = vtmp[2*c        ] -vtmp[2*(3*NC_+c)+1];
      vbdry[is +2*c+1      ] = vtmp[2*c+1      ] +vtmp[2*(3*NC_+c)  ];
      vbdry[is +2*(NC_+c)  ] = vtmp[2*(NC_+c)  ] -vtmp[2*(2*NC_+c)+1];
      vbdry[is +2*(NC_+c)+1] = vtmp[2*(NC_+c)+1] +vtmp[2*(2*NC_+c)  ];
    }
  }

  //Copy v1 from backward processor
  double vbcpy[Nih*Vsl];
  Communicator::instance()->transfer_fw(vbcpy,vbdry, Nih*Vsl,0);

  Xbdry = Nx-1;
  for(int bsite=0; bsite<Vsl; ++bsite){
    //write func in SiteIndex
    int site = Xbdry + Nx*bsite;
    int is = Nih*bsite;
    utmp=const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gm)(site),0));
    res =const_cast<Field*>(&fp)->getaddr(ff_->index_r(0,0,site));
    
    for(int c=0; c<NC_; ++c){
      v1[c][0] = 0.0; v1[c][1] = 0.0;
      v2[c][0] = 0.0; v2[c][1] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1[c][0]+=(utmp[2*(NC_*c+c1)  ]*vbcpy[is+2*c1      ] 
		  -utmp[2*(NC_*c+c1)+1]*vbcpy[is+2*c1+1    ]);
	v1[c][1]+=(utmp[2*(NC_*c+c1)+1]*vbcpy[is+2*c1      ] 
		  +utmp[2*(NC_*c+c1)  ]*vbcpy[is+2*c1 +1   ]);
	v2[c][0]+=(utmp[2*(NC_*c+c1)  ]*vbcpy[is+2*(NC_+c1)] 
		  -utmp[2*(NC_*c+c1)+1]*vbcpy[is+2*(NC_+c1)+1]);
	v2[c][1]+=(utmp[2*(NC_*c+c1)+1]*vbcpy[is+2*(NC_+c1)] 
		  +utmp[2*(NC_*c+c1)  ]*vbcpy[is+2*(NC_+c1)+1]);
      }
      res[        2*c  ] += v1[c][0];
      res[        2*c+1] += v1[c][1];
      res[2*(  NC_+c)  ] += v2[c][0];
      res[2*(  NC_+c)+1] += v2[c][1];
      res[2*(2*NC_+c)  ] += v2[c][1];
      res[2*(2*NC_+c)+1] -= v2[c][0];
      res[2*(3*NC_+c)  ] += v1[c][1];
      res[2*(3*NC_+c)+1] -= v1[c][0];
    }
  }
  
  //bulk part
  for(int bulk_site = 0; bulk_site < Vsl; ++bulk_site) {
    for(int ix = 0; ix < Nx-1; ++ix){
      //write func in SiteIndex
      int site = ix+1 + Nx*bulk_site;
      int current_idx = ix   + Nx*bulk_site;
      
      utmp = const_cast<Field*>(u_)->getaddr(gf_->index_r(0,0,(this->*gp)(current_idx),0));
      vtmp = const_cast<Field*>(&f)->getaddr(ff_->index_r(0,0,site));
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

#endif /*IMPROVED_WILSON*/

#ifndef IMPROVED_WILSON
void Dirac_Wilson::mult_xp(Field& fp, ShiftField* sfp) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site <Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c1][c2][0] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),0)];
	utmp[c1][c2][1] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),0)];
      }
    }
    bool on = sfp->on_bdry(site);
    double v1tmp[Nc][2], v2tmp[Nc][2];
    if (on) {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bdry(c,0,site) - sfp->im_on_bdry(c,3,site);
	v1tmp[c][1] = sfp->im_on_bdry(c,0,site) + sfp->re_on_bdry(c,3,site);
	v2tmp[c][0] = sfp->re_on_bdry(c,1,site) - sfp->im_on_bdry(c,2,site);
	v2tmp[c][1] = sfp->im_on_bdry(c,1,site) + sfp->re_on_bdry(c,2,site);
      } 
    } else {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bulk(c,0,site) - sfp->im_on_bulk(c,3,site);
	v1tmp[c][1] = sfp->im_on_bulk(c,0,site) + sfp->re_on_bulk(c,3,site);
	v2tmp[c][0] = sfp->re_on_bulk(c,1,site) - sfp->im_on_bulk(c,2,site);
	v2tmp[c][1] = sfp->im_on_bulk(c,1,site) + sfp->re_on_bulk(c,2,site);
      }
    }
    
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
		   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      fp.add(ff_->index_r(c,0,site),  v1[c][0]);
      fp.add(ff_->index_i(c,0,site),  v1[c][1]);
      fp.add(ff_->index_r(c,1,site),  v2[c][0]);
      fp.add(ff_->index_i(c,1,site),  v2[c][1]);
      fp.add(ff_->index_r(c,2,site),  v2[c][1]);
      fp.add(ff_->index_i(c,2,site), -v2[c][0]);
      fp.add(ff_->index_r(c,3,site),  v1[c][1]);
      fp.add(ff_->index_i(c,3,site), -v1[c][0]);
    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, ShiftField* sfp) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site <Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c1][c2][0] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),1)];
	utmp[c1][c2][1] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),1)];
      }
    }

    bool on = sfp->on_bdry(site);
    double v1tmp[Nc][2], v2tmp[Nc][2];
    if (on) {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bdry(c,0,site) + sfp->re_on_bdry(c,3,site);
	v1tmp[c][1] = sfp->im_on_bdry(c,0,site) + sfp->im_on_bdry(c,3,site);
	v2tmp[c][0] = sfp->re_on_bdry(c,1,site) - sfp->re_on_bdry(c,2,site);
	v2tmp[c][1] = sfp->im_on_bdry(c,1,site) - sfp->im_on_bdry(c,2,site);
      }
    } else { 
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bulk(c,0,site) + sfp->re_on_bulk(c,3,site);
	v1tmp[c][1] = sfp->im_on_bulk(c,0,site) + sfp->im_on_bulk(c,3,site);
	v2tmp[c][0] = sfp->re_on_bulk(c,1,site) - sfp->re_on_bulk(c,2,site);
	v2tmp[c][1] = sfp->im_on_bulk(c,1,site) - sfp->im_on_bulk(c,2,site);
      }
    }
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      fp.add(ff_->index_r(c,0,site),  v1[c][0]);
      fp.add(ff_->index_i(c,0,site),  v1[c][1]);
      fp.add(ff_->index_r(c,1,site),  v2[c][0]);
      fp.add(ff_->index_i(c,1,site),  v2[c][1]);
      fp.add(ff_->index_r(c,2,site), -v2[c][0]);
      fp.add(ff_->index_i(c,2,site), -v2[c][1]);
      fp.add(ff_->index_r(c,3,site),  v1[c][0]);
      fp.add(ff_->index_i(c,3,site),  v1[c][1]);
    }
  }
}

void Dirac_Wilson::mult_zp(Field& fp, ShiftField* sfp) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site <Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c1][c2][0] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),2)];
	utmp[c1][c2][1] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),2)];
      }
    }
    bool on = sfp->on_bdry(site);
    double v1tmp[Nc][2], v2tmp[Nc][2];
    if (on) {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bdry(c,0,site) - sfp->im_on_bdry(c,2,site);
	v1tmp[c][1] = sfp->im_on_bdry(c,0,site) + sfp->re_on_bdry(c,2,site);
	v2tmp[c][0] = sfp->re_on_bdry(c,1,site) + sfp->im_on_bdry(c,3,site);
	v2tmp[c][1] = sfp->im_on_bdry(c,1,site) - sfp->re_on_bdry(c,3,site);
      }
    } else {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bulk(c,0,site) - sfp->im_on_bulk(c,2,site);
	v1tmp[c][1] = sfp->im_on_bulk(c,0,site) + sfp->re_on_bulk(c,2,site);
	v2tmp[c][0] = sfp->re_on_bulk(c,1,site) + sfp->im_on_bulk(c,3,site);
	v2tmp[c][1] = sfp->im_on_bulk(c,1,site) - sfp->re_on_bulk(c,3,site);
      }
    }

    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }

    for (int c = 0; c < Nc; ++c) {
      fp.add(ff_->index_r(c,0,site),  v1[c][0]);
      fp.add(ff_->index_i(c,0,site),  v1[c][1]);
      fp.add(ff_->index_r(c,1,site),  v2[c][0]);
      fp.add(ff_->index_i(c,1,site),  v2[c][1]);
      fp.add(ff_->index_r(c,2,site),  v1[c][1]);
      fp.add(ff_->index_i(c,2,site), -v1[c][0]);
      fp.add(ff_->index_r(c,3,site), -v2[c][1]);
      fp.add(ff_->index_i(c,3,site),  v2[c][0]);
    }
  }
}

void Dirac_Wilson::mult_tp(Field& fp, ShiftField* sfp) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site <Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c1][c2][0] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),3)];
	utmp[c1][c2][1] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),3)];
      }
    }

    bool on = sfp->on_bdry(site);
    double v1tmp[Nc][2], v2tmp[Nc][2];
    if (on) {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bdry(c,2,site)*2.0;
	v1tmp[c][1] = sfp->im_on_bdry(c,2,site)*2.0;
	v2tmp[c][0] = sfp->re_on_bdry(c,3,site)*2.0;
	v2tmp[c][1] = sfp->im_on_bdry(c,3,site)*2.0;
      }
    } else {
      for (int c = 0; c < Nc; ++c) {
	v1tmp[c][0] = sfp->re_on_bulk(c,2,site)*2.0;
	v1tmp[c][1] = sfp->im_on_bulk(c,2,site)*2.0;
	v2tmp[c][0] = sfp->re_on_bulk(c,3,site)*2.0;
	v2tmp[c][1] = sfp->im_on_bulk(c,3,site)*2.0;
      }
    }
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      fp.add(ff_->index_r(c,2,site),  v1[c][0]);
      fp.add(ff_->index_i(c,2,site),  v1[c][1]);
      fp.add(ff_->index_r(c,3,site),  v2[c][0]);
      fp.add(ff_->index_i(c,3,site),  v2[c][1]);
    }
  }
}

void Dirac_Wilson::mult_xm(valarray<double>& w, const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site <Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c2][c1][0] =   (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),0)];
	utmp[c2][c1][1] = - (*u_)[gf_->index_i(c1,c2,(this->*gm)(site),0)];
      }
    }
    double v1tmp[Nc][2], v2tmp[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1tmp[c][0] = f[ff_->index_r(c,0,site)] + f[ff_->index_i(c,3,site)];
      v1tmp[c][1] = f[ff_->index_i(c,0,site)] - f[ff_->index_r(c,3,site)];
      v2tmp[c][0] = f[ff_->index_r(c,1,site)] + f[ff_->index_i(c,2,site)];
      v2tmp[c][1] = f[ff_->index_i(c,1,site)] - f[ff_->index_r(c,2,site)];
    }
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      w[ff_->index_r(c,0,site)] =  v1[c][0];
      w[ff_->index_i(c,0,site)] =  v1[c][1];
      w[ff_->index_r(c,1,site)] =  v2[c][0];
      w[ff_->index_i(c,1,site)] =  v2[c][1];
      w[ff_->index_r(c,2,site)] = -v2[c][1];
      w[ff_->index_i(c,2,site)] =  v2[c][0];
      w[ff_->index_r(c,3,site)] = -v1[c][1];
      w[ff_->index_i(c,3,site)] =  v1[c][0];
    }
  }
}

void Dirac_Wilson::mult_ym(valarray<double>& w, const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site <Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c2][c1][0] =   (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),1)];
	utmp[c2][c1][1] = - (*u_)[gf_->index_i(c1,c2,(this->*gm)(site),1)];
      }
    }
    double v1tmp[Nc][2], v2tmp[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1tmp[c][0] = f[ff_->index_r(c,0,site)] - f[ff_->index_r(c,3,site)];
      v1tmp[c][1] = f[ff_->index_i(c,0,site)] - f[ff_->index_i(c,3,site)];
      v2tmp[c][0] = f[ff_->index_r(c,1,site)] + f[ff_->index_r(c,2,site)];
      v2tmp[c][1] = f[ff_->index_i(c,1,site)] + f[ff_->index_i(c,2,site)];
    }
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      w[ff_->index_r(c,0,site)] =  v1[c][0];
      w[ff_->index_i(c,0,site)] =  v1[c][1];
      w[ff_->index_r(c,1,site)] =  v2[c][0];
      w[ff_->index_i(c,1,site)] =  v2[c][1];
      w[ff_->index_r(c,2,site)] =  v2[c][0];
      w[ff_->index_i(c,2,site)] =  v2[c][1];
      w[ff_->index_r(c,3,site)] = -v1[c][0];
      w[ff_->index_i(c,3,site)] = -v1[c][1];
    }
  }
}

void Dirac_Wilson::mult_zm(valarray<double>& w, const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site < Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c2][c1][0] =   (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),2)];
	utmp[c2][c1][1] = - (*u_)[gf_->index_i(c1,c2,(this->*gm)(site),2)];
      }
    }
    double v1tmp[Nc][2], v2tmp[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1tmp[c][0] = f[ff_->index_r(c,0,site)] + f[ff_->index_i(c,2,site)];
      v1tmp[c][1] = f[ff_->index_i(c,0,site)] - f[ff_->index_r(c,2,site)];
      v2tmp[c][0] = f[ff_->index_r(c,1,site)] - f[ff_->index_i(c,3,site)];
      v2tmp[c][1] = f[ff_->index_i(c,1,site)] + f[ff_->index_r(c,3,site)];
    }
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      w[ff_->index_r(c,0,site)] =  v1[c][0];
      w[ff_->index_i(c,0,site)] =  v1[c][1];
      w[ff_->index_r(c,1,site)] =  v2[c][0];
      w[ff_->index_i(c,1,site)] =  v2[c][1];
      w[ff_->index_r(c,2,site)] = -v1[c][1];
      w[ff_->index_i(c,2,site)] =  v1[c][0];
      w[ff_->index_r(c,3,site)] =  v2[c][1];
      w[ff_->index_i(c,3,site)] = -v2[c][0];
    }
  }
}

void Dirac_Wilson::mult_tm(valarray<double>& w, const Field& f) const{
  int Nc = CommonPrms::instance()->Nc();

  for(int site = 0; site < Nvol_; ++site){
    double utmp[Nc][Nc][2];

    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
	utmp[c2][c1][0] =   (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),3)];
	utmp[c2][c1][1] = - (*u_)[gf_->index_i(c1,c2,(this->*gm)(site),3)];
      }
    }
    double v1tmp[Nc][2], v2tmp[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1tmp[c][0] = f[ff_->index_r(c,0,site)]*2.0;
      v1tmp[c][1] = f[ff_->index_i(c,0,site)]*2.0;
      v2tmp[c][0] = f[ff_->index_r(c,1,site)]*2.0;
      v2tmp[c][1] = f[ff_->index_i(c,1,site)]*2.0;
    }
    double v1[Nc][2], v2[Nc][2];
    for (int c = 0; c < Nc; ++c) {
      v1[c][0] = 0.0; v1[c][1] = 0.0; v2[c][0] = 0.0; v2[c][1] = 0.0;
      for (int c1 = 0; c1 < Nc; ++c1) {
	v1[c][0] += (utmp[c][c1][0]*v1tmp[c1][0] 
                   - utmp[c][c1][1]*v1tmp[c1][1]);
	v1[c][1] += (utmp[c][c1][1]*v1tmp[c1][0] 
                   + utmp[c][c1][0]*v1tmp[c1][1]);
	v2[c][0] += (utmp[c][c1][0]*v2tmp[c1][0] 
                   - utmp[c][c1][1]*v2tmp[c1][1]);
	v2[c][1] += (utmp[c][c1][1]*v2tmp[c1][0] 
                   + utmp[c][c1][0]*v2tmp[c1][1]);
      }
    }
    for (int c = 0; c < Nc; ++c) {
      w[ff_->index_r(c,0,site)] =  v1[c][0];
      w[ff_->index_i(c,0,site)] =  v1[c][1];
      w[ff_->index_r(c,1,site)] =  v2[c][0];
      w[ff_->index_i(c,1,site)] =  v2[c][1];
      w[ff_->index_r(c,2,site)] =  0.0;
      w[ff_->index_i(c,2,site)] =  0.0;
      w[ff_->index_r(c,3,site)] =  0.0;
      w[ff_->index_i(c,3,site)] =  0.0;
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

void Dirac_Wilson::mult_a0(Field& w, const Field& f) const{

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
void Dirac_Wilson::mult_a1(Field& w, const Field& f) const{
  mult_a0(w,f);
  w += f;
}

#endif  /*no IMPROVED_WILSON*/

///////////////////////////////////////////////////////////////////////////////


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
#ifdef IMPROVED_WILSON
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
      fce.add(gf_->cslice(0,gsite,mu),f.getva());
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
      fce.add(gf_->cslice(0,gsite,mu),f.getva());
    }
  }
}
#endif

#ifndef IMPROVED_WILSON
void Dirac_Wilson::md_force_p(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmat_utils;

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();
  SUNmat f;

  for(int mu=0; mu<Ndim_; ++mu){
    Field xie(fsize_);
  
    //(this->*mult_p[mu])(xie, eta);
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
      int gsite = (this->*gp)(site);
      fce.add(gf_->cslice(0,gsite,mu),f.getva());
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

    //(this->*mult_p[mu])(xz5, zt5);
    
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
      fce.add(gf_->cslice(0,gsite,mu),f.getva());
    }
  }
}
#endif


const Field Dirac_Wilson::md_force_core(const Field& eta,const Field& zeta)const{

  Field fp(gf_->size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);

  fp *= -kpp_;
  return fp;
}

const Field Dirac_Wilson::md_force(const Field& eta,const Field& zeta)const{
  using namespace SUNmat_utils;
  GaugeField fp;
  
  SUNmat a_h;
  int gsite;

  fp.U = md_force_core(eta,zeta);
  
  for (int mu = 0; mu < Ndim_; ++mu){
    for (int site = 0; site < Nvol_; ++site){
      gsite = (this->*gp)(site);
      a_h = u(fp, gsite, mu);
      fp.U.set(gf_->cslice(0,gsite,mu),anti_hermite(a_h));
    }
  }
   
  return fp.U;
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
