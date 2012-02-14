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

  /// boundary part ///
  int Nx = CommonPrms::instance()->Nx();
  int Vsl = SiteIndex::instance()->Vdir(0);  /*!< @brief Ny*Nz*Nt */
  double vbd[Nih*Vsl]; 

  for(int bsite=0; bsite<Vsl; ++bsite) {
    int site = SiteIndex::instance()->bdsite(bsite,0,0);
    int is = Nih*bsite;
    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = vt[r0(c)]-vt[i3(c)]; vbd[i0(c)+is] = vt[i0(c)]+vt[r3(c)];
      vbd[r1(c)+is] = vt[r1(c)]-vt[i2(c)]; vbd[i1(c)+is] = vt[i1(c)]+vt[r2(c)];
    }
  }
  
  double vbc[Nih*Vsl];  //Copy vbd from backward processor
  Communicator::instance()->transfer_fw(vbc,vbd, Nih*Vsl,0);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     

  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,Nx-1,0);
    int is = Nih*bsite;
    int id = ff_->index_r(0,0,site);
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),0));

    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= ut[r(c,c1)]*vbc[r0(c1)+is] -ut[i(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= ut[i(c,c1)]*vbc[r0(c1)+is] +ut[r(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= ut[r(c,c1)]*vbc[r1(c1)+is] -ut[i(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= ut[i(c,c1)]*vbc[r1(c1)+is] +ut[r(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r0(c), v1r[c]);   fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);   fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c), v2i[c]);   fp.add(id+i2(c),-v2r[c]);
      fp.add(id+r3(c), v1i[c]);   fp.add(id+i3(c),-v1r[c]);
    }
  }
  /// bulk part ///
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  for(int x=0; x<Nx-1; ++x){
    for(int bsite=0; bsite<Vsl; ++bsite) {
      int site = SiteIndex::instance()->bdsite(bsite,x,0);
      int xp = SiteIndex::instance()->x_p(site,0);
      int id = ff_->index_r(0,0,site);
      
      const double* vt = f.getaddr(ff_->index_r(0,0,xp));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)] -vt[i3(c)];  v1ti[c] = vt[i0(c)] +vt[r3(c)];
	v2tr[c] = vt[r1(c)] -vt[i2(c)];  v2ti[c] = vt[i1(c)] +vt[r2(c)];
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gp)(site),0));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	  v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	  v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	  v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
	}
	fp.add(id+r0(c), v1r[c]);  fp.add(id+i0(c), v1i[c]);
	fp.add(id+r1(c), v2r[c]);  fp.add(id+i1(c), v2i[c]);
	fp.add(id+r2(c), v2i[c]);  fp.add(id+i2(c),-v2r[c]);
	fp.add(id+r3(c), v1i[c]);  fp.add(id+i3(c),-v1r[c]);
      }
    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */

  /// boundary part ///
  int Ny = CommonPrms::instance()->Ny();
  int Vsl = SiteIndex::instance()->Vdir(1); /*!< @brief Nx*Nz*Nt */
  double vbd[Nih*Vsl];

   for(int bsite=0; bsite<Vsl; ++bsite) {
    int site = SiteIndex::instance()->bdsite(bsite,0,1);
    int is = Nih*bsite;
    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = vt[r0(c)]+vt[r3(c)]; vbd[i0(c)+is] = vt[i0(c)]+vt[i3(c)];
      vbd[r1(c)+is] = vt[r1(c)]-vt[r2(c)]; vbd[i1(c)+is] = vt[i1(c)]-vt[i2(c)];
    }
   }
   
   double vbc[Nih*Vsl];  //Copy vbd from backward processor
   Communicator::instance()->transfer_fw(vbc,vbd, Nih*Vsl,1);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     

  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,Ny-1,1);
    int is = Nih*bsite;
    int id = ff_->index_r(0,0,site);
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),1));

    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= ut[r(c,c1)]*vbc[r0(c1)+is] -ut[i(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= ut[i(c,c1)]*vbc[r0(c1)+is] +ut[r(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= ut[r(c,c1)]*vbc[r1(c1)+is] -ut[i(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= ut[i(c,c1)]*vbc[r1(c1)+is] +ut[r(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r0(c), v1r[c]);   fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);   fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c),-v2r[c]);   fp.add(id+i2(c),-v2i[c]);
      fp.add(id+r3(c), v1r[c]);   fp.add(id+i3(c), v1i[c]);
    }
  }
  /// bulk part ///
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  for(int y=0; y<Ny-1; ++y){
    for(int bsite=0; bsite<Vsl; ++bsite) {
      int site = SiteIndex::instance()->bdsite(bsite,y,1);
      int yp = SiteIndex::instance()->x_p(site,1);
      int id = ff_->index_r(0,0,site);
      
      const double* vt = f.getaddr(ff_->index_r(0,0,yp));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)] +vt[r3(c)];  v1ti[c] = vt[i0(c)] +vt[i3(c)];
	v2tr[c] = vt[r1(c)] -vt[r2(c)];  v2ti[c] = vt[i1(c)] -vt[i2(c)];
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gp)(site),1));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	  v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	  v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	  v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
	}
	fp.add(id+r0(c), v1r[c]);  fp.add(id+i0(c), v1i[c]);
	fp.add(id+r1(c), v2r[c]);  fp.add(id+i1(c), v2i[c]);
	fp.add(id+r2(c),-v2r[c]);  fp.add(id+i2(c),-v2i[c]);
	fp.add(id+r3(c), v1r[c]);  fp.add(id+i3(c), v1i[c]);
      }
    }
  }
}

void Dirac_Wilson::mult_zp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */

  /// boundary part ///
  int Nz = CommonPrms::instance()->Nx();
  int Vsl = SiteIndex::instance()->Vdir(2);  /*!< @brief Nx*Ny*Nt */
  double vbd[Nih*Vsl]; 

  for(int bsite=0; bsite<Vsl; ++bsite) {
    int site = SiteIndex::instance()->bdsite(bsite,0,2);
    int is = Nih*bsite;
    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = vt[r0(c)]-vt[i2(c)]; vbd[i0(c)+is] = vt[i0(c)]+vt[r2(c)];
      vbd[r1(c)+is] = vt[r1(c)]+vt[i3(c)]; vbd[i1(c)+is] = vt[i1(c)]-vt[r3(c)];
    }
  }
  
  double vbc[Nih*Vsl];  //Copy vbd from backward processor
  Communicator::instance()->transfer_fw(vbc,vbd, Nih*Vsl,2);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     

  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,Nz-1,2);
    int is = Nih*bsite;
    int id = ff_->index_r(0,0,site);
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),2));

    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= ut[r(c,c1)]*vbc[r0(c1)+is] -ut[i(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= ut[i(c,c1)]*vbc[r0(c1)+is] +ut[r(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= ut[r(c,c1)]*vbc[r1(c1)+is] -ut[i(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= ut[i(c,c1)]*vbc[r1(c1)+is] +ut[r(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r0(c), v1r[c]);   fp.add(id+i0(c), v1i[c]);
      fp.add(id+r1(c), v2r[c]);   fp.add(id+i1(c), v2i[c]);
      fp.add(id+r2(c), v1i[c]);   fp.add(id+i2(c),-v1r[c]);
      fp.add(id+r3(c),-v2i[c]);   fp.add(id+i3(c), v2r[c]);
    }
  }
  /// bulk part ///
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  for(int z=0; z<Nz-1; ++z){
    for(int bsite=0; bsite<Vsl; ++bsite) {
      int site = SiteIndex::instance()->bdsite(bsite,z,2);
      int zp = SiteIndex::instance()->x_p(site,2);
      int id = ff_->index_r(0,0,site);
      
      const double* vt = f.getaddr(ff_->index_r(0,0,zp));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)] -vt[i2(c)];  v1ti[c] = vt[i0(c)] +vt[r2(c)];
	v2tr[c] = vt[r1(c)] +vt[i3(c)];  v2ti[c] = vt[i1(c)] -vt[r3(c)];
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gp)(site),2));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	  v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	  v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	  v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
	}
	fp.add(id+r0(c), v1r[c]);  fp.add(id+i0(c), v1i[c]);
	fp.add(id+r1(c), v2r[c]);  fp.add(id+i1(c), v2i[c]);
	fp.add(id+r2(c), v1i[c]);  fp.add(id+i2(c),-v1r[c]);
	fp.add(id+r3(c),-v2i[c]);  fp.add(id+i3(c), v2r[c]);
      }
    }
  }
}

void Dirac_Wilson::mult_tp(Field& fp, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */

  /// boundary part ///
  int Nt = CommonPrms::instance()->Nt();
  int Vsl = SiteIndex::instance()->Vdir(3);  /*!< @brief Nx*Ny*Nz */
  double vbd[Nih*Vsl]; 

  for(int bsite=0; bsite<Vsl; ++bsite) {
    int site = SiteIndex::instance()->bdsite(bsite,0,3);
    int is = Nih*bsite;
    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = vt[r2(c)]*2.0;  vbd[i0(c)+is] = vt[i2(c)]*2.0;
      vbd[r1(c)+is] = vt[r3(c)]*2.0;  vbd[i1(c)+is] = vt[i3(c)]*2.0;
    }
  }
  
  double vbc[Nih*Vsl];  //Copy vbd from backward processor
  Communicator::instance()->transfer_fw(vbc,vbd, Nih*Vsl,3);

  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];     

  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,Nt-1,3);
    int is = Nih*bsite;
    int id = ff_->index_r(0,0,site);
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),3));

    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0;
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c]+= ut[r(c,c1)]*vbc[r0(c1)+is] -ut[i(c,c1)]*vbc[i0(c1)+is];
	v1i[c]+= ut[i(c,c1)]*vbc[r0(c1)+is] +ut[r(c,c1)]*vbc[i0(c1)+is];
	v2r[c]+= ut[r(c,c1)]*vbc[r1(c1)+is] -ut[i(c,c1)]*vbc[i1(c1)+is];
	v2i[c]+= ut[i(c,c1)]*vbc[r1(c1)+is] +ut[r(c,c1)]*vbc[i1(c1)+is];
      }
      fp.add(id+r2(c), v1r[c]);   fp.add(id+i2(c), v1i[c]);
      fp.add(id+r3(c), v2r[c]);   fp.add(id+i3(c), v2i[c]);
    }
  }
  /// bulk part ///
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  for(int t=0; t<Nt-1; ++t){
    for(int bsite=0; bsite<Vsl; ++bsite) {
      int site = SiteIndex::instance()->bdsite(bsite,t,3);
      int tp = SiteIndex::instance()->x_p(site,3);
      int id = ff_->index_r(0,0,site);
      
      const double* vt = f.getaddr(ff_->index_r(0,0,tp));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r2(c)]*2.0;  v1ti[c] = vt[i2(c)]*2.0;
	v2tr[c] = vt[r3(c)]*2.0;  v2ti[c] = vt[i3(c)]*2.0;
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gp)(site),3));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0;
	v2r[c] = 0.0;  v2i[c] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	  v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	  v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	  v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
	}
	fp.add(id+r2(c), v1r[c]);  fp.add(id+i2(c), v1i[c]);
	fp.add(id+r3(c), v2r[c]);  fp.add(id+i3(c), v2i[c]);
      }
    }
  }
}

void Dirac_Wilson::mult_xm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_]; // auxiliary vectors
  int Nx = CommonPrms::instance()->Nx();
  int Vsl = SiteIndex::instance()->Vdir(0);  /*!< @brief Ny*Nz*Nt */

  /// boundary part ///
  double vbd[Nih*Vsl];

  for(int bsite=0; bsite<Vsl; ++bsite) {
    int site = SiteIndex::instance()->bdsite(bsite,Nx-1,0);
    int is = Nih*bsite;

    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    for(int c=0; c<NC_; ++c){
      v1tr[c] = vt[r0(c)] +vt[i3(c)];  v1ti[c] = vt[i0(c)] -vt[r3(c)];
      v2tr[c] = vt[r1(c)] +vt[i2(c)];  v2ti[c] = vt[i1(c)] -vt[r2(c)];
    }

    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),0));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	vbd[i0(c)+is] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	vbd[r1(c)+is] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	vbd[i1(c)+is] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
      }
    }
  }
  //Copy vbd from backward processor
  double vbc[Nih*Vsl];   
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Vsl,0);

  for(int bsite=0; bsite<Vsl; ++bsite) {
    int site = SiteIndex::instance()->bdsite(bsite,0,0);
    int is = Nih*bsite;
    int id = ff_->index_r(0,0,site);

    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]);   fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]);   fm.add(id+i1(c), vbc[i1(c)+is]);
      fm.add(id+r2(c),-vbc[i1(c)+is]);   fm.add(id+i2(c), vbc[r1(c)+is]);
      fm.add(id+r3(c),-vbc[i0(c)+is]);   fm.add(id+i3(c), vbc[r0(c)+is]);
    }
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_]; // auxiliary vectors
  for(int x=1; x<Nx; ++x){
    for(int bsite=0; bsite<Vsl; ++bsite){
      int site = SiteIndex::instance()->bdsite(bsite,x,0);
      int id = ff_->index_r(0,0,site);
      int xm = SiteIndex::instance()->x_m(site,0);
      const double* vt = f.getaddr(ff_->index_r(0,0,xm));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)] +vt[i3(c)];  v1ti[c] = vt[i0(c)] -vt[r3(c)];
	v2tr[c] = vt[r1(c)] +vt[i2(c)];  v2ti[c] = vt[i1(c)] -vt[r2(c)];
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(xm),0));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0; v1i[c] = 0.0; 
	v2r[c] = 0.0; v2i[c] = 0.0;

	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	  v1i[c] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	  v2r[c] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	  v2i[c] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
	}
	fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
	fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
	fm.add(id+r2(c),-v2i[c]);  fm.add(id+i2(c), v2r[c]);
	fm.add(id+r3(c),-v1i[c]);  fm.add(id+i3(c), v1r[c]);
      }
    }
  }
}

void Dirac_Wilson::mult_ym(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  int Ny = CommonPrms::instance()->Ny();
  int Vsl = SiteIndex::instance()->Vdir(1);

  // boundary part  
  double vbd[Nih*Vsl];

  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,Ny-1,1);
    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    for(int c=0; c<NC_; ++c){
      v1tr[c] = vt[r0(c)] -vt[r3(c)];  v1ti[c] = vt[i0(c)] -vt[i3(c)];
      v2tr[c] = vt[r1(c)] +vt[r2(c)];  v2ti[c] = vt[i1(c)] +vt[i2(c)];
    }
    int is = Nih*bsite;
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),1));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0;  vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;  vbd[i1(c)+is] = 0.0;
	
      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	vbd[i0(c)+is] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	vbd[r1(c)+is] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	vbd[i1(c)+is] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
      }
    }
  }
  //Copy v1 from backward processor
  double vbc[Nih*Vsl];
  Communicator::instance()->transfer_bk(vbc,vbd, Nih*Vsl,1);
  
  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,0,1);
    int id = ff_->index_r(0,0,site);
    int is = Nih*bsite;
    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]);   fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]);   fm.add(id+i1(c), vbc[i1(c)+is]);
      fm.add(id+r2(c), vbc[r1(c)+is]);   fm.add(id+i2(c), vbc[i1(c)+is]);
      fm.add(id+r3(c),-vbc[r0(c)+is]);   fm.add(id+i3(c),-vbc[i0(c)+is]);
    }
  }
  //bulk part
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];
  for(int y=1; y<Ny; ++y){
    for(int bsite=0; bsite<Vsl; ++bsite){
      int site = SiteIndex::instance()->bdsite(bsite,y,1);
      int ym = SiteIndex::instance()->x_m(site,1);
      int id = ff_->index_r(0,0,site);
	
      const double* vt = f.getaddr(ff_->index_r(0,0,ym));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)] -vt[r3(c)];  v1ti[c] = vt[i0(c)] -vt[i3(c)];
	v2tr[c] = vt[r1(c)] +vt[r2(c)];  v2ti[c] = vt[i1(c)] +vt[i2(c)];
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(ym),1));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0; v1i[c] = 0.0; 
	v2r[c] = 0.0; v2i[c] = 0.0;

	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	  v1i[c] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	  v2r[c] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	  v2i[c] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
	}
	fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
	fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
	fm.add(id+r2(c), v2r[c]);  fm.add(id+i2(c), v2i[c]);
	fm.add(id+r3(c),-v1r[c]);  fm.add(id+i3(c),-v1i[c]); 
      }
    }
  }
}

void Dirac_Wilson::mult_zm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_]; // auxiliary vectors
  int Nz = CommonPrms::instance()->Nz();
  int Vsl = SiteIndex::instance()->Vdir(2); /*! @brief Nx*Ny*Nt */

  // boundary part
  double vbd[Nih*Vsl];

  for(int bsite=0; bsite<Vsl; ++bsite){
    int site = SiteIndex::instance()->bdsite(bsite,Nz-1,2);
    
    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    for(int c=0; c<NC_; ++c){
      v1tr[c] = vt[r0(c)] +vt[i2(c)];  v1ti[c] = vt[i0(c)] -vt[r2(c)];
      v2tr[c] = vt[r1(c)] -vt[i3(c)];  v2ti[c] = vt[i1(c)] +vt[r3(c)];
    }
    int is = Nih*bsite;
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),2));
    for(int c=0; c<NC_; ++c){
      vbd[r0(c)+is] = 0.0; 	vbd[i0(c)+is] = 0.0;
      vbd[r1(c)+is] = 0.0;	vbd[i1(c)+is] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	vbd[r0(c)+is] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	vbd[i0(c)+is] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	vbd[r1(c)+is] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	vbd[i1(c)+is] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
      }
    }
  }
  //Copy v1 from backward processor
  double vbc[Nih*Vsl];
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Vsl,2);

  for(int bsite=0; bsite<Vsl; ++bsite){
    int is = Nih*bsite;
    int site = SiteIndex::instance()->bdsite(bsite,0,2);
    int id = ff_->index_r(0,0,site);

    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]);   fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]);   fm.add(id+i1(c), vbc[i1(c)+is]);
      fm.add(id+r2(c),-vbc[i0(c)+is]);   fm.add(id+i2(c), vbc[r0(c)+is]);
      fm.add(id+r3(c), vbc[i1(c)+is]);   fm.add(id+i3(c),-vbc[r1(c)+is]);
    }
  }
  //bulk part
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];
  for(int z=1; z<Nz; ++z){
    for(int bsite=0; bsite<Vsl; ++bsite){
      int site = SiteIndex::instance()->bdsite(bsite,z,2);
      int zm = SiteIndex::instance()->x_m(site,2);

      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(zm),2));
      const double* vt = f.getaddr(ff_->index_r(0,0,zm));
      int id = ff_->index_r(0,0,site);

      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)] +vt[i2(c)];  v1ti[c] = vt[i0(c)] -vt[r2(c)];
	v2tr[c] = vt[r1(c)] -vt[i3(c)];  v2ti[c] = vt[i1(c)] +vt[r3(c)];
      }
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0; v1i[c] = 0.0; 
	v2r[c] = 0.0; v2i[c] = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	  v1i[c] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	  v2r[c] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	  v2i[c] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
	}
	fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
	fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
	fm.add(id+r2(c),-v1i[c]);  fm.add(id+i2(c), v1r[c]);
	fm.add(id+r3(c), v2i[c]);  fm.add(id+i3(c),-v2r[c]); 
      }
    }
  }
}

void Dirac_Wilson::mult_tm(Field& fm, const Field& f) const{
  int Nih = Nd_*NC_; /*!< @brief num ob elements of a half spinor */
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_]; // auxiliary vectors
  int Nt = CommonPrms::instance()->Nt();
  int Vsl = SiteIndex::instance()->Vdir(3); /*!< @brief Nx*Ny*Nz*/

  /// boundary part ///
  double vbd[Nih*Vsl];

  for(int bsite=0; bsite<Vsl; ++bsite){
    int is = Nih*bsite;
    int site = SiteIndex::instance()->bdsite(bsite,Nt-1,3);

    const double* vt = f.getaddr(ff_->index_r(0,0,site));
    for(int c=0; c<NC_; ++c){
      v1tr[c] = vt[r0(c)]*2.0;   v1ti[c] = vt[i0(c)]*2.0;
      v2tr[c] = vt[r1(c)]*2.0;   v2ti[c] = vt[i1(c)]*2.0;
    }
    const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(site),3));
    for(int c=0; c<NC_; ++c){
      vbd[is+r0(c)] = 0.0;  vbd[is+i0(c)] = 0.0;
      vbd[is+r1(c)] = 0.0;  vbd[is+i1(c)] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	vbd[is+r0(c)] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	vbd[is+i0(c)] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	vbd[is+r1(c)] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	vbd[is+i1(c)] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
      }
    }
  }
  double vbc[Nih*Vsl];  //Copy vbd from backward processor
  Communicator::instance()->transfer_bk(vbc,vbd,Nih*Vsl,3);

  for(int bsite=0; bsite<Vsl; ++bsite){
    int is = Nih*bsite;
    int site = SiteIndex::instance()->bdsite(bsite,0,3);
    int id = ff_->index_r(0,0,site);
    for(int c=0; c<NC_; ++c){  
      fm.add(id+r0(c), vbc[r0(c)+is]); fm.add(id+i0(c), vbc[i0(c)+is]);
      fm.add(id+r1(c), vbc[r1(c)+is]); fm.add(id+i1(c), vbc[i1(c)+is]);
    }
  }
  /// bulk part ///
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int t=1; t<Nt; ++t){
    for(int bsite=0; bsite<Vsl; ++bsite){
      int site = SiteIndex::instance()->bdsite(bsite,t,3);
      int tm = SiteIndex::instance()->x_m(site,3);
      int id = ff_->index_r(0,0,site);       

      const double* vt = f.getaddr(ff_->index_r(0,0,tm));
      for(int c=0; c<NC_; ++c){
	v1tr[c] = vt[r0(c)]*2.0;   v1ti[c] = vt[i0(c)]*2.0;
	v2tr[c] = vt[r1(c)]*2.0;   v2ti[c] = vt[i1(c)]*2.0;
      }
      const double* ut = u_->getaddr(gf_->index_r(0,0,(this->*gm)(tm),3));
      for(int c=0; c<NC_; ++c){
	v1r[c] = 0.0;  v1i[c] = 0.0; 
	v2r[c] = 0.0;  v2i[c] = 0.0;

	for(int c1=0; c1<NC_; ++c1){
	  v1r[c] += ut[r(c1,c)]*v1tr[c1] +ut[i(c1,c)]*v1ti[c1];
	  v1i[c] -= ut[i(c1,c)]*v1tr[c1] -ut[r(c1,c)]*v1ti[c1];
	  v2r[c] += ut[r(c1,c)]*v2tr[c1] +ut[i(c1,c)]*v2ti[c1];
	  v2i[c] -= ut[i(c1,c)]*v2tr[c1] -ut[r(c1,c)]*v2ti[c1];
	}
	fm.add(id+r0(c), v1r[c]);  fm.add(id+i0(c), v1i[c]);
	fm.add(id+r1(c), v2r[c]);  fm.add(id+i1(c), v2i[c]);
      }
    }
  }
}

const Field Dirac_Wilson::gamma5(const Field& f) const{
  Field w(fsize_);
  
  for(int site=0; site<Nvol_; ++site){
    int id = ff_->index_r(0,0,site);
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
    int id = ff_->index_r(0,0,site);
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
    int id = ff_->index_r(0,0,site);
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

  sf_up_[0]->setf(f);
  double ut[2*NC_*NC_];
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site<Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c1,c2)] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),0)];
	ut[i(c1,c2)] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),0)];
      }
    }
    if(sf_up_[0]->on_bdry(site)){
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[0]->re_bdry(c,0,site) -sf_up_[0]->im_bdry(c,3,site);
	v1ti[c] = sf_up_[0]->im_bdry(c,0,site) +sf_up_[0]->re_bdry(c,3,site);
	v2tr[c] = sf_up_[0]->re_bdry(c,1,site) -sf_up_[0]->im_bdry(c,2,site);
	v2ti[c] = sf_up_[0]->im_bdry(c,1,site) +sf_up_[0]->re_bdry(c,2,site);
      } 
    }else{
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[0]->re_bulk(c,0,site) -sf_up_[0]->im_bulk(c,3,site);
	v1ti[c] = sf_up_[0]->im_bulk(c,0,site) +sf_up_[0]->re_bulk(c,3,site);
	v2tr[c] = sf_up_[0]->re_bulk(c,1,site) -sf_up_[0]->im_bulk(c,2,site);
	v2ti[c] = sf_up_[0]->im_bulk(c,1,site) +sf_up_[0]->re_bulk(c,2,site);
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      fp.add(ff_->index_r(c,0,site), v1r[c]);
      fp.add(ff_->index_i(c,0,site), v1i[c]);
      fp.add(ff_->index_r(c,1,site), v2r[c]);
      fp.add(ff_->index_i(c,1,site), v2i[c]);
      fp.add(ff_->index_r(c,2,site), v2i[c]);
      fp.add(ff_->index_i(c,2,site),-v2r[c]);
      fp.add(ff_->index_r(c,3,site), v1i[c]);
      fp.add(ff_->index_i(c,3,site),-v1r[c]);
    }
  }
}

void Dirac_Wilson::mult_yp(Field& fp, const Field& f) const{

  sf_up_[1]->setf(f);    
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site<Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c1,c2)] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),1)];
	ut[i(c1,c2)] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),1)];
      }
    }
    if(sf_up_[1]->on_bdry(site)){
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[1]->re_bdry(c,0,site) +sf_up_[1]->re_bdry(c,3,site);
	v1ti[c] = sf_up_[1]->im_bdry(c,0,site) +sf_up_[1]->im_bdry(c,3,site);
	v2tr[c] = sf_up_[1]->re_bdry(c,1,site) -sf_up_[1]->re_bdry(c,2,site);
	v2ti[c] = sf_up_[1]->im_bdry(c,1,site) -sf_up_[1]->im_bdry(c,2,site);
      }
    }else{ 
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[1]->re_bulk(c,0,site) +sf_up_[1]->re_bulk(c,3,site);
	v1ti[c] = sf_up_[1]->im_bulk(c,0,site) +sf_up_[1]->im_bulk(c,3,site);
	v2tr[c] = sf_up_[1]->re_bulk(c,1,site) -sf_up_[1]->re_bulk(c,2,site);
	v2ti[c] = sf_up_[1]->im_bulk(c,1,site) -sf_up_[1]->im_bulk(c,2,site);
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      fp.add(ff_->index_r(c,0,site), v1r[c]);
      fp.add(ff_->index_i(c,0,site), v1i[c]);
      fp.add(ff_->index_r(c,1,site), v2r[c]);
      fp.add(ff_->index_i(c,1,site), v2i[c]);
      fp.add(ff_->index_r(c,2,site),-v2r[c]);
      fp.add(ff_->index_i(c,2,site),-v2i[c]);
      fp.add(ff_->index_r(c,3,site), v1r[c]);
      fp.add(ff_->index_i(c,3,site), v1i[c]);
    }
  }
}

void Dirac_Wilson::mult_zp(Field& fp, const Field& f) const{

  sf_up_[2]->setf(f);
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site <Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c1,c2)] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),2)];
	ut[i(c1,c2)] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),2)];
      }
    }
    if(sf_up_[2]->on_bdry(site)){
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[2]->re_bdry(c,0,site) -sf_up_[2]->im_bdry(c,2,site);
	v1ti[c] = sf_up_[2]->im_bdry(c,0,site) +sf_up_[2]->re_bdry(c,2,site);
	v2tr[c] = sf_up_[2]->re_bdry(c,1,site) +sf_up_[2]->im_bdry(c,3,site);
	v2ti[c] = sf_up_[2]->im_bdry(c,1,site) -sf_up_[2]->re_bdry(c,3,site);
      }
    }else{
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[2]->re_bulk(c,0,site) -sf_up_[2]->im_bulk(c,2,site);
	v1ti[c] = sf_up_[2]->im_bulk(c,0,site) +sf_up_[2]->re_bulk(c,2,site);
	v2tr[c] = sf_up_[2]->re_bulk(c,1,site) +sf_up_[2]->im_bulk(c,3,site);
	v2ti[c] = sf_up_[2]->im_bulk(c,1,site) -sf_up_[2]->re_bulk(c,3,site);
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1< NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      fp.add(ff_->index_r(c,0,site), v1r[c]);
      fp.add(ff_->index_i(c,0,site), v1i[c]);
      fp.add(ff_->index_r(c,1,site), v2r[c]);
      fp.add(ff_->index_i(c,1,site), v2i[c]);
      fp.add(ff_->index_r(c,2,site), v1i[c]);
      fp.add(ff_->index_i(c,2,site),-v1r[c]);
      fp.add(ff_->index_r(c,3,site),-v2i[c]);
      fp.add(ff_->index_i(c,3,site), v2r[c]);
    }
  }
}

void Dirac_Wilson::mult_tp(Field& fp, const Field& f) const{

  sf_up_[3]->setf(f);
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site<Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c1,c2)] = (*u_)[gf_->index_r(c1,c2,(this->*gp)(site),3)];
	ut[i(c1,c2)] = (*u_)[gf_->index_i(c1,c2,(this->*gp)(site),3)];
      }
    }
    if(sf_up_[3]->on_bdry(site)){
      for (int c = 0; c < NC_; ++c){
	v1tr[c] = sf_up_[3]->re_bdry(c,2,site)*2.0;
	v1ti[c] = sf_up_[3]->im_bdry(c,2,site)*2.0;
	v2tr[c] = sf_up_[3]->re_bdry(c,3,site)*2.0;
	v2ti[c] = sf_up_[3]->im_bdry(c,3,site)*2.0;
      }
    }else{
      for(int c=0; c<NC_; ++c){
	v1tr[c] = sf_up_[3]->re_bulk(c,2,site)*2.0;
	v1ti[c] = sf_up_[3]->im_bulk(c,2,site)*2.0;
	v2tr[c] = sf_up_[3]->re_bulk(c,3,site)*2.0;
	v2ti[c] = sf_up_[3]->im_bulk(c,3,site)*2.0;
      }
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      fp.add(ff_->index_r(c,2,site), v1r[c]);
      fp.add(ff_->index_i(c,2,site), v1i[c]);
      fp.add(ff_->index_r(c,3,site), v2r[c]);
      fp.add(ff_->index_i(c,3,site), v2i[c]);
    }
  }
}

void Dirac_Wilson::mult_xm(Field& fm, const Field& f) const{

  valarray<double> w(fsize_);
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site<Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c2,c1)] =  (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),0)];
	ut[i(c2,c1)] = -(*u_)[gf_->index_i(c1,c2,(this->*gm)(site),0)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1tr[c] = f[ff_->index_r(c,0,site)] +f[ff_->index_i(c,3,site)];
      v1ti[c] = f[ff_->index_i(c,0,site)] -f[ff_->index_r(c,3,site)];
      v2tr[c] = f[ff_->index_r(c,1,site)] +f[ff_->index_i(c,2,site)];
      v2ti[c] = f[ff_->index_i(c,1,site)] -f[ff_->index_r(c,2,site)];
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      w[ff_->index_r(c,0,site)] = v1r[c];
      w[ff_->index_i(c,0,site)] = v1i[c];
      w[ff_->index_r(c,1,site)] = v2r[c];
      w[ff_->index_i(c,1,site)] = v2i[c];
      w[ff_->index_r(c,2,site)] =-v2i[c];
      w[ff_->index_i(c,2,site)] = v2r[c];
      w[ff_->index_r(c,3,site)] =-v1i[c];
      w[ff_->index_i(c,3,site)] = v1r[c];
    }
  }
  sf_dn_[0]->setf(w);
  fm += sf_dn_[0]->getva();
}

void Dirac_Wilson::mult_ym(Field& fm, const Field& f) const{

  valarray<double> w(fsize_);
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site = 0; site<Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c2,c1)] =  (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),1)];
	ut[i(c2,c1)] = -(*u_)[gf_->index_i(c1,c2,(this->*gm)(site),1)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1tr[c] = f[ff_->index_r(c,0,site)] -f[ff_->index_r(c,3,site)];
      v1ti[c] = f[ff_->index_i(c,0,site)] -f[ff_->index_i(c,3,site)];
      v2tr[c] = f[ff_->index_r(c,1,site)] +f[ff_->index_r(c,2,site)];
      v2ti[c] = f[ff_->index_i(c,1,site)] +f[ff_->index_i(c,2,site)];
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      w[ff_->index_r(c,0,site)] = v1r[c];
      w[ff_->index_i(c,0,site)] = v1i[c];
      w[ff_->index_r(c,1,site)] = v2r[c];
      w[ff_->index_i(c,1,site)] = v2i[c];
      w[ff_->index_r(c,2,site)] = v2r[c];
      w[ff_->index_i(c,2,site)] = v2i[c];
      w[ff_->index_r(c,3,site)] =-v1r[c];
      w[ff_->index_i(c,3,site)] =-v1i[c];
    }
  }
  sf_dn_[1]->setf(w);
  fm += sf_dn_[1]->getva();
}

void Dirac_Wilson::mult_zm(Field& fm, const Field& f) const{

  valarray<double> w(fsize_);
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site<Nvol_; ++site){
    for(int c1=0; c1<NC_; ++c1){
      for(int c2=0; c2<NC_; ++c2){
	ut[r(c2,c1)] =  (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),2)];
	ut[i(c2,c1)] = -(*u_)[gf_->index_i(c1,c2,(this->*gm)(site),2)];
      }
    }
    for(int c=0; c<NC_; ++c){
      v1tr[c] = f[ff_->index_r(c,0,site)] +f[ff_->index_i(c,2,site)];
      v1ti[c] = f[ff_->index_i(c,0,site)] -f[ff_->index_r(c,2,site)];
      v2tr[c] = f[ff_->index_r(c,1,site)] -f[ff_->index_i(c,3,site)];
      v2ti[c] = f[ff_->index_i(c,1,site)] +f[ff_->index_r(c,3,site)];
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;

      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      w[ff_->index_r(c,0,site)] = v1r[c];
      w[ff_->index_i(c,0,site)] = v1i[c];
      w[ff_->index_r(c,1,site)] = v2r[c];
      w[ff_->index_i(c,1,site)] = v2i[c];
      w[ff_->index_r(c,2,site)] =-v1i[c];
      w[ff_->index_i(c,2,site)] = v1r[c];
      w[ff_->index_r(c,3,site)] = v2i[c];
      w[ff_->index_i(c,3,site)] =-v2r[c];
    }
  }
  sf_dn_[2]->setf(w);
  fm += sf_dn_[2]->getva();
}

void Dirac_Wilson::mult_tm(Field& fm, const Field& f) const{

  valarray<double> w(fsize_);
  double ut[2*NC_*NC_];  
  double v1tr[NC_], v1ti[NC_], v2tr[NC_], v2ti[NC_];
  double v1r[NC_], v1i[NC_], v2r[NC_], v2i[NC_];

  for(int site=0; site<Nvol_; ++site){
    for (int c1=0; c1<NC_; ++c1){
      for (int c2=0; c2<NC_; ++c2){
	ut[r(c2,c1)] =  (*u_)[gf_->index_r(c1,c2,(this->*gm)(site),3)];
	ut[i(c2,c1)] = -(*u_)[gf_->index_i(c1,c2,(this->*gm)(site),3)];
      }
    }
    for (int c=0; c<NC_; ++c){
      v1tr[c] = f[ff_->index_r(c,0,site)]*2.0;
      v1ti[c] = f[ff_->index_i(c,0,site)]*2.0;
      v2tr[c] = f[ff_->index_r(c,1,site)]*2.0;
      v2ti[c] = f[ff_->index_i(c,1,site)]*2.0;
    }
    for(int c=0; c<NC_; ++c){
      v1r[c] = 0.0; v1i[c] = 0.0; 
      v2r[c] = 0.0; v2i[c] = 0.0;
      
      for(int c1=0; c1<NC_; ++c1){
	v1r[c] += ut[r(c,c1)]*v1tr[c1] -ut[i(c,c1)]*v1ti[c1];
	v1i[c] += ut[i(c,c1)]*v1tr[c1] +ut[r(c,c1)]*v1ti[c1];
	v2r[c] += ut[r(c,c1)]*v2tr[c1] -ut[i(c,c1)]*v2ti[c1];
	v2i[c] += ut[i(c,c1)]*v2tr[c1] +ut[r(c,c1)]*v2ti[c1];
      }
    }
    for(int c=0; c<NC_; ++c){
      w[ff_->index_r(c,0,site)] = v1r[c];
      w[ff_->index_i(c,0,site)] = v1i[c];
      w[ff_->index_r(c,1,site)] = v2r[c];
      w[ff_->index_i(c,1,site)] = v2i[c];
      w[ff_->index_r(c,2,site)] = 0.0;
      w[ff_->index_i(c,2,site)] = 0.0;
      w[ff_->index_r(c,3,site)] = 0.0;
      w[ff_->index_i(c,3,site)] = 0.0;
    }
  }
  sf_dn_[3]->setf(w);
  fm += sf_dn_[3]->getva();
}

const Field Dirac_Wilson::gamma5(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_->index_r(c,0,site), f[ff_->index_r(c,2,site)]);
      w.set(ff_->index_i(c,0,site), f[ff_->index_i(c,2,site)]);
      w.set(ff_->index_r(c,1,site), f[ff_->index_r(c,3,site)]);
      w.set(ff_->index_i(c,1,site), f[ff_->index_i(c,3,site)]);
      w.set(ff_->index_r(c,2,site), f[ff_->index_r(c,0,site)]);
      w.set(ff_->index_i(c,2,site), f[ff_->index_i(c,0,site)]);
      w.set(ff_->index_r(c,3,site), f[ff_->index_r(c,1,site)]);
      w.set(ff_->index_i(c,3,site), f[ff_->index_i(c,1,site)]);
    }
  }
  return w;
}

const Field Dirac_Wilson::proj_p(const Field& f) const{
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
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
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
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

const Field Dirac_Wilson::
md_force(const Field& eta,const Field& zeta)const{
  
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
