#include "map_parallelscalar.hpp"
#include "siteIndex.h"

Map::Map(const int dir, const ShiftSign sign):
  sites(CommonPrms::instance()->Nvol()),
  bdry_size(SiteIndex::instance()->Vdir(dir)),
  direction(dir)
{
  SiteIndex* Sindex  = SiteIndex::instance();
  node_site_map.resize(sites);
  bdry_site_map.resize(bdry_size);
  bdry_site.resize(sites);
  
  int Block[4];
  Block[0] =          CommonPrms::instance()->Nx();
  Block[1] = Block[0]*CommonPrms::instance()->Ny();
  Block[2] = Block[1]*CommonPrms::instance()->Nz();
  Block[3] = Block[2]*CommonPrms::instance()->Nt();  
    
  for (int site = 0, bdry_idx = 0; site < sites; ++site) {
    if (sign == Forward) {
      node_site_map[site] = Sindex->x_p(site,dir);
      if(Sindex->cmp(site,dir)==Sindex->Bdir(dir)){
	bdry_site[site] = true;
      }
      if(Sindex->cmp(site,dir)==0)
	bdry_site_map[bdry_idx++] = site;
    }
    if (sign == Backward) {
      node_site_map[site] = Sindex->x_m(site,dir);
      if(Sindex->cmp(site,dir)==0) {
	bdry_site[site] = true;
      }
      if(Sindex->cmp(site,dir)==Sindex->Bdir(dir))
	bdry_site_map[bdry_idx++] = site;
    }   
  }

  
  if (sign == Forward) {
    comm_transfer = &Communicator::transfer_fw;
  } else {
    comm_transfer = &Communicator::transfer_bk;
  }
  
}  

void Map::transfer(basic_type* receive, basic_type* send, size_t datasize) const{
  (Communicator::instance()->*comm_transfer)(receive, send, datasize,direction);
}
 
