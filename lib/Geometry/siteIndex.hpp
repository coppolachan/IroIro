/*!
  @file siteIndex.hpp
  @brief Defines the interface for site indexing 
  At compile time is decided the if lexicographic or even-odd.
*/
#ifndef SITEINDEX_INCLUDED
#define SITEINDEX_INCLUDED

#include "commonPrms.hpp"
#include<vector>

class SiteIndex{
private:
  // Local dimensions
  int Nx_,Ny_,Nz_,Nt_;
  int NxNy_,NxNyNz_;
  int Nvol_;

  int Nvolh_, Nxh_;

  // Global dimensions
  int Lx_,Ly_,Lz_,Lt_;
  int LxLy_,LxLyLz_;

  int Lvolh_, Lxh_;

  
  std::vector<int> Bdir_;
  std::vector<int> slsize_;
  std::vector<int> global_site_;

  // constructor
  SiteIndex();

  // hide copy constructors
  SiteIndex(const SiteIndex&){}  
  SiteIndex& operator=(const SiteIndex&);

  //void setup_lists();
  void setup_all();
  void setup_global();

public:
  static SiteIndex* instance();

  int site(int x,int y,int z,int t) const;

  // components
  int c_x(int site) const;
  int c_y(int site) const;
  int c_z(int site) const;
  int c_t(int site) const;

  int g_x(int gsite) const;
  int g_y(int gsite) const;
  int g_z(int gsite) const;
  int g_t(int gsite) const;
  
  int global_x(int x) const;
  int global_y(int y) const;
  int global_z(int z) const;
  int global_t(int t) const;

  // indices with a step forward/backward (for bulk sites)
  int p_x(int site) const;
  int p_y(int site) const;
  int p_z(int site) const;
  int p_t(int site) const;

  int m_x(int site) const;
  int m_y(int site) const;
  int m_z(int site) const;
  int m_t(int site) const;

  int slice_x(int x,int n) const;
  int slice_y(int y,int n) const;
  int slice_z(int z,int n) const;
  int slice_t(int t,int n) const;

  static int (SiteIndex::*slice_dir[])(int, int)const;

  int slsize(int x,int dir)const;

  int Bdir(int dir);
  const std::vector<int>& get_gsite() const;
  int get_gsite(int site) const;
};

#endif
