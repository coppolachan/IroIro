/*!
 * @file commonPrms.h
 * @brief Declares and defines CommonPrms class
 */
#ifndef COMMONPRMS_INCLUDED
#define COMMONPRMS_INCLUDED
#include "macros.hpp"
#include <vector>

enum{Ndim_max_= NDIM_};

enum TopBtm {Top,Btm};

enum site_dir{XDIR,YDIR,ZDIR,TDIR};

struct Lattice{
  bool stdinput;
  int Lx, Ly, Lz, Lt;
  int NPEx, NPEy, NPEz, NPEt;
};

/*!
 * @brief Defines a static class to store several global parameters
 * Lattice size and processors topology are stored in this class
 * Color and spinor dimensions are stored too.
 * Access to private parameters is provided through public functions
 */
class CommonPrms{
private:
  static int Lvol_, Nvol_, NP_;

  // global lattice size
  static std::vector<int> Lsize_;

  // Number of processors assigined in each direction
  static std::vector<int> Nnodes_;

  // local lattice size
  static std::vector<int> Nsize_;

  CommonPrms(const Lattice&);
  CommonPrms(const CommonPrms&){}
  CommonPrms& operator=(const CommonPrms&);

  static void setup(const Lattice& latt);
  static CommonPrms* instance_;

  ~CommonPrms(){delete instance_;}
public:
  static CommonPrms* instance(const Lattice & latt);
  static CommonPrms* instance();

  static int Lx(){return Lsize_[0];}
  static int Ly(){return Lsize_[1];}
  static int Lz(){return Lsize_[2];}
  static int Lt(){return Lsize_[3];}
  static int Lvol(){return Lvol_;}

  static int NPEx(){return Nnodes_[0];}
  static int NPEy(){return Nnodes_[1];}
  static int NPEz(){return Nnodes_[2];}
  static int NPEt(){return Nnodes_[3];}
  static int NP(){return NP_;}
  static int NPE(int dir){return Nnodes_[dir];}

  static int Nx(){return Nsize_[0];}
  static int Ny(){return Nsize_[1];}
  static int Nz(){return Nsize_[2];}
  static int Nt(){return Nsize_[3];}
  static int Nvol(){return Nvol_;}

  static int Nc(){return NC_;}
  static int Nd(){return ND_;}
  static int Ndim(){return NDIM_;}
  
  static int global_size(int mu){return Lsize_[mu];}
  static int local_size(int mu){return Nsize_[mu];}
  static int node_num(int mu){return Nnodes_[mu];}
};

#endif
