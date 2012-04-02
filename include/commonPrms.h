/*!
 * @file commonPrms.h
 * @brief Declares and defines CommonPrms class
 */
#ifndef COMMONPRMS_INCLUDED
#define COMMONPRMS_INCLUDED
#include "macros.hpp"

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
  // global lattice size
  static int Lx_,Ly_,Lz_,Lt_,Lvol_;

  // Number of processors assigined in each direction
  static int NPEx_,NPEy_,NPEz_,NPEt_,NP_;

  // local lattice size
  static int Nx_,Ny_,Nz_,Nt_,Nvol_;

  CommonPrms(const Lattice&);
  CommonPrms(const CommonPrms&){}
  CommonPrms& operator=(const CommonPrms&);

  static void setup(const Lattice& latt);
  static CommonPrms* instance_;

  ~CommonPrms(){delete instance_;}
public:
  static CommonPrms* instance(const Lattice & latt);
  static CommonPrms* instance();

  static int Lx(){return Lx_;}
  static int Ly(){return Ly_;}
  static int Lz(){return Lz_;}
  static int Lt(){return Lt_;}
  static int Lvol(){return Lvol_;}

  static int NPEx(){return NPEx_;}
  static int NPEy(){return NPEy_;}
  static int NPEz(){return NPEz_;}
  static int NPEt(){return NPEt_;}
  static int NP(){return NP_;}

  static int Nx(){return Nx_;}
  static int Ny(){return Ny_;}
  static int Nz(){return Nz_;}
  static int Nt(){return Nt_;}
  static int Nvol(){return Nvol_;}

  static int Nc(){return NC_;}
  static int Nd(){return ND_;}
  static int Ndim(){return NDIM_;}
};

#endif
