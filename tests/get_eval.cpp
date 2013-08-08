// For eigenvalues - hack C style

#include "test_smear.hpp"
#include "include/numerical_const.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <complex.h>


lapack_logical test1(const lapack_complex_double* in){
  return 1;
}

int compared( const void* a, const void* b){
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}
int write_eig(double* matrix){
  lapack_complex_double  cmat[9];
  for (int i = 0; i<9; i++){
    cmat[i] = matrix[2*i]+_Complex_I*matrix[2*i+1];
    //printf("%g %g \n",creal(cmat[i]),cimag(cmat[i]));
  }
  int sdim = 3;
  lapack_complex_double eigv[3];
  lapack_complex_double dummy;

  LAPACKE_zgees(LAPACK_ROW_MAJOR, 'N', 'N', test1, 3, cmat, 3, &sdim, eigv, &dummy, 3);
  double arg[3];
  for (int i = 0; i<3; i++){
    arg[i] = carg(eigv[i]);
    /*
    if (arg[i] < 0)
      arg[i] = 2*PI +arg[i];
    */
    printf("EV: %g %g %g \n",creal(eigv[i]),cimag(eigv[i]), arg[i]);
  }
  printf("THETA: %g %g %g \n", arg[0]/PI, arg[1]/PI, arg[2]/PI);
  printf("PL: %g %g\n",(cos(arg[0])+cos(arg[1])+cos(arg[2]))/3, (sin(arg[0])+sin(arg[1])+sin(arg[2]))/3);
  //qsort(arg, 3, sizeof(double), compared);
  double delta[3];
  delta[0] = fabs(arg[1]-arg[0]);
  delta[1] = fabs(arg[2]-arg[0]);
  delta[2] = fabs(arg[2]-arg[1]);
  //qsort(delta, 3, sizeof(double), compared);  
  printf("DELTA: %g %g %g\n", delta[0],delta[1],delta[2]);

}

