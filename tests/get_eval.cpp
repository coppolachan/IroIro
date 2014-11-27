// For eigenvalues and eigenvectors of unitary matrix - hack C style

#include "test_smear.hpp"
#include "include/numerical_const.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <complex.h>


void exchange_eigenvalues(double* ev, int id1, int id2){
  double temp= ev[id1];
  ev[id1] = ev[id2];
  ev[id2] = temp;
}

void exchange_eigenvectors(lapack_complex_double* ev, int id1, int id2){
  for (int i = 0; i < 3; i++){
    int idx1 = 3*i+id1;
    int idx2 = 3*i+id2;
    __complex__ double temp = ev[idx1];

    ev[idx1] = ev[idx2];
    ev[idx2] = temp;
  }
}

int compared( const void* a, const void* b){
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}


int write_eig(int gx, int gy, int gz, int gt, double* matrix, double* right_vector){
  int sdim = 3;
  lapack_complex_double cmat[sdim*sdim];
  lapack_complex_double eigv[sdim];
  lapack_complex_double evec[sdim*sdim];
  lapack_complex_double dummy;
  double arg[3];

  for (int i = 0; i<9; i++){
    cmat[i] = matrix[2*i]+_Complex_I*matrix[2*i+1];
  }

  LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N','V', sdim, cmat, sdim, eigv,&dummy, sdim, evec, sdim);

  for (int i = 0; i<3; i++){
    arg[i] = carg(eigv[i]);
    printf("EV:    [ %2d %2d %2d ] %15.10g %15.10g %15.10g \n",gx, gy, gz, 
	   creal(eigv[i]),cimag(eigv[i]), arg[i]);
  }

  
  // Reorder eigenvectors and eigenvalues
  if (arg[0] < arg[1]) {
    exchange_eigenvalues(arg, 0,1);
    exchange_eigenvectors(evec, 0,1);
    
  } 
  if (arg[1] < arg[2]) {
    exchange_eigenvalues(arg, 1,2);
    exchange_eigenvectors(evec, 1,2);
    
  } 
  if (arg[0] < arg[1]) {
    exchange_eigenvalues(arg, 0,1);
    exchange_eigenvectors(evec, 0,1);
    
  }  

  //reduce to one sector
  if (arg[0] < (-PI-arg[1])) {
    double temp = -arg[0];
    arg[0] = -arg[1];
    arg[1] = temp;
    exchange_eigenvectors(evec, 0,1);
  }
  
  
  // Deep copy
  for (int i = 0; i < sdim*sdim; i++){
    right_vector[2*i  ] = creal(evec[i]);
    right_vector[2*i+1] = cimag(evec[i]);
  }
  
  

  // check eigenmodes
  
  double result[6], norm;
  for (int i = 0; i < sdim; i++){
    result[2*i] = 0.0;
    result[2*i+1] = 0.0;
    
    norm += right_vector[6*i  ]*right_vector[6*i  ]+right_vector[6*i+1  ]*right_vector[6*i+1  ];

    for (int j = 0; j <sdim; j++){
      printf("( %f  ,  %f  )  ", matrix[2*(j+3*i)], matrix[2*(j+3*i)+1]);
      printf("( %f  ,  %f  )\n  ", creal(evec[3*j]), cimag(evec[3*j]));
      result[2*i]   += matrix[2*(j+3*i)  ]*right_vector[2*(3*j)  ]- matrix[2*(j+3*i)+1]*right_vector[2*(3*j)+1  ];
      result[2*i+1] += matrix[2*(j+3*i)+1]*right_vector[2*(3*j)  ]+ matrix[2*(j+3*i)  ]*right_vector[2*(3*j)+1  ];

      //result[2*i]   += matrix[2*(j+3*i)  ]*creal(evec[3*j])- matrix[2*(j+3*i)+1]*cimag(evec[3*j]);
      //result[2*i+1] += matrix[2*(j+3*i)+1]*creal(evec[3*j])+ matrix[2*(j+3*i)  ]*cimag(evec[3*j]);
    }
    printf("\n( %f  ,  %f  )\n", result[2*i], result[2*i+1]);

    //eigenvector * lambda
    printf("\n( %f  ,  %f  )  norm: %f\n",creal(eigv[0])*creal(evec[3*i])-cimag(eigv[0])*cimag(evec[3*i]),
	   creal(eigv[0])*cimag(evec[3*i])+cimag(eigv[0])*creal(evec[3*i]), norm);
    

  }
  

  printf("THETA: [ %2d %2d %2d ] %15.10g %15.10g %15.10g \n",gx, gy, gz,  arg[0]/PI, arg[1]/PI, arg[2]/PI);
  printf("PL:    [ %2d %2d %2d ] %15.10g %15.10g\n",gx, gy, gz, (cos(arg[0])+cos(arg[1])+cos(arg[2]))/3, 
	 (sin(arg[0])+sin(arg[1])+sin(arg[2]))/3);
  
  /*
  double delta[3];
  delta[0] = fabs(arg[1]-arg[0]);
  delta[1] = fabs(arg[2]-arg[0]);
  delta[2] = fabs(arg[2]-arg[1]);
  //qsort(delta, 3, sizeof(double), compared);  
  printf("DELTA: [ %d %d %d ] %g %g %g\n",gx, gy, gz,  delta[0],delta[1],delta[2]);
  */



}

