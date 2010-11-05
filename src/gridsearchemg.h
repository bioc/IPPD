#include <stdio.h>
#include <R.h>
#include <Rmath.h> 

static void gridsearchemg(double *x, double *y, double *alpha, double *sigma, double *mu, int *n, int *alphagridlength, int *sigmagridlength, int *mugridlength, double *alphawinner,
                   double *sigmawinner, double *muwinner){
  int d,i,j,k;
  double temprss;
  double currss = 10000;
  double factor = M_SQRT_PI * M_SQRT2;
  double p1,p2,p3, z;
  double alphai, sigmaj, muk;
  /* double *fitted = Calloc(sizeof(double) * (*n), double); */
  double residual;
  /* double maxfitted = 0; */

  for(i = 0; i < *alphagridlength; i++){
    for(j = 0; j < *sigmagridlength; j++){
      for(k = 0; k < *mugridlength; k++){
        alphai = alpha[i];
        sigmaj = sigma[j];
        muk = mu[k];
        temprss = 0;
	for(d = 0; d < *n; d++){
	  p2 = factor * sigmaj/alphai;
	  /* p2 = 1/alphai; */
	  p1 = exp(0.5 * sigmaj * sigmaj/(alphai * alphai) + muk/alphai - x[d]/alphai);
          z = sigmaj/alphai + muk/sigmaj - x[d]/sigmaj;
	  p3 = pnorm(z, 0, 1, 0, 0);
	  /* fitted[d] = p1 * p2 * p3;
	     if(fitted[d] > maxfitted) maxfitted = fitted[d];  */
	  residual = y[d] - p1 * p2 * p3;
	  temprss = temprss + residual*residual;
	}

        /*
	for(d = 0; d < *n; d++){
	 residual = y[d] - fitted[d]/maxfitted;
         temprss = temprss + residual*residual;
	 } */
	/*
        printf("alpha: %f \t", alphai);
        printf("sigma: %f \t", sigmaj);
        printf("mu: %f  \t", muk);
        printf("rss: %f  \n", temprss); */
        if(temprss < currss){
	  *alphawinner = alphai;
          *sigmawinner = sigmaj;
          *muwinner = muk;
          currss = temprss;   
	}    
      }
    }
  }
}
