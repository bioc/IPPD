# include <R.h>
# include <Rmath.h>

static void emgbasis(double *block, double *x, double *amplitudes,
              double *centers, int *npositions, int *n, int *maxnpeaks,
              double *alpha, double *sigma, double *mu, int *nonzero, double *eps, int *uppernonzero, double *lower, double *upper, double *colmax){  

  int i,j,k;

/* Cannot be defined globally anymore */
/* double p2 = M_SQRT_PI * M_SQRT2 * (*sigma)/(*alpha); */
/* double logp1part = 0.5 * (*sigma) * (*sigma)/((*alpha) * (*alpha)) + (*mu)/(*alpha); */
/* double zpart = (*sigma)/(*alpha)  + (*mu)/(*sigma); */


  double p1, p2, logp1part, zpart, p3, z, xcenter, temp, nestedtemp, colsup;
/* double lower, upper; */
int tempindex;
int kstart = 0;
 int jstart = 0;
int flag = 1;
 int flaginner = 1;

/* Note: 'nonzero' starts from nonzero = -1 */

 for(i = 0; i < *npositions; i++){
   colsup = 0;
   p2 = M_SQRT_PI * M_SQRT2 * (sigma[i])/(alpha[i]);
   logp1part = 0.5 * sigma[i] * sigma[i]/(alpha[i] * alpha[i]) + mu[i]/alpha[i];
   zpart = sigma[i]/alpha[i]  + mu[i]/sigma[i];
   flag = 1;
   jstart = 0;
   /* upper = centers[i] + (*fence); */
   /* lower = centers[i] - (*fence); */
    
    for(k= kstart; k < *n; k++){
       /* exploit that we have ordered x's */ 
      if(x[k] < lower[i]) continue;
      if(x[k] > upper[i]) break;
	temp = 0;
	/* Loop over different isotopes */
        flaginner = 1; 
	for(j = jstart; j < *maxnpeaks; j++){
          tempindex = i + (*npositions)*j; 
	  /* edge correction, see main R function */
	  if(centers[tempindex] != 0){
	     xcenter = x[k] - centers[tempindex];
	      p1 = exp(logp1part - xcenter/alpha[i]);
	      z = zpart - xcenter/sigma[i];
	      p3 = pnorm(z, 0, 1, 0, 0);
	      nestedtemp = amplitudes[tempindex] * p1 * p2 * p3;
	      if(ISNAN(nestedtemp) || !R_FINITE(nestedtemp))
		nestedtemp = 0;
	      temp = temp + nestedtemp;
	  }
	  if(temp > *eps){
            if(flaginner){
	      jstart = j;
              flaginner = 0;
	    }
	    /*  break; */
	  }           /* printf("%f \n", temp); */
	}
	  if(temp > *eps){
	    if(flag){ 
	      kstart = k;
	      flag = 0;
	    }
	    if(temp > colsup) colsup = temp;
	   *nonzero = *nonzero + 1;
           block[*nonzero] = k + 1;
           block[*nonzero + (*uppernonzero)] = i + 1;
           block[*nonzero + (2*(*uppernonzero))] = temp;
           
	  }
	}
    colmax[i] = colsup; 
       }
    }
 








