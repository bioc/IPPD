# include <R.h>
# include <Rmath.h>

static void gaussbasis(double *block, double *x, double *amplitudes,
              double *centers, int *npositions, int *n, int *maxnpeaks,
		double *sigma, int *nonzero, double *eps, int *uppernonzero, double *lower, double *upper, double *colmax){  

int i,j,k;
double sigmatilde; /*= sqrt(*sigma) * M_SQRT1_2; */
double norm; /* = M_SQRT_PI * M_SQRT2 * sigmatilde; */
double temp;
 double colsup;
/* double lower, upper; */
int tempindex;
int kstart = 0;
 int jstart = 0;
 int flag = 1;
 int flaginner = 1;

/* Note: 'nonzero' starts from nonzero = -1 */

 for(i = 0; i < *npositions; i++){
   colsup = 0;
   sigmatilde = sqrt(sigma[i]) * M_SQRT1_2;
   norm = M_SQRT_PI * M_SQRT2 * sigmatilde;
   flag = 1;
   jstart = 0;
   /* upper = centers[i] + fence[i]; */
   /* lower = centers[i] - fence[i]; */
    for(k= kstart; k < *n; k++){
       /* exploit that we have ordered x's */ 
      if(x[k] < lower[i]) continue;
      if(x[k] > upper[i]) break;
	temp = 0;
	 flaginner = 1; 
	/* Loop over different isotopes */
	for(j = jstart; j < *maxnpeaks; j++){
          tempindex = i + (*npositions)*j; 
	  /* edge correction, see main R function */
	  if(centers[tempindex] != 0){
	    temp = temp + amplitudes[tempindex] * norm * dnorm(x[k], centers[tempindex], sigmatilde, 0); 
	  }

           if(temp > *eps){
            if(flaginner){
	      jstart = j;
              flaginner = 0;
	    } 
	   }
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
 








