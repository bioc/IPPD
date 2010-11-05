#include <stdio.h>

static void interpolatepeakheights(double *peakheights, int *dismass, double *mass,
			    double *endpoints, double *amplitudes, int *maxnpeaks, int *numbermass, int *maxdismass){
    int i, j, tempint1, tempint2;
    double tempdouble1, tempdouble2, diff;
    for(i = 0; i < *numbermass; i++){
        if(dismass[i] > (*maxdismass)){
	    tempint1 = dismass[i] - 2;
           for(j = 0; j < *maxnpeaks; j++){
	       peakheights[i + (*numbermass)*j] = amplitudes[tempint1 + (*maxdismass)*j];
	   }
	}
        else if(dismass[i] == 1){ 
	    tempint1 = 0;
            for(j = 0; j < *maxnpeaks; j++){
               /* printf("%f \n", amplitudes[tempint1 + (*maxdismass)*j]); */ 
	       peakheights[i + (*numbermass)*j] = amplitudes[tempint1 + (*maxdismass)*j];
	   }
	} 
        else{  
        tempint1 = dismass[i] - 2;
        tempint2 = tempint1 + 1; 
        tempdouble1 = endpoints[tempint1];
        tempdouble2 = endpoints[tempint2];
        diff = tempdouble2 - tempdouble1;  
         for(j = 0; j < *maxnpeaks; j++){
	    peakheights[i + (*numbermass)*j] = amplitudes[tempint1 + (*maxdismass)*j] + (amplitudes[tempint2 + (*maxdismass)*j] - amplitudes[tempint1 + (*maxdismass)*j])/diff * (mass[i] - tempdouble1);	  
	}
      }
    }
}
