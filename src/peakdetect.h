#include <stdio.h>
static void peakdetect (double *spectrum, double *delta, int *l, int *window, double *threshold, int *peakind, int *counter){
  
  double curpos = (*window)*2; /* C indexing starts at 0 */
  int i = -1; /* C indexing starts at 0 */
  int j = 0;
  int flag = 0;

  /* Note: length of delta is one smaller than length of spectrum */
  while(curpos < (*l-1)){
    j = 0;
    flag = 0;
    if(spectrum[(*window)+i] > (*threshold)){
      for(j = 0; j < (*window - 1); j++){
        /* break if no ascent */ 
	if(delta[i + j + 1] < 0){
	    flag = 1;
            break;
	}
      }
	if(flag == 1){
	  i = i + 1;
          curpos = curpos + 1; 
          continue;
	}		       
      j = 0;
      /* Note: i starts from -1 */
      /* inclusive 'curpos' */
      /* break if no descent */
      /* Note: difference to original R implementation:
	 'i + window' instead of 'i + window + 1',
          '< curpos'  instead of '<= curpos' */
      for(j = (i+ (*window)); j < curpos; j++){   /* evtl. hier doch '<=' */       
	if(delta[j] > 0){ 
	  flag = 1;
	  break;
	}
      } 
      if(flag == 1){	
	i = i+1;
        curpos = curpos + 1;
        continue;
      }
      
      peakind[*counter] = *window + i; /* changed from 'i' to '*window + i' */
      *counter = *counter + 1;
      i = i+1; 
      curpos  = curpos+1;
    
    }
    else{
      i = i+1; 
      curpos  = curpos+1; 
    }
  }   
}
