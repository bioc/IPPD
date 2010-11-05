#include <stdio.h>
#include <math.h>

static void multiply(int *colptr, int *rowind, double *val, int *colptrphi,
              int *rowindphi, double *valphi, double *vec, int *p){

  int j, jprime, k, i, l, m;
  int flag, flag2; 

  int colptrj;
  int colptrjplusone;
  int colptrphij;
  int colptrphijplusone;
  int colptrphijprime;
  int colptrphijprimeplusone;
  int rowindphij;
  int rowindphijprime;
  int rowl;
  int rowptrj;
  
  double temp;

  for(j=0; j < *p; j++){
    colptrj = colptr[j];
    colptrjplusone = colptr[j+1];
    for(k = colptrj; k < colptrjplusone; k++){
      jprime = rowind[k];

      colptrphij = colptrphi[j];
      colptrphijplusone = colptrphi[j + 1];
      rowindphij = rowindphi[colptrphij];
      
      colptrphijprime = colptrphi[jprime];
      colptrphijprimeplusone = colptrphi[jprime + 1];
      rowindphijprime = rowindphi[colptrphijprime]; 
     
      if(rowindphijprime <= rowindphij){ 
     
      flag = 1;
      i = colptrphijprime;
      while(flag){
	if(rowindphi[i] < rowindphij){
	i = i + 1;
        continue;
	}
	else flag = 0;
      }

      flag2 = 1;
      l = 0;
      m = 0;
      temp = 0;
         
      while(flag2){
	if(((colptrphij + m) == colptrphijplusone) | ((i + l) == colptrphijprimeplusone))
	  flag2 = 0;
        else {
	  rowl = rowindphi[i + l]; 
          rowptrj = rowindphi[colptrphij + m];
          if(rowl == rowptrj){ 
	  temp = temp + valphi[i + l] * valphi[colptrphij + m] * vec[rowl];
          l = l + 1;
	  m = m + 1;
	  }
	  else{
	    while(rowl < rowptrj){
		l = l + 1;
	        if( (i + l) == colptrphijprimeplusone){
                  flag2 = 0; 
		  break;
		}
                rowl = rowindphi[i + l]; 
                if(rowl == rowptrj){
		  temp = temp + valphi[i + l] * valphi[colptrphij + m] * vec[rowl];
                  l = l + 1;
		  m = m + 1;
                  break;
		}
              }
	    while(rowptrj < rowl){
	      m = m + 1;
	      if((colptrphij + m) == colptrphijplusone){
		flag2 = 0;
		break;
	      }
	      rowptrj = rowindphi[colptrphij + m];
	      if(rowptrj == rowl){
		temp = temp + valphi[i + l] * valphi[colptrphij + m] * vec[rowl];
                l = l + 1;
		m = m + 1;
                break;
	      }
	     }
	  }
      }
      }
      }   
    else{
      flag = 1;
      i = colptrphij;
      while(flag){
	if(rowindphi[i] < rowindphijprime){
	i = i + 1;
        continue;
	}
	else flag = 0;
        
	  }

      flag2 = 1;
      l = 0;
      m = 0;
      temp = 0;

      while(flag2){
	if(((colptrphijprime + m) == colptrphijprimeplusone) | ((i + l) == colptrphijplusone))
	  flag2 = 0;
        else {
	  rowl = rowindphi[i + l]; 
	  rowptrj = rowindphi[colptrphijprime + m];
          if(rowl == rowptrj){ 
	  temp = temp + valphi[i + l] * valphi[colptrphijprime + m] * vec[rowl];
          l = l + 1;
	  m = m + 1;
	  }
	  else{
	    while(rowl < rowptrj){
		l = l + 1;
	        if( (i + l) == colptrphijplusone){
                  flag2 = 0; 
		  break;
		}
                rowl = rowindphi[i + l]; 
                if(rowl == rowptrj){
		  temp = temp + valphi[i + l] * valphi[colptrphijprime + m] * vec[rowl];
                  l = l + 1;
		  m = m + 1;
                  break;
		}
	    }
	    while(rowptrj < rowl){
	      m = m + 1;
	      if((colptrphij + m) == colptrphijprimeplusone){
		flag2 = 0;
		break;
	      }
	      rowptrj = rowindphi[colptrphijprime + m];
	      if(rowptrj == rowl){
		temp = temp + valphi[i + l] * valphi[colptrphijprime + m] * vec[rowl];
                l = l + 1;
		m = m + 1;
		break;
	      }
	     }
	  }
	}
      }
     }
    
      val[k] = temp;
    }
  }
}  
