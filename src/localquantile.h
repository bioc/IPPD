# include <R.h>
# include <Rmath.h>
# include <stdio.h>

static void localquantile(double *y,
		   double *sorty,
                   int *indexsort,
		   int *window,
		  int *halfwindow,
                   int *n,
                   int *flag, 
                   double *q,
                   int *qindex,
		   int *nquantiles){ 

  int i;
  int j;
  int k;
  int oldindex;
  int newindex;
  int ilo = *halfwindow;
  /* double temp[5000]; */
  /* int tempindex[5000]; */
  double *temp = Calloc(sizeof(double) * ((*window)-1), double);
  int *tempindex = Calloc(sizeof(int) * ((*window)-1), int); 
  int rhs;

  
  for(i = *window; i < *n ; i++){
    oldindex = indexsort[0];
    /* printf("oldindex: %d \n", oldindex); */
    /* build up temp */
    for(j = 0; j < oldindex; j++){
      temp[j] = sorty[j];
   }
    for(j = oldindex + 1; j < *window; j++){
      temp[j - 1] = sorty[j];
   }
    /* end: build up temp */

    /* core piece: findInterval(), and re-build sorted vector */
    /* printf("newy: %f \n", y[i]); */
    newindex = findInterval(temp, *window - 1, y[i], 0, 0, ilo, flag);         
    /* printf("newindex: %d \n", newindex); */
    /*  printf("flag: %d \n", *flag); */

    if(*flag == -1){
      sorty[0] = y[i];
      for(j = 0; j < (*window - 1); j++) 
	sorty[j + 1] = temp[j]; 
    }
    if(*flag == 1){
      sorty[*window - 1] = y[i];
      for(j = 0; j < (*window - 1); j++)
	sorty[j] = temp[j];
    }
    if(*flag == 0){
      sorty[newindex] = y[i];
      
      for(j = 0; j < newindex; j++)
	sorty[j] = temp[j];
      
      for(j = newindex + 1; j < *window; j++)
	sorty[j] = temp[j - 1];
    }

    /* 

    printf("256: %f \n", sorty[256]);
    printf("113: %f \n", sorty[113]);
    printf("13: %f \n", sorty[13]);
    printf("1113: %f \n", sorty[1113]);
    printf("5000: %f \n", sorty[5000]);
    printf("0: %f \n", sorty[0]);
    */

    /* end of core piece */

    /* 2nd piece: update indexsort */

    for(j = 0; j < (*window - 1); j++){
      rhs = indexsort[j+1];
      if(rhs < oldindex) rhs = rhs  + 1;
      if(rhs <= (newindex)) rhs = rhs - 1; 
      tempindex[j] = rhs;
    }


    indexsort[*window - 1] = newindex;
     for(j = 0; j < (*window - 1); j++)
       indexsort[j] = tempindex[j];
    
     /*
    printf("256: %d \n", indexsort[256]);
    printf("113: %d \n", indexsort[113]); 
    printf("13: %d \n", indexsort[13]);
    printf("1113: %d \n", indexsort[1113]);
    printf("5000: %d \n", indexsort[5000]);
    printf("0: %d \n", indexsort[0]); 
     */

    /* final piece: get quantiles */ 
	 
      for(k = 0; k < *nquantiles; k++){
	q[(i - *halfwindow) + (*n) * k] = sorty[qindex[k]]; 
    
    }
  }
} 
