# include <R.h>
# include <Rmath.h>

static void addiagonal(int *colptr, int *rowind, double *val, 
		double *add, int *p){


  int i, j, row, colptrj;

  for(j = 0; j < *p; j++){
    colptrj = colptr[j];

    i = colptrj;

    do{
      row = rowind[i];
      if(row == j){
	val[i] = val[i] + add[j];
	break;
      }
      i = i +1;
     } while(row <= j);
  }
}
