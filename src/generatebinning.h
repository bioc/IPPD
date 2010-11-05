#include <stdio.h>
#include <R.h>

static void generatebinning(double *x, double *y, double *locnoise,
		     double *intervals, double *binstart,
                     double *binend,  
                     int *nintervals, int *n, int *flag, int *nbins){

  double noise, left, right;
  int j = 0;
  int k; 
  int indexleft, indexright;
  int ilo = (int) (*n / 2);
  int flagbin = 1;
  int endedbin = 0;
  while(j < *nintervals){
    left = intervals[j];
    right = intervals[j + (*nintervals)];
    indexleft = findInterval(x, *n, left, 0, 0, ilo, flag);
    if(*flag == -1) indexleft = 0;
    else indexleft = indexleft - 1;
    indexright = findInterval(x, *n, right, 0, 0, ilo, flag);
    if(*flag == 1) indexright = *n - 1; 
    else indexright = indexright - 1;
    noise = locnoise[indexleft];
    for(k = indexleft; k <= indexright; k++){
      if(y[k] > noise){
        flagbin = 0;
	break; 
      }
    }
    if(!flagbin){
      if(endedbin){
	*nbins = *nbins + 1;
	binstart[*nbins] = left; 
      }
      endedbin = 0;
    }
    else{
      if(!endedbin){
	endedbin = 1;
	binend[*nbins] = left;
      }
    }
     flagbin = 1;
     j = j + 1;
  }
}
