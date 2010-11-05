#include <stdio.h>
#include <math.h>

static void mzfilter(double *positions, int *npositions, int *charge, int *filteredout){
		
double PEPTIDE_MASS_RULE_FACTOR = 0.000507f;
double PEPTIDE_MASS_RULE_BOUND = 	1./PEPTIDE_MASS_RULE_FACTOR;
double PEPTIDE_MASS_RULE_THEO_PPM_BOUND = 200;
double mass;
double correction_fac;
double old_frac_mass;
double new_mass;
double new_frac_mass;
 int j; 
 double ppms;


for(j = 0; j < *npositions; j++){
  mass = positions[j] * (*charge);
  correction_fac =  mass / PEPTIDE_MASS_RULE_BOUND;
   old_frac_mass = mass - (int)(mass);
   new_mass = ((int)(mass))* (1.+PEPTIDE_MASS_RULE_FACTOR)-(int)(correction_fac);
   new_frac_mass = new_mass - (int)(new_mass);
   if (new_frac_mass - old_frac_mass > 0.5)
	{
		new_mass -= 1.;
	}			

	if (new_frac_mass - old_frac_mass < -0.5)
	{
		new_mass += 1.;
	}
	ppms = fabs(new_mass - mass)/(0.5*(new_mass + mass))*1e6; /* getPPMs(new_mass, mz) */
	if (ppms > PEPTIDE_MASS_RULE_THEO_PPM_BOUND)
	  filteredout[j] = 1;
 }
}







