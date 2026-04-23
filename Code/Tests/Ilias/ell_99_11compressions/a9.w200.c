#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


int main(void) {



int n =  9 , i, auxa, nz;


int a1, a2, a3, a4, a5, a6, a7, a8, a9;

complex omegaPowers[9], omega, dft;

float PSD1, PSD2, PSD3, PSD4;

omega = 0.7660444431 + 0.6427876097*I;
omegaPowers[0] = 1.0 + 0.0*I;
omegaPowers[1] = omega;
for(i=2; i < n; i++) { omegaPowers[i] = omegaPowers[i-1]*omega; }



for(a1 = -11; a1 <= 11; a1+=2)
for(a2 = -11; a2 <= 11; a2+=2)
for(a3 = -11; a3 <= 11; a3+=2)
for(a4 = -11; a4 <= 11; a4+=2)
for(a5 = -11; a5 <= 11; a5+=2)
for(a6 = -11; a6 <= 11; a6+=2)
for(a7 = -11; a7 <= 11; a7+=2)
for(a8 = -11; a8 <= 11; a8+=2)
for(a9 = -11; a9 <= 11; a9+=2)
{


auxa = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9;
if(auxa != 1) {continue;}

dft = a1*omegaPowers[0] + a2*omegaPowers[1] + a3*omegaPowers[2] + a4*omegaPowers[3] + a5*omegaPowers[4] + a6*omegaPowers[5] + a7*omegaPowers[6] + a8*omegaPowers[7] + a9*omegaPowers[8];
PSD1 = creal(dft)*creal(dft) + cimag(dft)*cimag(dft);
if(PSD1 > 200) continue;

dft = a1*omegaPowers[0] + a2*omegaPowers[2] + a3*omegaPowers[4] + a4*omegaPowers[6] + a5*omegaPowers[8] + a6*omegaPowers[1] + a7*omegaPowers[3] + a8*omegaPowers[5] + a9*omegaPowers[7];
PSD2 = creal(dft)*creal(dft) + cimag(dft)*cimag(dft);
if(PSD2 > 200) continue;

dft = a1*omegaPowers[0] + a2*omegaPowers[3] + a3*omegaPowers[6] + a4*omegaPowers[0] + a5*omegaPowers[3] + a6*omegaPowers[6] + a7*omegaPowers[0] + a8*omegaPowers[3] + a9*omegaPowers[6];
PSD3 = creal(dft)*creal(dft) + cimag(dft)*cimag(dft);
if(PSD3 > 200) continue;

dft = a1*omegaPowers[0] + a2*omegaPowers[4] + a3*omegaPowers[8] + a4*omegaPowers[3] + a5*omegaPowers[7] + a6*omegaPowers[2] + a7*omegaPowers[6] + a8*omegaPowers[1] + a9*omegaPowers[5];
PSD4 = creal(dft)*creal(dft) + cimag(dft)*cimag(dft);
if(PSD4 > 200) continue;

printf("%d%d%d%d %d %d %d %d %d %d %d %d %d\n", (int)rint(PSD1), (int)rint(PSD2), (int)rint(PSD3), (int)rint(PSD4), a1, a2, a3, a4, a5, a6, a7, a8, a9);


}

return(0);
}
