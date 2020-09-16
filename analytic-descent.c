#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "WolframLibrary.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

//#define MAXDIM 500

DLLEXPORT mint WolframLibrary_getVersion( ) {
    return WolframLibraryVersion;
}

DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {
    return 0;
}

struct FitParams {
   int   nPrms;
   double a;
   double* b;
   double* c;
   double* d;
   double* prms;
};

double* compute_gradient_LOCAL(struct FitParams dat)
{  
   int   nPrms = dat.nPrms;
// if MAXDIM is not known at compile time
// then use dynamic variable arrays
#ifndef MAXDIM
  #define MAXDIM nPrms
#endif

   double ths[MAXDIM];
   double a;
   double b[MAXDIM];
   double c[MAXDIM];
   double d[MAXDIM*MAXDIM];
   double* result = calloc(nPrms+1,sizeof(double));
   
   //Copy data to stack
   a = dat.a;
   for(int i=0; i<nPrms; i++){
       ths[i] = dat.prms[i];
       b[i] = dat.b[i];
       c[i] = dat.c[i];
   }
   for(int i=0; i<nPrms*nPrms; i++){
            d[i] = dat.d[i];
   }
   
   double f0s[MAXDIM];
   double f1s[MAXDIM];
   double f2s[MAXDIM];
   for(int i=0; i<nPrms; i++){
       f0s[i] = (1 + cos(ths[i]))/2;
       f1s[i] = sin(ths[i]);
       f2s[i] = (1 - cos(ths[i]))/2;
   }
   
   //-------------------------------------------------------
   // CALCULATE ENERGY
   double aTerm, aProd;
   aProd = 1.;
   for(int i=0; i<nPrms; i++){
       aProd = aProd * f0s[i];
       }
   aTerm = aProd * a;
   
   double bTerm = 0.;
   double cTerm = 0.;
   for(int i=0; i<nPrms; i++){
       bTerm += aProd/f0s[i]*f1s[i]*b[i];
       cTerm += aProd/f0s[i]*f2s[i]*c[i];
   }
   
   double dTerm = 0.;
   for(int i=0; i<nPrms; i++){
       for(int j=i+1; j<nPrms; j++){
       dTerm += aProd/f0s[i]/f0s[j]*f1s[i]*f1s[j]*d[j+nPrms*i];
       }
   }
   
   double energy = aTerm + bTerm + cTerm + dTerm;
   result[0] = energy;
 //-------------------------------------------------------  
   
   double f0derivs[MAXDIM];
   double f1derivs[MAXDIM];
   double f2derivs[MAXDIM];
   for(int i=0; i<nPrms; i++){
       f0derivs[i] = -sin(ths[i])/2;
       f1derivs[i] = cos(ths[i]);
       f2derivs[i] = sin(ths[i])/2;
   }
   
   double aTermDerivs[MAXDIM];
   for(int i=0; i<nPrms; i++){
       aTermDerivs[i] = aProd/f0s[i]*f0derivs[i]*a;
   }
   

   double bTermDerivs[MAXDIM];
   double cTermDerivs[MAXDIM];
   for(int i=0; i<nPrms; i++){//gradient entries
       bTermDerivs[i] = 0.;
       cTermDerivs[i] = 0.;
       for(int j=0; j<nPrms; j++){//B_j index j
           if (i==j) {
             bTermDerivs[i] += aProd/f0s[i]*f1derivs[i]*b[i]; 
             cTermDerivs[i] += aProd/f0s[i]*f2derivs[i]*c[i]; 
           }
           else {
             bTermDerivs[i] += aProd/f0s[i]/f0s[j]*f1s[j]*f0derivs[i]*b[j];
             cTermDerivs[i] += aProd/f0s[i]/f0s[j]*f2s[j]*f0derivs[i]*c[j];
           }
       }
   }
   
   double dTermDerivs[MAXDIM];
   for(int i=0; i<nPrms; i++){//gradient entries
       dTermDerivs[i] = 0.;
       for(int j=0; j<nPrms; j++){
           for(int k=j+1; k<nPrms; k++){//D_{jk} index k//D_{jk} index k
               if(i==j) {
                  dTermDerivs[i] += aProd/f0s[j]/f0s[k]*f1derivs[j]*f1s[k]*d[k+nPrms*j]; 
               }
               else if (i==k) {
                   dTermDerivs[i] += aProd/f0s[j]/f0s[k]*f1s[j]*f1derivs[k]*d[k+nPrms*j];
               }
               else {
                   dTermDerivs[i] += aProd/f0s[i]/f0s[j]/f0s[k]*f0derivs[i]*f1s[j]*f1s[k]*d[k+nPrms*j];
               }
           }
        }
   }
   
   //--------------------------------
    // RETURN RESULT
    for(int i=0; i<nPrms; i++){
    result[1+i] = aTermDerivs[i] + bTermDerivs[i] + cTermDerivs[i] + dTermDerivs[i];
    }
   return(result);
    
}


DLLEXPORT int compute_gradient(WolframLibraryData libData, mint Argc,
                            MArgument *Args, MArgument Res){

     int err; // error code

     MTensor m1; // input tensor
     MTensor m2; // output tensor

     mreal *argument1; // actual data of the input tensor
     mreal *output; // data for the output tensor
     mint nPrms;

     m1 = MArgument_getMTensor(Args[0]);
     nPrms = MArgument_getMTensor(Args[1]);
     mint dimInput = libData->MTensor_getDimensions(m1);

     
     mint dimOutput[1];
     dimOutput[0] = nPrms + 1;
     err = libData->MTensor_new(MType_Real, 1, dimOutput,&m2);
     argument1 = libData->MTensor_getRealData(m1);
     output = libData->MTensor_getRealData(m2);
     
    int DoubleDataLength = nPrms*nPrms + 2*nPrms + 1 + nPrms;
	double* DoubleData = calloc(DoubleDataLength,sizeof(double));
    
    for(mint i = 0; i < DoubleDataLength; i++) {
        DoubleData[i] = argument1[i];
    }
    
    struct FitParams dat;
    dat.nPrms = (int) nPrms;
    dat.a = DoubleData[0];
    dat.b = &DoubleData[1];
    dat.c = &DoubleData[1+nPrms];
    dat.d = &DoubleData[1+2*nPrms];
    
    dat.prms = &DoubleData[nPrms*nPrms+1+2*nPrms];
    double* result = compute_gradient_LOCAL(dat);
    
    for(mint i = 0; i < nPrms + 1; i++) {
        output[i] = result[i];
    }   
     
     MArgument_setMTensor(Res, m2);
     free(DoubleData);
     free(result);
     return LIBRARY_NO_ERROR;
}
