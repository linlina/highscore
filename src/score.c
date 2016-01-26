//
//  score.c
//  
//
//  Created by Lina on 2015-04-25.
//
//
// C code derived from gconcord code (Khare et al. [2015]) with modifications
#include "score.h"
#include <sys/param.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <math.h>

// #endif /* defined(____Score_own__) */

inline double shrink(double a, double b) {
    if (b < fabs(a)) {
        if (a > 0) return(a-b);
        else       return(a+b);
    } else {
        return(0.0);
    }
}

inline int lindx(int r, int c, int p) {
    
    int rr,cc;
    
    if (r < c){ rr = r; cc = c; }
    else      { rr = c; cc = r; }
    
    return( (p*rr)+cc - ((rr*(rr+1))/2) );
    
}


/* for p < n case, gaussian data */
void Score1(int *pIn, double *S, double *lambdaIn, double *K,
            double *tol, int *maxit) {
    
    int i,j,k,r;
    int p = *pIn;
    double lambda = *lambdaIn; 
    int tp = p*(p+1)/2;
    int converged=0;
    double sSum,s1,s2;
    double maxdiff=1.;
    
    double *oldK;
    oldK = (double *)malloc(tp*sizeof(double));
    if (oldK == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    r = 0;
    while ((r<*maxit) && !converged) {
        maxdiff = 0;
        
        // update off-diagonal elements of K
        for (i=0; i<(p-1); i++){
            for (j=i+1; j<p; j++){
                s1 = 0;
                s2 = 0;
                for (k=0; k<p; k++){
                    s1 += K[i*p+k]*S[j*p+k];
                    s2 += K[k*p+j]*S[i*p+k];
                }
                s1 -= K[i*p+j]*S[j*p+j];
                s2 -= K[i*p+j]*S[i*p+i];
                
                sSum = S[i*p+i] + S[j*p+j];
                
                K[i*p+j] = shrink(-(s1+s2)/sSum, 2*lambda/sSum);
                K[j*p+i] = K[i*p+j];
            
                maxdiff +=  2*fabs(oldK[lindx(i,j,p)]-K[i*p+j]);
                oldK[lindx(i,j,p)] = K[i*p+j];
            }
        }
        
        // update diagonal elements of K
        for (i=0; i<p; i++){
            s1 = 0;
            for (k=0; k<p; k++){
                s1 += K[i*p+k]*S[i*p+k];
            }
            s1 -= K[i*p+i]*S[i*p+i];
            K[i*p+i] = (1-s1)/S[i*p+i];
            maxdiff += fabs(oldK[lindx(i,i,p)]-K[i*p+i]);
            oldK[lindx(i,i,p)] = K[i*p+i];
            
        }
        
        // check convergence
        if (maxdiff<*tol){
            converged = 1;
        }
        r++;
        
    }
    *maxit = r; // pass it back as number of iterations
    free(oldK);
}

/* using the residual trick; for n < p case, gaussian data */
void Score2(int *nIn, int *pIn, double *X, double *S, double *lambdaIn, double *K,
            double *tol, int *maxit) {
    
    int i,j,r;
    int one = 1;
    int n = *nIn;
    int p = *pIn;
    double lambda = *lambdaIn; 
    int np = n*p;
    int tp = p*(p+1)/2;
    int converged=0;
    double sSum,s1,s2;
    double maxdiff=1.;
    double tmpdiff=1.;
    
    double *oldK;
    oldK = (double *)malloc(tp*sizeof(double));
    if (oldK == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    for (i=0; i<tp; i++){
        oldK[i] = 0;
    }
    
    for (i=0; i<p; i++){
        oldK[lindx(i,i,p)] = 1.0;
    }
    
    double *resid;
    resid = (double *)malloc(np*sizeof(double));
    if (oldK == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    F77_NAME(dcopy)(&np, X, &one, resid, &one);
    r = 0;
    
    while ((r<*maxit) && !converged) {
        maxdiff = 0;
        
        // update off-diagonal elements of K
        for (i=0; i<(p-1); i++){
            for (j=i+1; j<p; j++){
                s1 = -K[i*p+j]*S[j*p+j] +
                K[i*p+i] * F77_NAME(ddot)(&n, &X[j*n], &one, &resid[i*n], &one)/n;
                s2 = -K[j*p+i]*S[i*p+i] +
                K[j*p+j] * F77_NAME(ddot)(&n, &X[i*n], &one, &resid[j*n], &one)/n;
                sSum = S[i*p+i]+S[j*p+j];
                K[i*p+j] = shrink(-(s1+s2)/sSum, 2*lambda/sSum );
                K[j*p+i] = K[i*p+j];
                tmpdiff = K[i*p+j] - oldK[lindx(i,j,p)];
                if (tmpdiff != 0) {
                    // update residuals
                    sSum = tmpdiff/K[i*p+i];
                    F77_NAME(daxpy)(&n, &sSum, &X[j*n], &one, &resid[i*n], &one);
                    sSum = tmpdiff/K[j*p+j];
                    F77_NAME(daxpy)(&n, &sSum, &X[i*n], &one, &resid[j*n], &one);
                    maxdiff +=  2*fabs(oldK[lindx(i,j,p)]-K[i*p+j]);
                    oldK[lindx(i,j,p)] = K[i*p+j];
                }
            }
        }
        
        // update diagonal elements of K
        for (i=0; i<p; i++){
            s1 = -K[i*p+i]*S[i*p+i] +
            K[i*p+i] * F77_NAME(ddot)(&n, &X[i*n], &one, &resid[i*n], &one)/n;
            K[i*p+i] = (1-s1)/S[i*p+i];
            sSum = oldK[lindx(i,i,p)]/K[i*p+i];
            F77_NAME(dscal)(&n, &sSum, &resid[i*n], &one);
            sSum = 1 - sSum;
            F77_NAME(daxpy)(&n, &sSum, &X[i*n], &one, &resid[i*n], &one);
            maxdiff += fabs(oldK[lindx(i,i,p)]-K[i*p+i]);
            oldK[lindx(i,i,p)] = K[i*p+i];
        }
        
        // check convergence
        if (maxdiff<*tol){
            converged = 1;
        }
        r++;
        
    }
    *maxit = r; // pass it back as number of iterations
    free(oldK);
    
}


void ScoreNN1(int *pIn, double *S1, double *S2, double *lambdaIn, double *K, double *tol, int *maxit) {
    
    int i,j,k,r;
    int p = *pIn;
    double lambda = *lambdaIn;
    int tp = p*(p+1)/2;
    int converged=0;
    double sSum, s1, s2;
    double maxdiff=1.;
    
    double *oldK;
    oldK = (double *)malloc(tp*sizeof(double));
    if (oldK == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    r = 0;    
    while ((r<*maxit) && !converged) {
        maxdiff = 0;
        
        // update off-diagonal elemetns of K
        for (i=0; i<(p-1); i++){
            for (j=i+1;j<p;j++) {
                s1 = 0;
                s2 = 0;
                for (k=0;k<p;k++) {
                    s1 += K[i*p+k]*S2[i*(p*p)+j*p+k];
                    s2 += K[j*p+k]*S2[j*(p*p)+i*p+k];
                }
                s1 = 2*S1[i*p+j] - s1 + S2[i*(p*p)+j*p+j] * K[i*p+j];
                s2 = 2*S1[i*p+j] - s2 + S2[j*(p*p)+i*p+i] * K[i*p+j];
                sSum = S2[i*(p*p)+j*p+j] + S2[j*(p*p)+i*p+i];
                K[i*p+j] = shrink((s1+s2)/sSum, 2*lambda/sSum);
                K[j*p+i] = K[i*p + j];
                maxdiff += 2*fabs(oldK[lindx(i, j, p)] - K[i*p+j]);
                oldK[lindx(i, j, p)] = K[i*p + j];
            }
        }
        // update diagonal entries
        for (i=0; i<p; i++) {
            s1 = 0;
            sSum = S2[i*(p*p)+i*p+i];
            for (k = 0; k<p; k++) {
                s1 += K[i*p+k]*S2[i*(p*p)+i*p+k];
            }
            s1 = 3*S1[i*p+i] - s1 + K[i*p+i]*S2[i*(p*p)+i*p+i];
            K[i*p + i] = s1/sSum;
            
            maxdiff += fabs(oldK[lindx(i,i,p)]-K[i*p + i]);
            oldK[lindx(i,i,p)] = K[i*p+i];
        }
        
        
        // check convergence
        if (maxdiff<*tol){
            converged = 1;
        }
        r++;
    }
    
    *maxit = r; // pass it back as number of iterations
    free(oldK);
}


void ScoreNN2(int *pIn, double *Xbar, double *S1, double *S2, double *S3, double *lambdaIn, double *K, double *b, double *tol, int *maxit){
  
  int i, j, k, r;
  int p = *pIn;
  double lambda = *lambdaIn;
  int tp = p*(p+1)/2;
  int converged=0;
  double sSum, s1, s2;
  double maxdiff=1.;

  double *oldK;
  oldK = (double*)malloc(tp*sizeof(double));
  if (oldK == 0){
    Rprintf("Out of Memory!\n");
    return;
  }
  
  double *oldb;
  oldb = (double*)malloc(p*sizeof(double));
  if (oldb == 0){
    Rprintf("Out of Memory!\n");
    return;
  }
  
  r = 0;
  
  while ((r<*maxit) && !converged) {
    maxdiff = 0;
    
    for (i = 0; i < (p-1); i++) {
      for (j = i+1; j < p; j++){
        s1 = 0;
        s2 = 0;
        for (k = 0; k < p; k++){
          s1 += K[i*p + k]*(S3[i*(p*p) + j*p + k]);
          s2 += K[j*p + k]*(S3[j*(p*p) + i*p + k]);
        }
        s1 = 2*S1[i*p + j] - s1 + S3[i*(p*p)+j*p+j]*K[i*p+j] + b[i]*S2[j*p+i];
        s2 = 2*S1[i*p + j] - s2 + S3[j*(p*p)+i*p+i]*K[j*p+i] + b[j]*S2[i*p+j];
        sSum = S3[i*(p*p) + j*p + j] + S3[j*(p*p) + i*p + i];
        K[i*p+j] = shrink((s1 + s2)/sSum, 2*lambda/sSum);
        K[j*p+i] = K[i*p+j];
        maxdiff += 2*fabs(oldK[lindx(i,j,p)] - K[i*p+j]);
        oldK[lindx(i,j,p)] = K[i*p+j];
      }
    }
    
    for (i = 0; i < p; i++){
      s1 = 0;
      sSum = S3[i*(p*p)+i*p+i];
      for (k = 0; k<p; k++){
        s1 += K[i*p+k]*S3[i*(p*p)+i*p+k];
      }
      s1 = 3*S1[i*p+i] - s1 + K[i*p+i] * S3[i*(p*p)+i*p+i] + b[i]*S2[i*p+i];
      K[i*p+i] = s1/sSum;
      maxdiff += fabs(oldK[lindx(i,i,p)] - K[i*p+i]);
      oldK[lindx(i,i,p)] = K[i*p+i];
    }
    
    for (i = 0; i < p; i++){
      s1 = 0;
      sSum = S1[i*p+i];
      for (k = 0; k<p; k++){
        s1 += K[i*p+k]*S2[k*p+i];
      }
      s1 += -2*Xbar[i];
      b[i] = s1/sSum;
      maxdiff += fabs(oldb[i] - b[i]);
      oldb[i] = b[i];
    }
    
    if (maxdiff<*tol){
      converged = 1;
    }
    r++;
  }
  *maxit = r;
  free(oldK);
  free(oldb);
}


// Conditional

void Score_cond(int *pIn, double *Xbar, double *S, double *SS, double *SD, double *SDQ, double *ST, double *lambdaIn, double *A, double *d, double *e, double *tol, int *maxit) {
    
    int i, j, k, r;
    int p = *pIn;
    double lambda = *lambdaIn; 
    int tp = p*(p+1)/2;
    int converged = 0;
    double sSum, s1, s2;
    double maxdiff = 1.;
    
    double *oldA;
    oldA = (double*)malloc(tp*sizeof(double));
    if (oldA == 0){
        Rprintf("Old of Memory!\n");
        return;
    }
    
    double *oldd;
    oldd = (double *)malloc(p*sizeof(double));
    if (oldd == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    double *olde;
    olde = (double *)malloc(p*sizeof(double));
    if (olde == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    r = 0;
    while ((r<*maxit) && !converged) {
        maxdiff = 0;
        
        for (i=0; i<(p-1); i++){
            for (j=i+1; j<p; j++){
                sSum = 16*(SDQ[i*p+j] + SDQ[j*p+i]);
                s1 = 0;
                s2 = 0;
                for (k = 0; k<p; k++){
                    s1 += 16 * ST[i*(p*p)+j*p+k]*A[i*p+k];
                    s2 += 16 * ST[j*(p*p)+i*p+k]*A[j*p+k];
                }
                s1 += 8*d[i]*SD[i*p+j] + 4*e[i]*SS[i*p+j] + 4*S[i*p+i] - 16*ST[i*(p*p)+j*p+j]*A[i*p+j];
                s2 += 8*d[j]*SD[j*p+i] + 4*e[j]*SS[j*p+i] + 4*S[j*p+j] - 16*ST[j*(p*p)+i*p+i]*A[j*p+i];
                A[i*p+j] = shrink(-(s1 + s2)/sSum, 2*lambda/sSum);
                A[j*p+i] = A[i*p+j];
                maxdiff += 2*fabs(oldA[lindx(i, j, p)] - A[i*p+j]);
                oldA[lindx(i, j, p)] = A[i*p+j];
            }
        }
        for (i=0; i<p; i++){
            sSum = 4*S[i*p+i];
            s1 = 0;
            for (k=0; k<p; k++){
                s1 += 8*SD[i*p+k]*A[i*p+k];
            }
            s1 += 2*e[i]*Xbar[i] + 2;
            d[i] = -s1/sSum;
            maxdiff += 2*fabs(oldd[i] - d[i]);
            oldd[i] = d[i];
        }
        for (i=0;i<p;i++){
            sSum = 1.0;
            s1 = 0;
            for (k=0; k<p; k++){
                s1 += 4*SS[i*p+k]*A[i*p+k];
            }
            s1 += 2*d[i]*Xbar[i];
            e[i] = -s1/sSum;
            maxdiff += 2*fabs(olde[i] - e[i]);
            olde[i] = e[i];
        }
        if (maxdiff<*tol){
            converged = 1;
        }
        r++;
    }
    *maxit = r; // pass it back as number of iterations
    free(oldA);
    free(oldd);
    free(olde);
}


// with x_ix_j cross terms

void Score_cond2(int *pIn, double *Xbar, double *S, double *SS, double *SSS, double *SD, double *SDQ, double *ST, double *lambdaIn, double *A, double *B, double *d, double *e, double *tol, int *maxit) {
    
    int i, j, k, r;
    int p = *pIn;
    double lambda = *lambdaIn;
    int tp = p*(p+1)/2;
    int converged = 0;
    double sSum, s1, s2;
    double maxdiff = 1.;
    
    double *oldA;
    oldA = (double*)malloc(tp*sizeof(double));
    if (oldA == 0){
        Rprintf("Old of Memory!\n");
        return;
    }
    
    double *oldB;
    oldB = (double*)malloc(tp*sizeof(double));
    if (oldB == 0){
        Rprintf("Old of Memory!\n");
        return;
    }
    
    double *oldd;
    oldd = (double *)malloc(p*sizeof(double));
    if (oldd == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    double *olde;
    olde = (double *)malloc(p*sizeof(double));
    if (olde == 0){
        Rprintf("Out of Memory!\n");
        return;
    }
    
    r = 0;
    while ((r<*maxit) && !converged) {
        maxdiff = 0;
        
        for (i=0; i<(p-1); i++){
            for (j=i+1; j<p; j++){
                sSum = 16*(SDQ[i*p+j] + SDQ[j*p+i]);
                s1 = 0;
                s2 = 0;
                for (k = 0; k<p; k++){
                    s1 += 16*ST[i*(p*p)+j*p+k] * A[i*p+k] + 8*SSS[j*(p*p)+i*p+k]*B[i*p+k];
                    s2 += 16*ST[j*(p*p)+i*p+k] * A[j*p+k] + 8*SSS[i*(p*p)+j*p+k]*B[j*p+k];
                }
                s1 += 8*d[i]*SD[i*p+j] + 4*e[i]*SS[i*p+j] + 4*S[i*p+i] - 16*ST[i*(p*p)+j*p+j]*A[i*p+j];
                s2 += 8*d[j]*SD[j*p+i] + 4*e[j]*SS[j*p+i] + 4*S[j*p+j] - 16*ST[j*(p*p)+i*p+i]*A[j*p+i];
                A[i*p+j] = shrink(-(s1 + s2)/sSum, 2*lambda/sSum);
                A[j*p+i] = A[i*p+j];
                maxdiff += 2*fabs(oldA[lindx(i, j, p)] - A[i*p+j]);
                oldA[lindx(i, j, p)] = A[i*p+j];
            }
        }
        for (i=0; i<(p-1); i++){
            for (j=i+1; j<p; j++){
                sSum = 4*(S[j*p+j] + S[i*p+i]);
                s1 = 0;
                s2 = 0;
                for (k = 0; k<p; k++){
                    s1 += 4*S[j*p+k] * B[i*p+k] + 8*SSS[k*(p*p)+i*p+j]*A[i*p+k];
                    s2 += 4*S[i*p+k] * B[j*p+k] + 8*SSS[k*(p*p)+j*p+i]*A[j*p+k];
                }
                s1 += 4*d[i]*S[i*p+j] + 2*e[i]*Xbar[j] - 4*S[j*p+j]*B[i*p+j];
                s2 += 4*d[j]*S[j*p+i] + 2*e[j]*Xbar[i] - 4*S[i*p+i]*B[j*p+i];
                B[i*p+j] = shrink(-(s1+s2)/sSum, 2*lambda/sSum);
                B[j*p+i] = B[i*p+j];
                maxdiff += 2*fabs(oldB[lindx(i,j,p)] - B[i*p+j]);
                oldB[lindx(i,j,p)] = B[i*p+j];
            }
        }
        for (i=0; i<p; i++){
            sSum = 4*S[i*p+i];
            s1 = 0;
            for (k=0; k<p; k++){
                s1 += 8*SD[i*p+k]*A[i*p+k] + 4*S[i*p+k]*B[i*p+k];
            }
            s1 += 2*e[i]*Xbar[i] + 2;
            d[i] = -s1/sSum;
            maxdiff += 2*fabs(oldd[i] - d[i]);
            oldd[i] = d[i];
        }
        for (i=0;i<p;i++){
            sSum = 1.0;
            s1 = 0;
            for (k=0; k<p; k++){
                s1 += 4*SS[i*p+k]*A[i*p+k] + 2*Xbar[k]*B[i*p+k];
            }
            s1 += 2*d[i]*Xbar[i];
            e[i] = -s1/sSum;
            maxdiff += 2*fabs(olde[i] - e[i]);
            olde[i] = e[i];
        }
        if (maxdiff<*tol){
            converged = 1;
        }
        r++;
    }
    *maxit = r; // pass it back as number of iterations
    free(oldA);
    free(oldB);
    free(oldd);
    free(olde);
}


