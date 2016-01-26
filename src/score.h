//
//  score.h
//  
//
//  Created by Lina on 2015-04-25.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif


double shrink(double a, double b);

int lindx(int r, int c, int p);

void Score1(int *pIn, double *S, double *lambdaIn, double *K, double *tol, int *maxit);

void Score2(int *nIn, int *pIn, double *Y, double *S, double *lambdaIn, double *K,
            double *tol, int *maxit);

void ScoreNN1(int *pIn,  double *S1, double *S2, double *lambdaIn, double *K, double *tol, int *maxit);

void ScoreNN2(int *pIn, double *Xbar, double *S1, double *S2, double *S3, double *lambdaIn, double *K, double *b, double *tol, int *maxit);

void Score_cond(int *pIn, double *Xbar, double *S, double *SS, double *SD, double *SDQ, double *ST, double *lambdaIn, double *A, double *d, double *e, double *tol, int *maxit);

void Score_cond2(int *pIn, double *Xbar, double *S, double *SS, double *SSS, double *SD, double *SDQ, double *ST, double *lambdaIn, double *A, double *B, double *d, double *e, double *tol, int *maxit);