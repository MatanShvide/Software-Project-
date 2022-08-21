#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef SPKMEANS_H
#define SPKMEANS_H

void print2Darray(double**, int, int); /* delete later !!!!!!!!!!!!!!!!!!!!! */
double** mem2D(int, int);
void free2D(double **);
double **mainFuncCapi(double **, int , int , int );
double** sortMat(double**);
int eigenGap(double**);
double **kmeans(double **tMat, int *initCentroids, int, int);

#endif
