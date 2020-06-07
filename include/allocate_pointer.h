#pragma once
#include <stdio.h>
#include <stdlib.h>
///////malloc memory
int ****D4i(int r, int c, int p, int s); // malloc 2d array of Double
int ***D3i(int r, int c, int p);
double **D2d(int r, int c); // malloc 2d array of Double
int **D2i(int r, int c); // malloc 2d array of Integer
int *D1i(int x);		 // malloc 1d array of Integer
double *D1d(int x); // malloc 1d array of Double
