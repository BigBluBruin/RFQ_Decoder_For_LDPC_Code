///////malloc memory

#include "allocate_pointer.h"
int**** D4i(int r, int c, int p, int s)   // malloc 2d array of Double
{
	int**** theArray;           
	int index1, index2, index3;
	theArray = (int****)malloc(r * sizeof(int***));
	for (index1 = 0; index1 < r; index1++)
	{
		theArray[index1] = (int ***)malloc(c * sizeof(int**));
		for (index2 = 0; index2 < c; index2++)
		{
			theArray[index1][index2] = (int**)malloc(p * sizeof(int*));
			for (index3 = 0; index3 < p; index3++)
			{
				theArray[index1][index2][index3] = (int*)malloc(s * sizeof(int));
			}
		}
	}
	return theArray;
}
int*** D3i(int r, int c, int p)
{
	int*** theArray;
	int index1, index2;
	theArray = (int***)malloc(r * sizeof(int**));
	for (index1 = 0;index1<r;index1++)
	{
		theArray[index1] = (int **)malloc(c * sizeof(int*));
		for (index2 = 0;index2<c;index2++)
		{
			theArray[index1][index2] = (int*)malloc(p * sizeof(int));
		}
	}
	return theArray;

}
double** D2d(int r, int c)   // malloc 2d array of Double
{
	double** theArray;
	int i;
	theArray = (double**)malloc(r * sizeof(double*));
	for (i = 0; i < r; i++)
		theArray[i] = (double*)malloc(c * sizeof(double));
	return theArray;
}

int** D2i(int r, int c)   // malloc 2d array of Integer
{
	int** theArray;
	int i;
	theArray = (int**)malloc(r * sizeof(int*));
	for (i = 0; i < r; i++)
		theArray[i] = (int*)malloc(c * sizeof(int));    
	return theArray;
}

int* D1i(int x)   // malloc 1d array of Integer
{
	int* theArray;
	theArray = (int*)malloc(x * sizeof(int));  //array pointer
	return theArray;
}

double* D1d(int x)    // malloc 1d array of Double
{
	double* theArray;
	theArray = (double*)malloc(x * sizeof(double)); //array pointer
	return theArray;
}