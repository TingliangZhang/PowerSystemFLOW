#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning( disable : 4996)
#define f1(i) (i-1)  
#define f2(i,j,n) ((i-1)*(n)+j-1)
int ReadLineParameter(float* lineZ_REAL, float* lineZ_IMAG, float* lineY_2_IMAG);
int ReadLineNumber(int line_number, int* line_START, int* line_END);
int ReadNodeProperty(int* node_totalnumber, int* node_PQnumber, float*node_P, float* node_Q, float* node_U, float* node_D, float* node_B_ground, float* node_Qmin, float* node_Qmax);
int ReadUnbalance_B(int node_number, int line_number, float* node_B_unbalance);
void ybus(int n, int l, int m, float *g, float *b, float *g1, float *b1, float *c1, float *c, float *co, int k, int *s1, int *e1);
void dpqc(float *p, float *q, float *p0, float *q0, float *v, float *v0, int m, int n, float *e, float *f, int k, float *g, float *b, float *dd, int *node_PVtoPQ);
void jmcc(int m, int n, int n0, float *e, float *f, float *g, float *b, float *jm, int k, int * node_PVtoPQ);
void sevc(float a[], int n0, int k, int n1);
void plsc(int n, int l, int m, float g[], float b[], float e[], float f[], int e1[], int s1[], float g1[], float b1[], float c1[], float c[], \
	float co[], float p1[], float q1[], float p2[], float q2[], float p3[], float q3[], float p[], float q[], float v[], float angle[], int k1, int* node_PVtoPQ);