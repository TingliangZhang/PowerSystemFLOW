#include<stdio.h>
#include "stdlib.h"
#include "string.h"
#pragma warning( disable : 4996)
#define f1(i) (i-1)  
int ReadNodeProperty(int* node_totalnumber, int* node_PQnumber, float*node_P, float* node_Q, float* node_U, float* node_D, float* node_B_ground, float* node_Qmin, float* node_Qmax) {
	FILE *fp = NULL;
	int i = 0, temp = 0;
	char node_property[10] = { 0 };
	int cnt = 0;
	int cnt_PQ = 0;
	int cnt_PV = 0;
	int cnt_balance = 0;
	int pos1 = 0;

	if ((node_totalnumber == NULL) || (node_PQnumber == NULL) || (node_P == NULL) || (node_Q == NULL) || (node_U == NULL) || (node_D == NULL) || (node_B_ground == NULL) || (node_Qmax == NULL) || (node_Qmin == NULL))
	{
		printf("No enough memory !\n");
		exit(-1);
	}
	if ((fp = fopen("data_nodeprop.txt", "r")) == NULL)
	{
		printf("The file is not opened\n");
		exit(0);
	}
	else
	{

		for (i = 0; feof(fp) == 0; i++)
		{
			fscanf(fp, "%d %s", &temp, node_property);
			if (strcmp(node_property, "PQ") == 0)//PQ node
			{
				cnt++; cnt_PQ++;
				fscanf(fp, "%f %f %f", node_P + i, node_Q + i, node_B_ground + i);
				*(node_U + i) = 1;
				*(node_D + i) = 0;
			}
			else if (strcmp(node_property, "PV") == 0)//PV node
			{
				cnt++; cnt_PV++;
				fscanf(fp, "%f %f %f %f %f", node_P + i, node_U + i, node_B_ground + i, node_Qmin + i, node_Qmax + i);
				*(node_Q + i) = 0;
				*(node_D + i) = 0;
			}
			else if (strcmp(node_property, "BALANCE") == 0)//balance node
			{
				cnt++; cnt_balance++;
				fscanf(fp, "%f %f %f", node_U + i, node_D + i, node_B_ground + i);
				*(node_P + i) = 0;
				*(node_Q + i) = 0;
			}
			else
				return 1;
		}

		fclose(fp);
		fp = NULL;
	}
	if (cnt_balance != 1)
		return 1;
	else
	{
		*node_totalnumber = cnt;
		*node_PQnumber = cnt_PQ;
		return 0;
	}
}