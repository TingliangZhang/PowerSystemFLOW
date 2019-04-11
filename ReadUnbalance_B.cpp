#include<stdio.h>
#include "stdlib.h"
#pragma warning( disable : 4996)
#define f1(i) (i-1)  
/* 把习惯的一阶矩阵的下标转化为C语言数组下标*/

#define f2(i,j,n) ((i-1)*(n)+j-1)
/* 把习惯的二阶矩阵的下标转化为C语言数组下标*/
int ReadUnbalance_B(int node_number, int line_number, float* node_B_unbalance) {
	FILE *fp = NULL;
	int i = 0, j = 0;
	int pos2 = 0;
	int node = 0;
	int line = 0;
	float value = 0;
	if (node_B_unbalance == NULL)
	{
		printf("No enough memory !\n");
		exit(-1);
	}
	if ((fp = fopen("data_unbalance.txt", "r")) == NULL)
	{
		printf("The file is not opened\n");
		exit(0);
	}
	else
	{
		for (i = 1; i <= node_number; i++)
		{
			for (j = 1; j <= line_number; j++)
			{
				pos2 = f2(i, j, line_number);
				node_B_unbalance[pos2] = 0;
			}
		}
		for (i = 0; feof(fp) == 0; i++)
		{
			fscanf(fp, "%d %d %f", &node, &line, &value);
			pos2 = f2(node, line, line_number);
			node_B_unbalance[pos2] = value;
		}
		fclose(fp);
		fp = NULL;
		return 0;
	}
}