#include<stdio.h>
#include "stdlib.h"
#pragma warning( disable : 4996)
#define f1(i) (i-1)  
int ReadLineNumber(int line_number, int* line_START, int* line_END) {
	FILE *fp = NULL;
	int i = 0, temp = 0;

	if ((line_START == NULL) || (line_END == NULL))
	{
		printf("No enough memory !\n");
		exit(-1);
	}
	if ((fp = fopen("data_connect.txt", "r")) == NULL)
	{
		printf("The file is not opened\n");
		exit(-1);
	}
	else
	{
		for (i = 0; i < line_number; i++)
		{
			fscanf(fp, "%d %d %d", &temp, line_START + i, line_END + i);
		}
		fclose(fp);
		fp = NULL;
	}

	return 0;
}