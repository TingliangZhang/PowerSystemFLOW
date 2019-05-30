#include<stdio.h>
#include "stdlib.h"
#pragma warning( disable : 4996)
int ReadLineParameter(float* lineZ_REAL,float* lineZ_IMAG,float* lineY_2_IMAG) {
	FILE *fp = NULL;
	char sign = '+';
	int i;

	int num_line = 0;

	if ((lineZ_REAL == NULL) || (lineZ_IMAG == NULL)|| (lineY_2_IMAG == NULL))
	{
		printf("No enough memory !\n");
		exit(-1);
	}

	if ((fp = fopen("data_line.txt", "r")) == NULL)
	{
		printf("The file is not opened\n");
		exit(-1);
	}
	else
	{
		for (i = 0; feof(fp) == 0; i++)
		{
			fscanf(fp, "%d", &num_line);
			fscanf(fp, "%f%cj%f", lineZ_REAL + i, &sign, lineZ_IMAG + i);
			if (sign == '-')
				*(lineZ_IMAG+i) = *(lineZ_IMAG+i)*(-1);
			fscanf(fp, " %cj%f", &sign, lineY_2_IMAG + i);
			if (sign == '-')
				*(lineY_2_IMAG+i) = *(lineY_2_IMAG+i)*(-1);

		}
		num_line = i;
		fclose(fp);
		fp = NULL;
	}

	return num_line;
}
