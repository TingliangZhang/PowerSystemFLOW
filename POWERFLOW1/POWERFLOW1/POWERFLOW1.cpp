// POWERFLOW1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "function.h"
#define isPrint 0
#define Threshold 1e-5
#define startPVtoPQmonitor 1e-1
#define Initial_E_value 1
#define Initial_F_value 0
#define isPOLAR 0
#define MAX_iteration_number 1e4
#define MAX_ARRAY_NUMBER 100
#define PI 3.14159265
int main() {
	int test = 0;
	int line_number = 0;
	int node_totalnumber = 0;
	int node_PQnumber = 0;
	int i = 0;
	int j = 0;
	int pos1 = 0;
	int pos2 = 0;
	int pos22 = 0;
	int gamma = 0;
	float delta = 1;
	float* lineZ_REAL = NULL; lineZ_REAL = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* lineZ_IMAG = NULL; lineZ_IMAG = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* lineY_2_IMAG = NULL; lineY_2_IMAG = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_B_ground = NULL; node_B_ground = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	int* line_START = NULL; line_START = (int*)malloc(sizeof(int) * MAX_ARRAY_NUMBER);
	int* line_END = NULL; line_END = (int*)malloc(sizeof(int) * MAX_ARRAY_NUMBER);
	float* node_P = NULL; node_P = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_Q = NULL; node_Q = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_U = NULL; node_U = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_D = NULL; node_D = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_Qmax = NULL; node_Qmax = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_Qmin = NULL; node_Qmin = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_Pt = NULL; node_Pt = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_Qt = NULL; node_Qt = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	float* node_Ut = NULL; node_Ut = (float*)malloc(sizeof(float) * MAX_ARRAY_NUMBER);
	/*动态分配内存*/
	for (i = 1; i <= MAX_ARRAY_NUMBER; i++)
	{
		pos1 = f1(i);
		lineZ_REAL[pos1] = 0; lineZ_IMAG[pos1] = 0;
		line_START[pos1] = 0; line_END[pos1] = 0;
		lineY_2_IMAG[pos1] = 0;
		node_B_ground[pos1] = 0;
		node_P[pos1] = 0; node_Q[pos1] = 0;
		node_U[pos1] = 1; node_D[pos1] = 0;
		node_Qmax[pos1] = 999; node_Qmin[pos1] = -999;
		node_Pt[pos1] = 0; node_Qt[pos1] = 0;
		node_Ut[pos1] = 1;
	}
	/*读参数*/

	//读线路参数
	if ((line_number = ReadLineParameter(lineZ_REAL, lineZ_IMAG, lineY_2_IMAG)) == -1)
		exit(-1);
	//读节点连接
	if ((test = ReadLineNumber(line_number, line_START, line_END)) == -1)
		exit(-1);
	//读节点性质
	if ((test = ReadNodeProperty(&node_totalnumber, &node_PQnumber, node_P, node_Q, node_U, node_D, node_B_ground, node_Qmin, node_Qmax)) != 0)
		exit(-1);
	float* node_B_unbalance = NULL; node_B_unbalance = (float*)malloc(sizeof(float) * node_totalnumber * line_number);
	//读不对称节点
	if ((test = ReadUnbalance_B(node_totalnumber, line_number, node_B_unbalance)) != 0)
		exit(-1);

	float* node_G = NULL; node_G = (float*)malloc(sizeof(float) * (node_totalnumber * node_totalnumber));
	float* node_B = NULL; node_B = (float*)malloc(sizeof(float) * (node_totalnumber * node_totalnumber));
	float* node_E = NULL; node_E = (float*)malloc(sizeof(float) * (node_totalnumber));
	float* node_F = NULL; node_F = (float*)malloc(sizeof(float) * (node_totalnumber));
	float* node_deltaP = NULL; node_deltaP = (float*)malloc(sizeof(float) * (node_totalnumber));
	float* node_deltaQ = NULL; node_deltaQ = (float*)malloc(sizeof(float) * (node_totalnumber));
	float* node_deltaU = NULL; node_deltaU = (float*)malloc(sizeof(float) * (node_totalnumber));
	int Jacobi_dim = node_totalnumber * 2 - 2;
	float* Jacobi = NULL; Jacobi = (float*)malloc(sizeof(float) * (Jacobi_dim) * (Jacobi_dim));
	float* aug_matrix = NULL; aug_matrix = (float*)malloc(sizeof(float) * (Jacobi_dim) * (Jacobi_dim + 1));
	float* node_P1 = NULL; node_P1 = (float*)malloc(sizeof(float) * line_number);
	float* node_Q1 = NULL; node_Q1 = (float*)malloc(sizeof(float) * line_number);
	float* node_P2 = NULL; node_P2 = (float*)malloc(sizeof(float) * line_number);
	float* node_Q2 = NULL; node_Q2 = (float*)malloc(sizeof(float) * line_number);
	float* node_P3 = NULL; node_P3 = (float*)malloc(sizeof(float) * line_number);
	float* node_Q3 = NULL; node_Q3 = (float*)malloc(sizeof(float) * line_number);
	int* node_PVtoPQ = NULL; node_PVtoPQ = (int*)malloc(sizeof(int) * line_number);
	if ((node_G == NULL) || (node_B == NULL) || (node_E == NULL) || (node_F == NULL)
		|| (node_deltaP == NULL) || (node_deltaQ == NULL) || (node_deltaU == NULL)
		|| (node_PVtoPQ == NULL) || (node_P1 == NULL) || (node_P2 == NULL)
		|| (node_P3 == NULL) || (node_Q1 == NULL) || (node_Q2 == NULL)
		|| (node_Q3 == NULL) || (Jacobi == NULL) || (aug_matrix == NULL))
	{
		printf("No enough memory !\n");
		exit(-1);
	}
	/*节点P1,P2,P3,Q1,Q2,Q3初始化*/
	for (i = 1; i <= line_number; i++)
	{
		pos1 = f1(i);
		node_P1[pos1] = 0; node_Q1[pos1] = 0;
		node_P2[pos1] = 0; node_Q2[pos1] = 0;
		node_P3[pos1] = 0; node_Q3[pos1] = 0;
		node_PVtoPQ[pos1] = 0;
	}
	/*雅可比矩阵初始化*/
	for (i = 1; i <= Jacobi_dim; i++)
	{
		for (j = 1; j <= Jacobi_dim; j++)
		{
			pos2 = f2(i, j, Jacobi_dim);
			Jacobi[pos2] = 0;
		}
	}
	/*note_E,node_F 初始化*/
	for (i = 1; i <= node_totalnumber; i++)
	{
		pos1 = f1(i); node_E[pos1] = Initial_E_value; node_F[pos1] = Initial_F_value;
		node_deltaP[pos1] = 0; node_deltaQ[pos1] = 0; node_deltaU[pos1] = 0;
	}
	pos1 = f1(node_totalnumber);
	node_E[pos1] = node_U[pos1] * cos(node_D[pos1] * PI / 180);
	node_F[pos1] = node_U[pos1] * sin(node_D[pos1] * PI / 180);
	/*计算导纳矩阵*/
	ybus(node_totalnumber, line_number, node_PQnumber, node_G, node_B,
		lineZ_REAL, lineZ_IMAG, lineY_2_IMAG, node_B_unbalance, node_B_ground,
		isPrint, line_START, line_END);
	while ((fabs(delta) > Threshold) && (gamma < MAX_iteration_number))
	{

		if (fabs(delta) < startPVtoPQmonitor)//达到一定误差范围后，开始进行PV转PQ检测
		{
			plsc(node_totalnumber, line_number, node_PQnumber,
				node_G, node_B, node_E, node_F, line_END, line_START,
				lineY_2_IMAG, lineZ_IMAG, lineY_2_IMAG, node_B_unbalance,
				node_B_ground, node_P1, node_Q1, node_P2, node_Q2, node_P3,
				node_Q3, node_Pt, node_Qt, node_Ut, node_D, isPOLAR, node_PVtoPQ);
			for (i = 1; i < (node_totalnumber - node_PQnumber); i++)
			{
				pos1 = f1(node_PQnumber + i);
				if (node_Qt[pos1] > node_Qmax[pos1])
				{
					node_PVtoPQ[pos1] = 1;
					node_deltaQ[pos1] = node_Qmax[pos1] - node_Qt[pos1];
					node_deltaU[pos1] = 0;
					node_Q[pos1] = node_Qmax[pos1];
					node_Qt[pos1] = node_Qmax[pos1];
					gamma = 0;
				}
				else if (node_Qt[pos1] < node_Qmin[pos1])
				{
					node_PVtoPQ[pos1] = -1;
					node_deltaQ[pos1] = node_Qmin[pos1] - node_Qt[pos1];
					node_deltaU[pos1] = 0;
					node_Q[pos1] = node_Qmin[pos1];
					node_Qt[pos1] = node_Qmin[pos1];
					gamma = 0;
				}
				else
					;
			}
		}

	AA://标识符
		/*计算deltaP和deltaQ*/
		dpqc(node_P, node_Q, node_deltaP, node_deltaQ, node_U, node_deltaU,
			node_PQnumber, node_totalnumber, node_E, node_F, isPrint,
			node_G, node_B, &delta, node_PVtoPQ);
		/*计算雅可比矩阵*/
		jmcc(node_PQnumber, node_totalnumber, Jacobi_dim, node_E, node_F,
			node_G, node_B, Jacobi, isPrint, node_PVtoPQ);
		/*构造增广矩阵*/
		{
			for (i = 1; i <= (node_PQnumber * 2); i++)
			{
				for (j = 1; j <= Jacobi_dim; j++)
				{
					pos2 = f2(i, j, Jacobi_dim);
					pos22 = f2(i, j, Jacobi_dim + 1);
					aug_matrix[pos22] = Jacobi[pos2];
				}
				pos1 = f1((i + 1) / 2);
				pos22 = f2(i, j, Jacobi_dim + 1);
				aug_matrix[pos22] = node_deltaQ[pos1];
				i++;
				for (j = 1; j <= Jacobi_dim; j++)
				{
					pos2 = f2(i, j, Jacobi_dim);
					pos22 = f2(i, j, Jacobi_dim + 1);
					aug_matrix[pos22] = Jacobi[pos2];
				}
				pos1 = f1((i + 1) / 2); pos22 = f2(i, j, Jacobi_dim + 1);
				aug_matrix[pos22] = node_deltaP[pos1];
			}
			for (i = 2 * node_PQnumber + 1; i <= Jacobi_dim; i++)
			{
				for (j = 1; j <= Jacobi_dim; j++)
				{
					pos2 = f2(i, j, Jacobi_dim);
					pos22 = f2(i, j, Jacobi_dim + 1);
					aug_matrix[pos22] = Jacobi[pos2];
				}
				pos1 = f1((i + 1) / 2); pos22 = f2(i, j, Jacobi_dim + 1);
				if (node_PVtoPQ[pos1] == 0)
					aug_matrix[pos22] = node_deltaU[pos1];
				else
					aug_matrix[pos22] = node_deltaQ[pos1];
				i++;
				for (j = 1; j <= Jacobi_dim; j++)
				{
					pos2 = f2(i, j, Jacobi_dim);
					pos22 = f2(i, j, Jacobi_dim + 1);
					aug_matrix[pos22] = Jacobi[pos2];
				}
				pos1 = f1((i + 1) / 2); pos22 = f2(i, j, Jacobi_dim + 1);
				aug_matrix[pos22] = node_deltaP[pos1];
			}
		}
		/*求解线性方程组*/
		sevc(aug_matrix, Jacobi_dim, isPrint, Jacobi_dim + 1);
		/*回代作差*/
		for (i = 1; i <= Jacobi_dim; i++)
		{
			pos1 = f1((i + 1) / 2);
			pos22 = f2(i, (Jacobi_dim + 1), (Jacobi_dim + 1));
			if ((i % 2) == 1)
			{
				node_E[pos1] = node_E[pos1] - aug_matrix[pos22];
			}
			else
			{
				node_F[pos1] = node_F[pos1] - aug_matrix[pos22];
			}
		}

		gamma++;


	}
	/*求解功率分布*/
	plsc(node_totalnumber, line_number, node_PQnumber, node_G, node_B, node_E,
		node_F, line_END, line_START, lineY_2_IMAG, lineZ_IMAG, lineY_2_IMAG,
		node_B_unbalance, node_B_ground, node_P1, node_Q1, node_P2, node_Q2,
		node_P3, node_Q3, node_Pt, node_Qt, node_Ut, node_D, isPOLAR, node_PVtoPQ);
	/*检查PV结点是否越界*/

	for (i = 1; i <= (node_totalnumber - node_PQnumber); i++)
	{
		pos1 = f1(node_PQnumber + i);
		if (node_Qt[pos1] > node_Qmax[pos1])
		{
			node_PVtoPQ[pos1] = 1;
			node_deltaQ[pos1] = node_Qmax[pos1] - node_Qt[pos1];
			node_Q[pos1] = node_Qmax[pos1];
			node_Qt[pos1] = node_Qmax[pos1];
			gamma = 0;
			goto AA;
		}
		else if (node_Qt[pos1] < node_Qmin[pos1])
		{
			node_PVtoPQ[pos1] = -1;
			node_deltaQ[pos1] = node_Qmin[pos1] - node_Qt[pos1];
			node_Q[pos1] = node_Qmin[pos1];
			node_Qt[pos1] = node_Qmin[pos1];
			gamma = 0;
			goto AA;
		}
		else
			;
	}
	plsc(node_totalnumber, line_number, node_PQnumber,
		node_G, node_B, node_E, node_F, line_END, line_START,
		lineY_2_IMAG, lineZ_IMAG, lineY_2_IMAG, node_B_unbalance,
		node_B_ground, node_P1, node_Q1, node_P2, node_Q2, node_P3,
		node_Q3, node_P, node_Q, node_U, node_D, isPOLAR, node_PVtoPQ);

	if (delta <= Threshold)
		printf("SUCCESS!\nPower system flow calculation has been successfully completed!\n\n");
	else
		printf("ATTETION!\nThe maximum number of iterations has been reached and the precision target is not reached.\n\n");

	printf("The set threshold is %e.\n", Threshold);
	printf("The maximum error is %e, %d iterations are performed.\nThe result is saved at \"output.txt\".\n\n", delta, gamma);
	/*释放空间*/
	free(node_P); free(node_Q); free(node_U); free(node_D); free(line_START);
	free(line_END); free(lineZ_REAL); free(lineZ_IMAG); free(lineY_2_IMAG);
	free(node_B_ground); free(node_B_unbalance); free(node_Qmax); free(node_Qmin);
	free(node_G); free(node_B); free(node_E); free(node_F); free(node_deltaP);
	free(node_deltaQ); free(node_deltaU); free(node_PVtoPQ); free(node_P1);
	free(node_P2); free(node_P3); free(node_Q1); free(node_Q2); free(node_Q3);
	free(Jacobi); free(aug_matrix); free(node_Pt); free(node_Qt); free(node_Ut);
	system("pause");
	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
