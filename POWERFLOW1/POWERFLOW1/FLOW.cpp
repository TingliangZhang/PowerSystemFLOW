/**************************FLOW.C*********************************/

/*******************************************************************
 *     这里提供的是电力系统潮流计算机解法的五个子程序,采用的方法是 *
 *  Newton_Raphson法.                                              *
 *     程序中所用的变量说明如下:                                   *
 *       N:网络节点总数.         M:网络的PQ节点数.                 *
 *       L:网络的支路总数.       N0:雅可比矩阵的行数.              *
 *       N1:N0+1                 K:打印开关.K=1,则打印;否则,不打印.*
 *       K1:子程序PLSC中判断输入电压的形式.K1=1,则为极座标形式.否则*
 *          为直角坐标形式.                                        *
 *       D:有功及无功功率误差的最大值（绝对值）.                   *
 *       G(I,J):Ybus的电导元素(实部).                              *
 *       B(I,J):Ybus的电纳元素(虚部).                              *
 *       G1(I) :第I支路的串联电导.      B1(I):第I支路的串联电纳.   *
 *       C1(I) :第I支路的pie型对称接地电纳.                        *
 *       C(I,J):第I节点J支路不对称接地电纳.                        *
 *       CO(I) :第I节点的接地电纳.                                 *
 *       S1(I) :第I支路的起始节点号.    E1(I):第I支路的终止节点号. *
 *       P(I)  :第I节点的注入有功功率.  Q(I):第I节点的注入无功功率.*
 *       P0(I) :第I节点有功功率误差.    Q0(I):第I节点无功功率误差. *
 *       V0(I) :第I节点(PV节点)的电压误差(平方误差).               *
 *       V(I)  :第I节点的电压幅值.                                 *
 *       E(I)  :第I节点的电压的实部.    F(I):第I节点的电压的虚部.  *
 *      JM(I,J):Jacoby矩阵的第I行J列元素.                          *
 *       A(I,J):修正方程的增广矩阵,三角化矩阵的第I行J列元素,运算结 *
 *              束后A矩阵的最后一列存放修正的解.                   *
 *       P1(I) :第I支路由S1(I)节点注入的有功功率.                  *
 *       Q1(I) :第I支路由S1(I)节点注入的无功功率.                  *
 *       P2(I) :第I支路由E1(I)节点注入的有功功率.                  *
 *       Q2(I) :第I支路由E1(I)节点注入的无功功率.                  *
 *       P3(I) :第I支路的有功功率损耗.                             *
 *       Q3(I) :第I支路的无功功率损耗.                             *
 *     ANGLE(I):第I节点电压的角度.                                 *
 *     节点编号顺序：N个节点中，前M个为PQ节点，M+1至N-1个节点为    *
 *                   PV节点，第N个为平衡节点                       *
*******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#pragma warning( disable : 4996)
#define f1(i) (i-1)  
/* 把习惯的一阶矩阵的下标转化为C语言数组下标*/

#define f2(i,j,n) ((i-1)*(n)+j-1)
/* 把习惯的二阶矩阵的下标转化为C语言数组下标*/

/****************************************************
 *    本子程序根据所给的支路导纳及有关信息,形成结点 *
 * 导纳矩阵,如打印参数K=1,则输出电导矩阵G和电纳矩B  *
 ****************************************************/
void ybus(int n, int l, int m, float *g, float *b, float *g1, float *b1, float *c1, \
	float *c, float *co, int k, int *s1, int *e1)
{
	FILE *fp;
	int i, j, io, i0;
	int pos1, pos2;
	int st, en;
	if ((fp = fopen("output.txt", "w+")) == NULL)
		return;
	else
	{

		/* 初始化矩阵G,B */
		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= n; j++)
			{
				pos2 = f2(i, j, n);
				g[pos2] = 0; b[pos2] = 0;
			}
		}

		/* 计算支路导纳 */
		for (i = 1; i <= l; i++)
		{
			/* 计算对角元 */
			pos1 = f1(i);
			st = s1[pos1]; en = e1[pos1];
			pos2 = f2(st, st, n);
			g[pos2] += g1[pos1];
			b[pos2] += b1[pos1] + c1[pos1];
			pos2 = f2(en, en, n);
			g[pos2] += g1[pos1];
			b[pos2] += b1[pos1] + c1[pos1];

			/* 计算非对角元 */
			pos2 = f2(st, en, n);
			g[pos2] -= g1[pos1];
			b[pos2] -= b1[pos1];
			g[f2(en, st, n)] = g[f2(st, en, n)];
			b[f2(en, st, n)] = b[f2(st, en, n)];
		}

		/* 计算接地支路导纳 */
		for (i = 1; i <= n; i++)
		{
			/*　对称部分　*/
			b[f2(i, i, n)] += co[f1(i)];

			/* 非对称部分　*/
			for (j = 1; j <= l; j++)
			{
				b[f2(i, i, n)] += c[f2(i, j, l)];
			}
		}
		if (k != 1)
		{
			fclose(fp);
			fp = NULL;
			return; /* 如果K不为 1,则返回;否则,打印导纳矩阵 */
		}
		fprintf(fp, "\n          BUS ADMITTANCE MATRIX Y(BUS):");
		fprintf(fp, "\n ******************* ARRAY G ********************");
		for (io = 1; io <= n; io += 5)
		{
			i0 = (io + 4) > n ? n : (io + 4);
			fprintf(fp, "\n");
			for (j = io; j <= i0; j++)
			{
				fprintf(fp, "%13d", j);
			}
			for (i = 1; i <= n; i++)
			{
				fprintf(fp, "\n%2d", i);
				for (j = io; j <= i0; j++)
				{
					fprintf(fp, "%13.6f", g[f2(i, j, n)]);
				}
			}
			fprintf(fp, "\n");
		}

		fprintf(fp, "\n ******************* ARRAY B ********************");
		for (io = 1; io <= n; io += 5)
		{
			i0 = (io + 4) > n ? n : (io + 4);
			fprintf(fp, "\n");
			for (j = io; j <= i0; j++)
			{
				fprintf(fp, "%13d", j);
			}
			for (i = 1; i <= n; i++)
			{
				fprintf(fp, "\n%2d", i);
				for (j = io; j <= i0; j++)
				{
					fprintf(fp, "%13.6f", b[f2(i, j, n)]);
				}
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n************************************************");
		fclose(fp);
		fp = NULL;
	}
}

/*******************************************
 *     本子程序根据所给的功率及电压等数据  *
 * 求出功率及电压误差量,并返回最大有功功率 *
 * 以用于与给定误差比较.如打印参数K=1,则输 *
 * 出P0,Q0(对PQ结点),V0(对PV结点).         *
 * 对应附录一P177式(4-86)(4-87)			   *
 *******************************************/

void dpqc(float *p, float *q, float *p0, float *q0, float *v, float *v0, int m, \
	int n, float *e, float *f, int k, float *g, float *b, float *dd, int *node_PVtoPQ)
{
	FILE *fp;
	int i, j, l;
	int pos1, pos2;
	float a1, a2, d1, d;
	if ((fp = fopen("output.txt", "a+")) == NULL)
		return;
	else
	{
		l = n - 1;
		if (k == 1)
		{
			fprintf(fp, "\n        CHANGE OF P0,V**2,P0(I),Q0(I),V0(I) ");
			fprintf(fp, "\n  I       P0(I)           Q0(I)");
		}
		for (i = 1; i <= l; i++)
		{
			a1 = 0; a2 = 0;
			pos1 = f1(i);
			for (j = 1; j <= n; j++)
			{
				/* a1,a2对应附录一P177式(4-86)中括号内的式子 */
				pos2 = f2(i, j, n);
				a1 += g[pos2] * e[f1(j)] - b[pos2] * f[f1(j)];
				a2 += g[pos2] * f[f1(j)] + b[pos2] * e[f1(j)];
			}
			/* 计算式(4-86)(4-87)中的deltaPi　*/
			p0[pos1] = p[pos1] - e[pos1] * a1 - f[pos1] * a2;
			if (i <= m)
			{	/* 计算PQ结点中的deltaQi　*/
				q0[pos1] = q[pos1] - f[pos1] * a1 + e[pos1] * a2;
			}
			else if ((node_PVtoPQ[pos1] == 1) || (node_PVtoPQ[pos1] == -1))
			{	/*计算PVtoPQ节点的deltaQi*/
				q0[pos1] = q[pos1] - f[pos1] * a1 + e[pos1] * a2;
			}
			else
			{	/* 计算PV结点中的deltaVi平方　*/
				v0[pos1] = v[pos1] * v[pos1] - e[pos1] * e[pos1] - f[pos1] * f[pos1];
			}

			/* 输出结果 */
			if (k == 1)
			{
				if (i < m)
				{
					fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], q0[pos1]);
				}
				else if (i == m)
				{
					fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], q0[pos1]);
					fprintf(fp, "\n  I       P0(I)           V0(I)");
				}
				else if (node_PVtoPQ[pos1] != 0)
				{
					fprintf(fp, "\n  I       P0(I)           Q0(I)");
					fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], q0[pos1]);
				}
				else
				{
					fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], v0[pos1]);
				}
			}
		}

		/* 找到deltaP和deltaQ中的最大者，作为收敛指标, 存在dd中 */
		d = 0;
		for (i = 1; i <= l; i++)
		{
			pos1 = f1(i);
			d1 = p0[pos1] > 0 ? p0[pos1] : -p0[pos1];
			if (d < d1)
			{
				d = d1;
			}
			if ((i <= m) || (node_PVtoPQ[pos1] != 0))
			{
				d1 = q0[pos1] > 0 ? q0[pos1] : -q0[pos1];
				if (d < d1)
				{
					d = d1;
				}
			}
		}
		(*dd) = d;
		fclose(fp);
		fp = NULL;
	}
}

/***************************************************
 *    本子程序根据节点导纳及电压求Jacoby矩阵,用于求*
 *  电压修正量,如打印参数K=1,则输出Jacoby矩阵.     *
 *  对应于附录一P178式(4-89)(4-90)				   *
 *    值得注意的是，程序中Jacobi阵中H N J L的排列顺*
 *  序与式（4－88）略有不同，程序中H N在偶数行     *
 *  （2*i），J L在奇数行（2*i-1）				   *
 ***************************************************/

void jmcc(int m, int n, int n0, float *e, float *f, float *g, float *b, float *jm, int k, int* node_PVtoPQ)
{
	FILE *fp;
	int i, j, i1, io, i0, ns;
	int pos1, pos2, pos11;
	if ((fp = fopen("output.txt", "a+")) == NULL)
		return;
	else
	{

		/* 初始化矩阵jm */
		for (i = 1; i <= n0; i++)
		{
			for (j = 1; j <= n0; j++)
			{
				jm[f2(i, j, n0)] = 0;
			}
		}

		ns = n - 1; /* 去掉一个平衡结点 */

		/* 计算式(4-89)(4-90) */
		for (i = 1; i <= ns; i++)
		{
			pos11 = f1(i);
			/* 计算式(4-90) */
			for (i1 = 1; i1 <= n; i1++)
			{
				/* pos1是式(4-90)中的j */
				pos1 = f1(i1);

				/* pos2是式(4-90)中的ij */
				pos2 = f2(i, i1, n);

				if ((i <= m) || (node_PVtoPQ[pos11] != 0)) /* i是PQ结点 */
				{
					/* 计算式(4-90)中的Jii等式右侧第一部分 */
					jm[f2(2 * i - 1, 2 * i - 1, n0)] += g[pos2] * f[pos1] + b[pos2] * e[pos1];

					/* 计算式(4-90)中的Lii等式右侧第一部分 */
					jm[f2(2 * i - 1, 2 * i, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];
				}

				/* 计算式(4-90)中的Hii等式右侧第一部分 */
				jm[f2(2 * i, 2 * i - 1, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];

				/* 计算式(4-90)中的Nii等式右侧第一部分 */
				jm[f2(2 * i, 2 * i, n0)] += -g[pos2] * f[pos1] - b[pos2] * e[pos1];
			}

			/* pos2是式(4-90)中的ii */
			pos2 = f2(i, i, n);

			/* pos1是式(4-90)中的i */
			pos1 = f1(i);

			if ((i <= m) || (node_PVtoPQ[pos1] != 0))/* i是PQ结点 */
			{
				/* 计算式(4-90)中的Jii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

				/* 计算式(4-90)中的Lii */
				jm[f2(2 * i - 1, 2 * i, n0)] += g[pos2] * e[pos1] + b[pos2] * f[pos1];
			}

			/* 计算式(4-90)中的Hii */
			jm[f2(2 * i, 2 * i - 1, n0)] += -g[pos2] * e[pos1] - b[pos2] * f[pos1];

			/* 计算式(4-90)中的Jii */
			jm[f2(2 * i, 2 * i, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

			if ((i > m) && (node_PVtoPQ[pos1] == 0))/* PV结点 */
			{
				/* 计算式(4-90)中的Rii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] = -2 * e[pos1];

				/* 计算式(4-90)中的Sii */
				jm[f2(2 * i - 1, 2 * i, n0)] = -2 * f[pos1];
			}

			/* 计算式(4-89) */
			for (j = 1; j <= ns; j++)
			{
				if (j != i)
				{
					/* pos1是式(4-89)中的i */
					pos1 = f1(i);

					/* pos2是式(4-89)中的ij */
					pos2 = f2(i, j, n);

					/* 计算式(4-89)中的Nij */
					jm[f2(2 * i, 2 * j, n0)] = b[pos2] * e[pos1] - g[pos2] * f[pos1];

					/* 计算式(4-89)中的Hij */
					jm[f2(2 * i, 2 * j - 1, n0)] = -g[pos2] * e[pos1] - b[pos2] * f[pos1];

					if ((i <= m) || (node_PVtoPQ[pos1] != 0)) /* i是PQ结点 */
					{
						/* 计算式(4-89)中的Lij (=-Hij) */
						jm[f2(2 * i - 1, 2 * j, n0)] = -jm[f2(2 * i, 2 * j - 1, n0)];

						/* 计算式(4-89)中的Jij (=Nij) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = jm[f2(2 * i, 2 * j, n0)];
					}
					else	/* i是PV结点 */
					{
						/* 计算式(4-89)中的Rij (=0) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = 0;

						/* 计算式(4-89)中的Sij (=0) */
						jm[f2(2 * i - 1, 2 * j, n0)] = 0;
					}
				}
			}
		}
		if (k != 1)
		{
			fclose(fp);
			fp = NULL;
			return;
		}

		/* 输出Jacoby矩阵 */
		fprintf(fp, "\n J                 MATRIX(C)");
		for (io = 1; io <= n0; io += 5)
		{
			i1 = (io + 4) > n0 ? n0 : (io + 4);
			fprintf(fp, "\n");
			for (j = io; j <= i1; j++)
			{
				fprintf(fp, "%10d", j);
			}
			for (i = 1; i <= n0; i++)
			{
				fprintf(fp, "\n%2d", i);
				for (j = io; j <= i1; j++)
				{
					fprintf(fp, "%12.6f", jm[f2(i, j, n0)]);
				}
			}
		}
		fprintf(fp, "\n");
		fclose(fp);
		fp = NULL;
	}
}

/**********************************************
 *     本子程序用选列主元素的高斯消元法求解组 *
 * 性方程组求各结点电压修正量,如打印参数K=1,则*
 * 输出增广矩阵变换中的上三角及电压修正量.如果*
 * 无唯一解,则给出信息,并停止程序运行.        *
 **********************************************/
void sevc(float a[], int n0, int k, int n1)
{

	FILE *fp;
	int i, j, l, n2, n3, n4, i0, io, j1, i1;
	float t0, t, c;
	if ((fp = fopen("output.txt", "a+")) == NULL)
		return;
	else
	{
		for (i = 1; i <= n0; i++)
		{
			l = i;
			for (j = i; j <= n0; j++)
			{
				if (fabs(a[f2(j, i, n1)]) > fabs(a[f2(l, i, n1)]))
				{
					l = j; /* 找到这列中的最大元 */
				}
			}
			if (l != i)
			{	/* 行交换 */
				for (j = i; j <= n1; j++)
				{
					t = a[f2(i, j, n1)];
					a[f2(i, j, n1)] = a[f2(l, j, n1)];
					a[f2(l, j, n1)] = t;
				}
			}
			if (fabs(a[f2(i, i, n1)] - 0) < 1e-10)
			{	/* 对角元近似于0, 无解 */
				printf("\nNo Solution\n");
				fclose(fp);
				fp = NULL;
				exit(-1);
			}

			t0 = a[f2(i, i, n1)];
			for (j = i; j <= n1; j++)
			{
				/* 除对角元 */
				a[f2(i, j, n1)] /= t0;
			}
			if (i == n0)
			{   /* 最后一行，不用消元 */
				continue;
			}

			/* 消元 */
			j1 = i + 1;
			for (i1 = j1; i1 <= n0; i1++)
			{
				c = a[f2(i1, i, n1)];
				for (j = i; j <= n1; j++)
				{
					a[f2(i1, j, n1)] -= a[f2(i, j, n1)] * c;
				}
			}
		}

		if (k == 1)
		{	/* 输出上三角矩阵 */
			fprintf(fp, "\nTrianglar Angmentex Matrix ");
			for (io = 1; io <= n1; io += 5)
			{
				i0 = (io + 4) > n1 ? n1 : (io + 4);
				fprintf(fp, "\n");
				fprintf(fp, "       ");
				for (i = io; i <= i0; i++)
				{
					fprintf(fp, "%12d", i);
				}
				for (i = 1; i <= n0; i++)
				{
					fprintf(fp, "\n");
					fprintf(fp, "%2d", i);
					for (j = io; j <= i0; j++)
					{
						fprintf(fp, "%15.6f", a[f2(i, j, n1)]);
					}
				}
			}
		}

		/* 回代求方程解 */
		n2 = n1 - 2;
		for (i = 1; i <= n2; i++)
		{
			n3 = n1 - i;
			for (i1 = n3; i1 <= n0; i1++)
			{
				n4 = n0 - i;
				a[f2(n4, n1, n1)] -= a[f2(i1, n1, n1)] * a[f2(n4, i1, n1)];
			}
		}

		if (k != 1)
		{
			fclose(fp);
			fp = NULL;
			return;
		}

		/* 输出电压修正值 */
		fprintf(fp, "\nVoltage correction E(i), F(i) :");
		for (io = 1; io <= n0; io += 4)
		{
			i1 = (io + 1) / 2;
			i0 = ((io + 3) / 2) > (n0 / 2) ? (n0 / 2) : ((io + 3) / 2);
			fprintf(fp, "\n");
			for (j = i1; j <= i0; j++)
			{
				fprintf(fp, "%16d%16d", j, j);
			}
			i1 = 2 * i0;
			fprintf(fp, "\n");
			for (i = io; i <= i1; i++)
			{
				fprintf(fp, "%15.6f", a[f2(i, n1, n1)]);
			}
		}
		fclose(fp);
		fp = NULL;
	}
}

/****************************************************
 *   本子程序计算线路功率,平衡节点功率,PV节点无功功 *
 * 率及线路的功率损耗并输出.如选择参数K1=1,则表示输 *
 * 入为极座标.                                      *
 ****************************************************/
#define Pi 3.1415927/180
void plsc(int n, int l, int m, float g[], float b[], float e[], float f[], \
	int e1[], int s1[], float g1[], float b1[], float c1[], float c[], \
	float co[], float p1[], float q1[], float p2[], float q2[], float p3[], \
	float q3[], float p[], float q[], float v[], float angle[], int k1, int* node_PVtoPQ)
{

	FILE *fp;
	float t1, t2, cm, x, y, z, x1, x2, y1, y2;
	int i, i1, j, m1, ns, pos1, pos2, km, st, en;
	ns = n - 1;
	if ((fp = fopen("output.txt", "w+")) == NULL)
		return;
	else
	{
		fprintf(fp, "\nTHE RESULT ARE:");
		if (k1 == 1)
		{
			for (i = 0; i < n; i++)
			{
				angle[i] *= Pi;
				e[i] = v[i] * cos(angle[i]);
				f[i] = v[i] * sin(angle[i]);
			}
		}
		t1 = 0.0; t2 = 0.0;
		for (i = 1; i <= n; i++)
		{
			pos1 = f1(i); pos2 = f2(n, i, n);
			t1 += g[pos2] * e[pos1] - b[pos2] * f[pos1];
			t2 += g[pos2] * f[pos1] + b[pos2] * e[pos1];
		}
		pos1 = f1(n);
		p[pos1] = t1 * e[pos1];
		q[pos1] = -t2 * e[pos1];
		m1 = m + 1;
		for (i1 = m1; i1 <= ns; i1++)
		{
			pos1 = f1(i1);
			if (node_PVtoPQ[pos1] == 0)//排除已转为PQ节点的PV节点
			{
				t1 = 0; t2 = 0;
				for (i = 1; i <= n; i++)
				{
					pos1 = f1(i); pos2 = f2(i1, i, n);
					t1 += g[pos2] * e[pos1] - b[pos2] * f[pos1];
					t2 += g[pos2] * f[pos1] + b[pos2] * e[pos1];
				}
				pos1 = f1(i1);
				q[pos1] = f[pos1] * t1 - e[pos1] * t2;
			}
		}
		for (i = 0; i < n; i++)
		{
			cm = co[i];
			if (cm != 0)
			{
				q[i] -= (e[i] * e[i] + f[i] * f[i])*cm;
			}
		}
		fprintf(fp, "\nBUS DATA");
		fprintf(fp, "\nBUS     VOLTAGE      ANGLE(DEGS.)      BUS P          BUS Q");
		for (i = 0; i < n; i++)
		{
			pos1 = f1(i1);
			v[i] = sqrt(e[i] * e[i] + f[i] * f[i]);
			x = e[i];
			y = f[i];
			z = y / x;
			angle[i] = atan(z);
			angle[i] /= Pi;
			fprintf(fp, "\n%3d%13.5e%15.5f%15.5e%15.5e", i + 1, v[i], angle[i], p[i], q[i]);
		}
		fprintf(fp, "\n LINE FLOW ");
		for (i = 1; i <= l; i++)
		{
			pos1 = f1(i);
			st = s1[pos1];
			en = e1[pos1];
			x1 = e[f1(st)] * e[f1(st)] + f[f1(st)] * f[f1(st)];
			x2 = e[f1(en)] * e[f1(en)] + f[f1(en)] * f[f1(en)];
			y1 = e[f1(st)] * e[f1(en)] + f[f1(st)] * f[f1(en)];
			y2 = f[f1(st)] * e[f1(en)] - e[f1(st)] * f[f1(en)];
			p1[pos1] = (x1 - y1)*g1[pos1] - y2 * b1[pos1];
			q1[pos1] = -x1 * (c1[pos1] + b1[pos1]) + y1 * b1[pos1] - y2 * g1[pos1];
			p2[pos1] = (x2 - y1)*g1[pos1] + y2 * b1[pos1];
			q2[pos1] = -x2 * (c1[pos1] + b1[pos1]) + y1 * b1[pos1] + y2 * g1[pos1];
			for (j = 1; j <= n; j++)
			{
				cm = c[f2(j, i, l)];
				if (cm != 0.0)
				{
					km = 1;
					if (en == j)
					{
						km = 2;
					}
					if (km == 1)
					{
						q1[pos1] -= (e[f1(j)] * e[f1(j)] + f[f1(j)] * f[f1(j)])*cm;
					}
					else
					{
						q2[pos1] -= (e[f1(j)] * e[f1(j)] + f[f1(j)] * f[f1(j)])*cm;
					}
				}
			}
			p3[pos1] = p1[pos1] + p2[pos1];
			q3[pos1] = q1[pos1] + q2[pos1];
			fprintf(fp, "\n%2d%8d%11d%13.6e%13.6e%13.6e%13.6e%17d%11d%13.6e%13.6e", \
				i, s1[pos1], e1[pos1], p1[pos1], q1[pos1], p3[pos1], q3[pos1], \
				e1[pos1], s1[pos1], p2[pos1], q2[pos1]);

		}
		fclose(fp);
		fp = NULL;
	}
}
