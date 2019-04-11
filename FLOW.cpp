/**************************FLOW.C*********************************/

/*******************************************************************
 *     �����ṩ���ǵ���ϵͳ����������ⷨ������ӳ���,���õķ����� *
 *  Newton_Raphson��.                                              *
 *     ���������õı���˵������:                                   *
 *       N:����ڵ�����.         M:�����PQ�ڵ���.                 *
 *       L:�����֧·����.       N0:�ſɱȾ��������.              *
 *       N1:N0+1                 K:��ӡ����.K=1,���ӡ;����,����ӡ.*
 *       K1:�ӳ���PLSC���ж������ѹ����ʽ.K1=1,��Ϊ��������ʽ.����*
 *          Ϊֱ��������ʽ.                                        *
 *       D:�й����޹������������ֵ������ֵ��.                   *
 *       G(I,J):Ybus�ĵ絼Ԫ��(ʵ��).                              *
 *       B(I,J):Ybus�ĵ���Ԫ��(�鲿).                              *
 *       G1(I) :��I֧·�Ĵ����絼.      B1(I):��I֧·�Ĵ�������.   *
 *       C1(I) :��I֧·��pie�ͶԳƽӵص���.                        *
 *       C(I,J):��I�ڵ�J֧·���Գƽӵص���.                        *
 *       CO(I) :��I�ڵ�Ľӵص���.                                 *
 *       S1(I) :��I֧·����ʼ�ڵ��.    E1(I):��I֧·����ֹ�ڵ��. *
 *       P(I)  :��I�ڵ��ע���й�����.  Q(I):��I�ڵ��ע���޹�����.*
 *       P0(I) :��I�ڵ��й��������.    Q0(I):��I�ڵ��޹��������. *
 *       V0(I) :��I�ڵ�(PV�ڵ�)�ĵ�ѹ���(ƽ�����).               *
 *       V(I)  :��I�ڵ�ĵ�ѹ��ֵ.                                 *
 *       E(I)  :��I�ڵ�ĵ�ѹ��ʵ��.    F(I):��I�ڵ�ĵ�ѹ���鲿.  *
 *      JM(I,J):Jacoby����ĵ�I��J��Ԫ��.                          *
 *       A(I,J):�������̵��������,���ǻ�����ĵ�I��J��Ԫ��,����� *
 *              ����A��������һ�д�������Ľ�.                   *
 *       P1(I) :��I֧·��S1(I)�ڵ�ע����й�����.                  *
 *       Q1(I) :��I֧·��S1(I)�ڵ�ע����޹�����.                  *
 *       P2(I) :��I֧·��E1(I)�ڵ�ע����й�����.                  *
 *       Q2(I) :��I֧·��E1(I)�ڵ�ע����޹�����.                  *
 *       P3(I) :��I֧·���й��������.                             *
 *       Q3(I) :��I֧·���޹��������.                             *
 *     ANGLE(I):��I�ڵ��ѹ�ĽǶ�.                                 *
 *     �ڵ���˳��N���ڵ��У�ǰM��ΪPQ�ڵ㣬M+1��N-1���ڵ�Ϊ    *
 *                   PV�ڵ㣬��N��Ϊƽ��ڵ�                       *
*******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#pragma warning( disable : 4996)
#define f1(i) (i-1)  
/* ��ϰ�ߵ�һ�׾�����±�ת��ΪC���������±�*/

#define f2(i,j,n) ((i-1)*(n)+j-1)
/* ��ϰ�ߵĶ��׾�����±�ת��ΪC���������±�*/

/****************************************************
 *    ���ӳ������������֧·���ɼ��й���Ϣ,�γɽ�� *
 * ���ɾ���,���ӡ����K=1,������絼����G�͵��ɾ�B  *
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

		/* ��ʼ������G,B */
		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= n; j++)
			{
				pos2 = f2(i, j, n);
				g[pos2] = 0; b[pos2] = 0;
			}
		}

		/* ����֧·���� */
		for (i = 1; i <= l; i++)
		{
			/* ����Խ�Ԫ */
			pos1 = f1(i);
			st = s1[pos1]; en = e1[pos1];
			pos2 = f2(st, st, n);
			g[pos2] += g1[pos1];
			b[pos2] += b1[pos1] + c1[pos1];
			pos2 = f2(en, en, n);
			g[pos2] += g1[pos1];
			b[pos2] += b1[pos1] + c1[pos1];

			/* ����ǶԽ�Ԫ */
			pos2 = f2(st, en, n);
			g[pos2] -= g1[pos1];
			b[pos2] -= b1[pos1];
			g[f2(en, st, n)] = g[f2(st, en, n)];
			b[f2(en, st, n)] = b[f2(st, en, n)];
		}

		/* ����ӵ�֧·���� */
		for (i = 1; i <= n; i++)
		{
			/*���ԳƲ��֡�*/
			b[f2(i, i, n)] += co[f1(i)];

			/* �ǶԳƲ��֡�*/
			for (j = 1; j <= l; j++)
			{
				b[f2(i, i, n)] += c[f2(i, j, l)];
			}
		}
		if (k != 1)
		{
			fclose(fp);
			fp = NULL;
			return; /* ���K��Ϊ 1,�򷵻�;����,��ӡ���ɾ��� */
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
 *     ���ӳ�����������Ĺ��ʼ���ѹ������  *
 * ������ʼ���ѹ�����,����������й����� *
 * ��������������Ƚ�.���ӡ����K=1,���� *
 * ��P0,Q0(��PQ���),V0(��PV���).         *
 * ��Ӧ��¼һP177ʽ(4-86)(4-87)			   *
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
				/* a1,a2��Ӧ��¼һP177ʽ(4-86)�������ڵ�ʽ�� */
				pos2 = f2(i, j, n);
				a1 += g[pos2] * e[f1(j)] - b[pos2] * f[f1(j)];
				a2 += g[pos2] * f[f1(j)] + b[pos2] * e[f1(j)];
			}
			/* ����ʽ(4-86)(4-87)�е�deltaPi��*/
			p0[pos1] = p[pos1] - e[pos1] * a1 - f[pos1] * a2;
			if (i <= m)
			{	/* ����PQ����е�deltaQi��*/
				q0[pos1] = q[pos1] - f[pos1] * a1 + e[pos1] * a2;
			}
			else if ((node_PVtoPQ[pos1] == 1) || (node_PVtoPQ[pos1] == -1))
			{	/*����PVtoPQ�ڵ��deltaQi*/
				q0[pos1] = q[pos1] - f[pos1] * a1 + e[pos1] * a2;
			}
			else
			{	/* ����PV����е�deltaViƽ����*/
				v0[pos1] = v[pos1] * v[pos1] - e[pos1] * e[pos1] - f[pos1] * f[pos1];
			}

			/* ������ */
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

		/* �ҵ�deltaP��deltaQ�е�����ߣ���Ϊ����ָ��, ����dd�� */
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
 *    ���ӳ�����ݽڵ㵼�ɼ���ѹ��Jacoby����,������*
 *  ��ѹ������,���ӡ����K=1,�����Jacoby����.     *
 *  ��Ӧ�ڸ�¼һP178ʽ(4-89)(4-90)				   *
 *    ֵ��ע����ǣ�������Jacobi����H N J L������˳*
 *  ����ʽ��4��88�����в�ͬ��������H N��ż����     *
 *  ��2*i����J L�������У�2*i-1��				   *
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

		/* ��ʼ������jm */
		for (i = 1; i <= n0; i++)
		{
			for (j = 1; j <= n0; j++)
			{
				jm[f2(i, j, n0)] = 0;
			}
		}

		ns = n - 1; /* ȥ��һ��ƽ���� */

		/* ����ʽ(4-89)(4-90) */
		for (i = 1; i <= ns; i++)
		{
			pos11 = f1(i);
			/* ����ʽ(4-90) */
			for (i1 = 1; i1 <= n; i1++)
			{
				/* pos1��ʽ(4-90)�е�j */
				pos1 = f1(i1);

				/* pos2��ʽ(4-90)�е�ij */
				pos2 = f2(i, i1, n);

				if ((i <= m) || (node_PVtoPQ[pos11] != 0)) /* i��PQ��� */
				{
					/* ����ʽ(4-90)�е�Jii��ʽ�Ҳ��һ���� */
					jm[f2(2 * i - 1, 2 * i - 1, n0)] += g[pos2] * f[pos1] + b[pos2] * e[pos1];

					/* ����ʽ(4-90)�е�Lii��ʽ�Ҳ��һ���� */
					jm[f2(2 * i - 1, 2 * i, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];
				}

				/* ����ʽ(4-90)�е�Hii��ʽ�Ҳ��һ���� */
				jm[f2(2 * i, 2 * i - 1, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];

				/* ����ʽ(4-90)�е�Nii��ʽ�Ҳ��һ���� */
				jm[f2(2 * i, 2 * i, n0)] += -g[pos2] * f[pos1] - b[pos2] * e[pos1];
			}

			/* pos2��ʽ(4-90)�е�ii */
			pos2 = f2(i, i, n);

			/* pos1��ʽ(4-90)�е�i */
			pos1 = f1(i);

			if ((i <= m) || (node_PVtoPQ[pos1] != 0))/* i��PQ��� */
			{
				/* ����ʽ(4-90)�е�Jii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

				/* ����ʽ(4-90)�е�Lii */
				jm[f2(2 * i - 1, 2 * i, n0)] += g[pos2] * e[pos1] + b[pos2] * f[pos1];
			}

			/* ����ʽ(4-90)�е�Hii */
			jm[f2(2 * i, 2 * i - 1, n0)] += -g[pos2] * e[pos1] - b[pos2] * f[pos1];

			/* ����ʽ(4-90)�е�Jii */
			jm[f2(2 * i, 2 * i, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

			if ((i > m) && (node_PVtoPQ[pos1] == 0))/* PV��� */
			{
				/* ����ʽ(4-90)�е�Rii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] = -2 * e[pos1];

				/* ����ʽ(4-90)�е�Sii */
				jm[f2(2 * i - 1, 2 * i, n0)] = -2 * f[pos1];
			}

			/* ����ʽ(4-89) */
			for (j = 1; j <= ns; j++)
			{
				if (j != i)
				{
					/* pos1��ʽ(4-89)�е�i */
					pos1 = f1(i);

					/* pos2��ʽ(4-89)�е�ij */
					pos2 = f2(i, j, n);

					/* ����ʽ(4-89)�е�Nij */
					jm[f2(2 * i, 2 * j, n0)] = b[pos2] * e[pos1] - g[pos2] * f[pos1];

					/* ����ʽ(4-89)�е�Hij */
					jm[f2(2 * i, 2 * j - 1, n0)] = -g[pos2] * e[pos1] - b[pos2] * f[pos1];

					if ((i <= m) || (node_PVtoPQ[pos1] != 0)) /* i��PQ��� */
					{
						/* ����ʽ(4-89)�е�Lij (=-Hij) */
						jm[f2(2 * i - 1, 2 * j, n0)] = -jm[f2(2 * i, 2 * j - 1, n0)];

						/* ����ʽ(4-89)�е�Jij (=Nij) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = jm[f2(2 * i, 2 * j, n0)];
					}
					else	/* i��PV��� */
					{
						/* ����ʽ(4-89)�е�Rij (=0) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = 0;

						/* ����ʽ(4-89)�е�Sij (=0) */
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

		/* ���Jacoby���� */
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
 *     ���ӳ�����ѡ����Ԫ�صĸ�˹��Ԫ������� *
 * �Է������������ѹ������,���ӡ����K=1,��*
 * ����������任�е������Ǽ���ѹ������.���*
 * ��Ψһ��,�������Ϣ,��ֹͣ��������.        *
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
					l = j; /* �ҵ������е����Ԫ */
				}
			}
			if (l != i)
			{	/* �н��� */
				for (j = i; j <= n1; j++)
				{
					t = a[f2(i, j, n1)];
					a[f2(i, j, n1)] = a[f2(l, j, n1)];
					a[f2(l, j, n1)] = t;
				}
			}
			if (fabs(a[f2(i, i, n1)] - 0) < 1e-10)
			{	/* �Խ�Ԫ������0, �޽� */
				printf("\nNo Solution\n");
				fclose(fp);
				fp = NULL;
				exit(-1);
			}

			t0 = a[f2(i, i, n1)];
			for (j = i; j <= n1; j++)
			{
				/* ���Խ�Ԫ */
				a[f2(i, j, n1)] /= t0;
			}
			if (i == n0)
			{   /* ���һ�У�������Ԫ */
				continue;
			}

			/* ��Ԫ */
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
		{	/* ��������Ǿ��� */
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

		/* �ش��󷽳̽� */
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

		/* �����ѹ����ֵ */
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
 *   ���ӳ��������·����,ƽ��ڵ㹦��,PV�ڵ��޹��� *
 * �ʼ���·�Ĺ�����Ĳ����.��ѡ�����K1=1,���ʾ�� *
 * ��Ϊ������.                                      *
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
			if (node_PVtoPQ[pos1] == 0)//�ų���תΪPQ�ڵ��PV�ڵ�
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
