#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <process.h>
#include <time.h>
#include <math.h>
#define pi  3.1415926535897932384626433832795


unsigned long int nn1 = 1;
unsigned long int nn2 = 1;
unsigned long int nn3 = 1;
unsigned long int nn4 = 1;
srand1(unsigned long int s) { nn1 = s; }
srand2(unsigned long int s) { nn2 = s; }
srand3(unsigned long int s) { nn3 = s; }
srand4(unsigned long int s) { nn4 = s; }

unsigned long int *tm1;
unsigned long int tm0, tm2, tm2_old;

int* _cal_1dint(int *a, int N)
{
	if ((a = (int*)calloc(N, sizeof(int))) == NULL)
	{
		printf("Memory is not sufficient");
		exit(-1);
	}
	return(a);
}

double rand1()
{
	nn1 = (1664525 * nn1 + 1013904223);
	return nn1 / 4294967296.0;
}

double rand2()
{
	nn2 = (171 * nn2) % 30269;
	nn3 = (172 * nn3) % 30307;
	nn4 = (170 * nn4) % 30323;
	return(fmod(nn2 / 30269.0 + nn3 / 30307.0 + nn4 / 30323.0, 1.0));
}



double* _cal_1d(double *a, int N)
{
	if ((a = (double*)calloc(N, sizeof(double))) == NULL)
	{
		printf("Memory is not sufficient");
		exit(-1);
	}
	return(a);
}

double** _cal_2d(double **a, int N1, int N2)
{
	int k;
	if ((a = (double**)calloc(N1, sizeof(double*))) == NULL)
	{
		printf("Memory is not sufficient");
		exit(-1);
	}

	for (k = 0; k < N1; k++)
	{
		a[k] = _cal_1d(a[k], N2);
	}
	return(a);
}


double cbrt1(double a)
{
	double r, x;
	x = 1.0; r = 10.0;
	while (fabs(x - r) > 1.0e-12) { x = r; r = (a / (x*x) + 2.0*x) / 3.0; }
	return(r);
}


double min3d(double a, double b, double c, double d)
{
	double m;
	if (fabs(a - d) <= fabs(b - d) && fabs(a - d) <= fabs(c - d)) { m = a; }
	if (fabs(c - d) <= fabs(b - d) && fabs(c - d) <= fabs(a - d)) { m = c; }
	if (fabs(b - d) <= fabs(a - d) && fabs(b - d) <= fabs(c - d)) { m = b; }
	return(m);
}

double root_cube(double B, double C, double D, double CC, double CCC, double CCCC)
{
	double p, q, fi, delta, tmp;
	double root_cardano, root_trig, root_trig1, root_trig2, root_trig3;

	p = C - B * B / 3; q = 2.0*B*B*B / 27.0 - B * C / 3.0 + D;
	delta = 0.25*q*q + p * p*p / 27.0;
	printf("\n delta=%le", delta);
	if (delta > 0)
	{
		tmp = sqrt(delta);
		root_cardano = cbrt1(-0.5*q + tmp) + cbrt1(-0.5*q - tmp) - B / 3.0;
		printf("\n root_cardano=%le", root_cardano);
		return(root_cardano);
	}

	if (delta <= 0)
	{
		fi = asin(9.0*q*sqrt(-4.0*p / 3.0) / (4.0*p*p));
		tmp = sqrt(-4.0*p / 3.0);
		root_trig1 = tmp * sin(fi / 3.0) - B / 3.0;
		root_trig2 = tmp * sin(fi / 3.0 + 2.0*pi / 3.0) - B / 3.0;
		root_trig3 = tmp * sin(fi / 3.0 + 4.0*pi / 3.0) - B / 3.0;
		if (2.0*CC + 6.0*CCC*root_trig1 + 12.0*CCCC*root_trig1*root_trig1 > 0) { root_trig = root_trig1; }
		if (2.0*CC + 6.0*CCC*root_trig2 + 12.0*CCCC*root_trig2*root_trig2 > 0) { root_trig = root_trig2; }
		if (2.0*CC + 6.0*CCC*root_trig3 + 12.0*CCCC*root_trig3*root_trig3 > 0) { root_trig = root_trig3; }
		printf("\n root_trig=%le", root_trig);
		return(root_trig);
	}
}

double root_cube1(double B, double C, double D, double K)
{
	double p, q, fi, delta, tmp;
	double root_cardano, root_trig, root_trig1, root_trig2, root_trig3;

	p = C - B * B / 3; q = 2.0*B*B*B / 27.0 - B * C / 3.0 + D;
	delta = 0.25*q*q + p * p*p / 27.0;
	printf("\n delta=%le", delta);
	if (delta > 0)
	{
		tmp = sqrt(delta);
		root_cardano = cbrt1(-0.5*q + tmp) + cbrt1(-0.5*q - tmp) - B / 3.0;
		printf("\n root_cardano=%le", root_cardano);
		return(root_cardano);
	}

	if (delta <= 0)
	{
		fi = asin(9.0*q*sqrt(-4.0*p / 3.0) / (4.0*p*p));
		tmp = sqrt(-4.0*p / 3.0);
		root_trig1 = tmp * sin(fi / 3.0) - B / 3.0;
		root_trig2 = tmp * sin(fi / 3.0 + 2.0*pi / 3.0) - B / 3.0;
		root_trig3 = tmp * sin(fi / 3.0 + 4.0*pi / 3.0) - B / 3.0;
		root_trig = min3d(root_trig1, root_trig2, root_trig3, K);
		printf("\n root_trig=%le", root_trig);
		return(root_trig);
	}
}



double determinant55(double **a)
{
	return(
		(a[0][4] * a[1][3] * a[2][2] * a[3][1] - a[0][3] * a[1][4] * a[2][2] * a[3][1] -
			a[0][4] * a[1][2] * a[2][3] * a[3][1] +
			a[0][2] * a[1][4] * a[2][3] * a[3][1] +
			a[0][3] * a[1][2] * a[2][4] * a[3][1] -
			a[0][2] * a[1][3] * a[2][4] * a[3][1] -
			a[0][4] * a[1][3] * a[2][1] * a[3][2] +
			a[0][3] * a[1][4] * a[2][1] * a[3][2] +
			a[0][4] * a[1][1] * a[2][3] * a[3][2] -
			a[0][1] * a[1][4] * a[2][3] * a[3][2] -
			a[0][3] * a[1][1] * a[2][4] * a[3][2] +
			a[0][1] * a[1][3] * a[2][4] * a[3][2] +
			a[0][4] * a[1][2] * a[2][1] * a[3][3] -
			a[0][2] * a[1][4] * a[2][1] * a[3][3] -
			a[0][4] * a[1][1] * a[2][2] * a[3][3] +
			a[0][1] * a[1][4] * a[2][2] * a[3][3] +
			a[0][2] * a[1][1] * a[2][4] * a[3][3] -
			a[0][1] * a[1][2] * a[2][4] * a[3][3] -
			a[0][3] * a[1][2] * a[2][1] * a[3][4] +
			a[0][2] * a[1][3] * a[2][1] * a[3][4] +
			a[0][3] * a[1][1] * a[2][2] * a[3][4] -
			a[0][1] * a[1][3] * a[2][2] * a[3][4] -
			a[0][2] * a[1][1] * a[2][3] * a[3][4] +
			a[0][1] * a[1][2] * a[2][3] * a[3][4])* a[4][0] - (a[0][4] * a[1][3] * a[2][2] * a[3][0] -
				a[0][3] * a[1][4] * a[2][2] * a[3][0] -
				a[0][4] * a[1][2] * a[2][3] * a[3][0] +
				a[0][2] * a[1][4] * a[2][3] * a[3][0] +
				a[0][3] * a[1][2] * a[2][4] * a[3][0] -
				a[0][2] * a[1][3] * a[2][4] * a[3][0] -
				a[0][4] * a[1][3] * a[2][0] * a[3][2] +
				a[0][3] * a[1][4] * a[2][0] * a[3][2] +
				a[0][4] * a[1][0] * a[2][3] * a[3][2] -
				a[0][0] * a[1][4] * a[2][3] * a[3][2] -
				a[0][3] * a[1][0] * a[2][4] * a[3][2] +
				a[0][0] * a[1][3] * a[2][4] * a[3][2] +
				a[0][4] * a[1][2] * a[2][0] * a[3][3] -
				a[0][2] * a[1][4] * a[2][0] * a[3][3] -
				a[0][4] * a[1][0] * a[2][2] * a[3][3] +
				a[0][0] * a[1][4] * a[2][2] * a[3][3] +
				a[0][2] * a[1][0] * a[2][4] * a[3][3] -
				a[0][0] * a[1][2] * a[2][4] * a[3][3] -
				a[0][3] * a[1][2] * a[2][0] * a[3][4] +
				a[0][2] * a[1][3] * a[2][0] * a[3][4] +
				a[0][3] * a[1][0] * a[2][2] * a[3][4] -
				a[0][0] * a[1][3] * a[2][2] * a[3][4] -
				a[0][2] * a[1][0] * a[2][3] * a[3][4] +
				a[0][0] * a[1][2] * a[2][3] * a[3][4])* a[4][1] + (a[0][4] * a[1][3] * a[2][1] * a[3][0] -
					a[0][3] * a[1][4] * a[2][1] * a[3][0] -
					a[0][4] * a[1][1] * a[2][3] * a[3][0] +
					a[0][1] * a[1][4] * a[2][3] * a[3][0] +
					a[0][3] * a[1][1] * a[2][4] * a[3][0] -
					a[0][1] * a[1][3] * a[2][4] * a[3][0] -
					a[0][4] * a[1][3] * a[2][0] * a[3][1] +
					a[0][3] * a[1][4] * a[2][0] * a[3][1] +
					a[0][4] * a[1][0] * a[2][3] * a[3][1] -
					a[0][0] * a[1][4] * a[2][3] * a[3][1] -
					a[0][3] * a[1][0] * a[2][4] * a[3][1] +
					a[0][0] * a[1][3] * a[2][4] * a[3][1] +
					a[0][4] * a[1][1] * a[2][0] * a[3][3] -
					a[0][1] * a[1][4] * a[2][0] * a[3][3] -
					a[0][4] * a[1][0] * a[2][1] * a[3][3] +
					a[0][0] * a[1][4] * a[2][1] * a[3][3] +
					a[0][1] * a[1][0] * a[2][4] * a[3][3] -
					a[0][0] * a[1][1] * a[2][4] * a[3][3] -
					a[0][3] * a[1][1] * a[2][0] * a[3][4] +
					a[0][1] * a[1][3] * a[2][0] * a[3][4] +
					a[0][3] * a[1][0] * a[2][1] * a[3][4] -
					a[0][0] * a[1][3] * a[2][1] * a[3][4] -
					a[0][1] * a[1][0] * a[2][3] * a[3][4] +
					a[0][0] * a[1][1] * a[2][3] * a[3][4])* a[4][2] - (a[0][4] * a[1][2] * a[2][1] * a[3][0] -
						a[0][2] * a[1][4] * a[2][1] * a[3][0] -
						a[0][4] * a[1][1] * a[2][2] * a[3][0] +
						a[0][1] * a[1][4] * a[2][2] * a[3][0] +
						a[0][2] * a[1][1] * a[2][4] * a[3][0] -
						a[0][1] * a[1][2] * a[2][4] * a[3][0] -
						a[0][4] * a[1][2] * a[2][0] * a[3][1] +
						a[0][2] * a[1][4] * a[2][0] * a[3][1] +
						a[0][4] * a[1][0] * a[2][2] * a[3][1] -
						a[0][0] * a[1][4] * a[2][2] * a[3][1] -
						a[0][2] * a[1][0] * a[2][4] * a[3][1] +
						a[0][0] * a[1][2] * a[2][4] * a[3][1] +
						a[0][4] * a[1][1] * a[2][0] * a[3][2] -
						a[0][1] * a[1][4] * a[2][0] * a[3][2] -
						a[0][4] * a[1][0] * a[2][1] * a[3][2] +
						a[0][0] * a[1][4] * a[2][1] * a[3][2] +
						a[0][1] * a[1][0] * a[2][4] * a[3][2] -
						a[0][0] * a[1][1] * a[2][4] * a[3][2] -
						a[0][2] * a[1][1] * a[2][0] * a[3][4] +
						a[0][1] * a[1][2] * a[2][0] * a[3][4] +
						a[0][2] * a[1][0] * a[2][1] * a[3][4] -
						a[0][0] * a[1][2] * a[2][1] * a[3][4] -
						a[0][1] * a[1][0] * a[2][2] * a[3][4] +
						a[0][0] * a[1][1] * a[2][2] * a[3][4])* a[4][3] + (a[0][3] * a[1][2] * a[2][1] * a[3][0] -
							a[0][2] * a[1][3] * a[2][1] * a[3][0] -
							a[0][3] * a[1][1] * a[2][2] * a[3][0] +
							a[0][1] * a[1][3] * a[2][2] * a[3][0] +
							a[0][2] * a[1][1] * a[2][3] * a[3][0] -
							a[0][1] * a[1][2] * a[2][3] * a[3][0] -
							a[0][3] * a[1][2] * a[2][0] * a[3][1] +
							a[0][2] * a[1][3] * a[2][0] * a[3][1] +
							a[0][3] * a[1][0] * a[2][2] * a[3][1] -
							a[0][0] * a[1][3] * a[2][2] * a[3][1] -
							a[0][2] * a[1][0] * a[2][3] * a[3][1] +
							a[0][0] * a[1][2] * a[2][3] * a[3][1] +
							a[0][3] * a[1][1] * a[2][0] * a[3][2] -
							a[0][1] * a[1][3] * a[2][0] * a[3][2] -
							a[0][3] * a[1][0] * a[2][1] * a[3][2] +
							a[0][0] * a[1][3] * a[2][1] * a[3][2] +
							a[0][1] * a[1][0] * a[2][3] * a[3][2] -
							a[0][0] * a[1][1] * a[2][3] * a[3][2] -
							a[0][2] * a[1][1] * a[2][0] * a[3][3] +
							a[0][1] * a[1][2] * a[2][0] * a[3][3] +
							a[0][2] * a[1][0] * a[2][1] * a[3][3] -
							a[0][0] * a[1][2] * a[2][1] * a[3][3] -
							a[0][1] * a[1][0] * a[2][2] * a[3][3] +
							a[0][0] * a[1][1] * a[2][2] * a[3][3])*a[4][4]);
}
double **chtriha(double **a, int k)
{
	double  **b = NULL;
	b = _cal_2d(b, 5, 5);
	int i, j;
	for (i = k + 1; i <= 4; i++)
	{
		for (j = k + 1; j <= 4; j++)
		{
			b[i][j] = a[i][j] - a[i][k] * a[k][j] / a[k][k];
		}

	}
	return (b);
}
double *chtrihb(double **a, double *b, int k)
{
	double  *c = NULL;
	c = _cal_1d(c, 5);
	int i, j;

	for (j = k + 1; j <= 4; j++)
	{
		c[j] = b[j] - a[j][k] * b[k] / a[k][k];
	}
	return (c);
}
double *solve(double d0, double d1, double d2, double d3, double d4, double **matr)
{
	double **matr1 = NULL, **matr2 = NULL, **matr3 = NULL, **matr4 = NULL;
	matr1 = _cal_2d(matr1, 5, 5);
	matr2 = _cal_2d(matr2, 5, 5);
	matr3 = _cal_2d(matr3, 5, 5);
	matr4 = _cal_2d(matr4, 5, 5);

	double *d00 = NULL, *d11 = NULL, *d22 = NULL, *d33 = NULL, *d44 = NULL, *Csolve = NULL;
	d00 = _cal_1d(d00, 5);
	d11 = _cal_1d(d11, 5);
	d22 = _cal_1d(d22, 5);
	d33 = _cal_1d(d33, 5);
	d44 = _cal_1d(d44, 5);
	Csolve = _cal_1d(Csolve, 5);

	double c44, c33, c22, c11, c00;
	d00[0] = d0; d00[1] = d1; d00[2] = d2; d00[3] = d3; d00[4] = d4;
	d11 = chtrihb(matr, d00, 0);
	matr1 = chtriha(matr, 0);
	d22 = chtrihb(matr1, d11, 1);
	matr2 = chtriha(matr1, 1);
	d33 = chtrihb(matr2, d22, 2);
	matr3 = chtriha(matr2, 2);
	d44 = chtrihb(matr3, d33, 3);
	matr4 = chtriha(matr3, 3);
	c44 = d44[4] / matr4[4][4];
	c33 = (1.0 / matr3[3][3])*(d33[3] - matr3[3][4] * c44);
	c22 = (1.0 / matr2[2][2])*(d22[2] - matr2[2][3] * c33 - matr2[2][4] * c44);
	c11 = (1.0 / matr1[1][1])*(d11[1] - matr1[1][2] * c22 - matr1[1][3] * c33 - matr1[1][4] * c44);
	c00 = (1.0 / matr[0][0])*(d00[0] - matr[0][1] * c11 - matr[0][2] * c22 - matr[0][3] * c33 - matr[0][4] * c44);
	Csolve[0] = c00; Csolve[1] = c11; Csolve[2] = c22; Csolve[3] = c33; Csolve[4] = c44;
	return (Csolve);

}
#define CN 8
static double cof[CN] = {
	   2.5066282746310005,
	   1.0000000000190015,
	   76.18009172947146,
	  -86.50532032941677,
	   24.01409824083091,
	  -1.231739572450155,
	   0.1208650973866179e-2,
	  -0.5395239384953e-5,
};

/* logarithm of gamma-function */
double GammLn(double x)
{
	double y, ser, *co;
	int j;
	/* calculate the series */
	ser = cof[1]; y = x; co = cof + 2;
	for (j = 2; j < CN; j++) {
		y += 1.; ser += (*co) / y; co++;
	}
	/* and the other parts of the function */
	y = x + 5.5;
	y -= (x + 0.5)*log(y);
	return(-y + log(cof[0] * ser / x));
}

/* the gamma-function itself */
double Gamma(double x)
{
	return(exp(GammLn(x)));
}

/* the beta function */
double Beta(double x, double y)
{
	return(exp(GammLn(x) + GammLn(y) - GammLn(x + y)));
}
double Ymaxx(double q, double B, double G)
{
	double y;
	y = pow((pow(2.0, (q - 1)) - 1.0), (0.5)) / (G*B);
	return(y);
}
void main()
{
	int i, j, NN, check, smoothing, noise_type, noise_var, i1_max, i1_min, Nsim, is_min, is_max;
	double tmp, tmp1, tmp2, tmp3, tmp4, gamma;
	double  q, G, G2, a, b, x, x0, x1, B, dB, B0;
	double *Y = NULL, *Y1 = NULL, Y2;
	double *Y1c = NULL, *Y2c = NULL, *XXs = NULL;
	double qopt = 0, G_opt;
	double errorlast = 0.0, Ymax_opt=0.0, B0_opt=0.0;
	int Nsimopt=0;
	double Y_old, Y1_old, Ymax, Y2_old;
	double Y_max, Y1_max, Y1_min, x1_min, x1_max;
	double I1, I10, App, delta_Bpp, Y_x1, error_App, error_q, error_q_min, error_q_old, q_error_min, dq,dG;
	double koef_Bscan, koef_Bscan_min, koef_Bscan_max, d_koef_Bscan;
	double *Y1s = NULL, *sk = NULL, *XX = NULL, *Xr = NULL;
	double Ymax_min, x0_min, G_min;
	double c0, c1, c2, c3, c4, f1, f2, f3, f4, f5, f6, f7, f8, d0, d1, d2, d3, d4, Y1_max_s, Xr_max, x1_max_s, Y1_min_s, Xr_min, x1_min_s;
	double noise_ampl;
	int i1_scan, j1_min;
	double err_parab_min, err_parab;
	double error_Xr, Xr_tmp, deter, deter0, deter1, deter2, deter3, deter4, **matr = NULL;
	double *errorQ = NULL, *QN = NULL,*XXish=NULL,*Y1ish=NULL, *GN = NULL;
	errorQ = _cal_1d(errorQ, 400);
	QN = _cal_1d(QN, 400);
	GN = _cal_1d(QN, 400);

	FILE *fp = NULL, *fp1 = NULL ,*fp3 = NULL, *fp2 = NULL, *fp4 = NULL;
	fp1 = fopen("Polinom.dat", "w");
	fp3 = fopen("resultQEXP.dat", "w");
	fp2 = fopen("resultspectr.dat", "w");
	fp4 = fopen("resultQoptEXP.dat", "w");
	double q0,dBish;
	int k = 0;
	//NN = 3999;
	NN = 8000;
	smoothing = 0;
	int flag = 0;
	Nsim = 60; i1_scan = 2;


	sk = _cal_1d(sk, 5);

	if (smoothing == 3) { sk[0] = 1.0 / 3.0; sk[1] = 1.0 / 3.0; sk[2] = 1.0 / 3.0; }

	if (smoothing == 5) { sk[0] = 0.2; sk[1] = 0.2; sk[2] = 0.2; sk[3] = 0.2; sk[4] = 0.2; }

	if (smoothing == 6) { sk[0] = -0.086; sk[1] = 0.343; sk[2] = 0.486; sk[3] = 0.343; sk[4] = -0.086; }

	Y = _cal_1d(Y, NN);
	Y1 = _cal_1d(Y1, NN);
	Y1c = _cal_1d(Y1c, NN);
	Y2c = _cal_1d(Y2c, NN);
	Y1s = _cal_1d(Y1s, NN);
	XXs = _cal_1d(XXs, NN);
	XX = _cal_1d(XX, NN);
	Xr = _cal_1d(Xr, NN);
	matr = _cal_2d(matr, 5, 5);

	XXish = _cal_1d(XXish, NN);
	Y1ish = _cal_1d(Y1ish, NN);

	double *Csolve = NULL;
	Csolve = _cal_1d(Csolve, 5);
	double aa, bb;
	fp = fopen("RT-fc=3300-ma=0,5-rg=500-pw=1.dat", "r");
	//fp = fopen("BDPA-3253-20-pw=1-rg=50-ma=0,08.dat", "r");
	//fp = fopen("MyNoiseTestSample.dat", "r");
	for (i = 0; i < NN; i++)
	{
		fscanf(fp, "%le%le", &aa, &bb);
		XX[i] = aa; Y1[i] = bb;
		XXish[i] = aa; Y1ish[i] = bb;
	}
	fclose(fp);




	if (smoothing == 5 || smoothing == 6) {
		Y1s[0] = Y1[0]; Y1s[1] = Y1[1];
		for (i = 2; i < NN - 2; i++)
		{
			Y1s[i] = sk[0] * Y1[i - 2] + sk[1] * Y1[i - 1] + sk[2] * Y1[i] + sk[3] * Y1[i + 1] + sk[4] * Y1[i + 2];
		}
		Y1s[NN - 2] = Y1[NN - 2]; Y1s[NN - 1] = Y1[NN - 1];

		for (i = 0; i < NN; i++)
		{
			Y1[i] = Y1s[i];
		}
	}

	if ((smoothing == 3))
	{
		Y1s[0] = Y1[0];
		for (i = 1; i < NN - 1; i++)
		{
			Y1s[i] = sk[0] * Y1[i - 1] + sk[1] * Y1[i] + sk[2] * Y1[i + 1];
		}
		Y1s[NN - 1] = Y1[NN - 1];

		for (i = 0; i < NN; i++)
		{
			Y1[i] = Y1s[i];
		}

	}
	if (flag == 1)
	{
		Y1s[0] = Y1[0];
		for (i = 1; i < NN - 1; i = i + 3)
		{
			Y1s[(i + 2) / 3] = Y1[i];
			XXs[(i + 2) / 3] = XX[i];
		}
		Y1s[NN - 1] = Y1[NN - 1];
		NN = NN / 3;
		for (i = 0; i < NN; i++)
		{
			Y1[i] = Y1s[i + 1];
			XX[i] = XXs[i + 1];
		}
	}
	
	double B0teor, I10_old;
	I10 = 0.0; Y1_max = -1.0e+100; Y1_min = 1.0e+100; I10_old = 0.0;
	B = XX[0]; Y1_old = Y1[0];
	for (i = 1; i < NN; i++)
	{
		B = XX[i];
		/*			if (B <= B0)
						I10 += 0.5*(Y1[i] + Y1_old)*dB;*/
		I10 += 0.5*(Y1[i] + Y1_old)*(XX[i] - XX[i - 1]);
		if (I10 > I10_old) { B0teor = B; I10_old = I10; }
		//fprintf(fp, "\n%le\t%le\t%le", B, Y1[i], XX[i]);
		if (Y1[i] > Y1_max) { Y1_max = Y1[i];  i1_max = i; }
		if (Y1[i] < Y1_min) { Y1_min = Y1[i];  i1_min = i; }
		Y1_old = Y1[i];
	}
	B0 = B0teor;
	I10 = I10_old;
	for (i = 0; i < NN; i++)
	{
		Xr[i] = XX[i] - B0;
		//fprintf(fp2, "\n %le\t %le\t%le\t  ", XX[i], Y1[i]);
	}
	/*
	if (smoothing != 0) {
		fp2 = fopen("test6asa55-B,x,Y1_spectra_smooth.dat", "w");
		for (i = 1; i < NN; i++) { fprintf(fp2, "\n%le\t%le", XX[i], Y1[i]); }
		fclose(fp2);
	}
	*/
	int *NsimM = NULL, *i1_scanM = NULL;;
	NsimM = _cal_1dint(NsimM, 200);
	i1_scanM = _cal_1dint(i1_scanM, 200);
	int t = 0, m = 0;
	NsimM[0] = 3; i1_scanM[0] = 0;

	while (NsimM[t] < 180)
	{

			Nsim = NsimM[t];
			i1_scan = 4;
			if (Nsim > 0)
			{

				err_parab_min = 1.0e+100;
				j = i1_max - i1_scan; j1_min = i1_max;
				while (j < i1_max + i1_scan) {

					is_min = j - Nsim; is_max = j + Nsim;
					f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; f7 = 0.0; f8 = 0.0;
					d0 = 0.0; d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
					for (i = is_min; i <= is_max; i++)
					{
						tmp2 = Xr[i] * Xr[i]; tmp3 = tmp2 * Xr[i]; tmp4 = tmp2 * tmp2;
						f1 += Xr[i]; f2 += tmp2; f3 += tmp3; f4 += tmp4;
						f5 += tmp4 * Xr[i]; f6 += tmp4 * tmp2; f7 += tmp3 * tmp4; f8 += tmp4 * tmp4;
						d0 += Y1[i]; d1 += Y1[i] * Xr[i]; d2 += Y1[i] * tmp2; d3 += Y1[i] * tmp3; d4 += Y1[i] * tmp4;
					}
					f1 /= (2 * Nsim + 1); f2 /= (2 * Nsim + 1); f3 /= (2 * Nsim + 1); f4 /= (2 * Nsim + 1);
					f5 /= (2 * Nsim + 1); f6 /= (2 * Nsim + 1); f7 /= (2 * Nsim + 1); f8 /= (2 * Nsim + 1);
					d0 /= (2 * Nsim + 1); d1 /= (2 * Nsim + 1); d2 /= (2 * Nsim + 1); d3 /= (2 * Nsim + 1); d4 /= (2 * Nsim + 1);

					matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
					matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
					matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
					matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
					matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

					deter = determinant55(matr);

					matr[0][0] = d0; matr[1][0] = d1; matr[2][0] = d2; matr[3][0] = d3; matr[4][0] = d4;
					deter0 = determinant55(matr);
					matr[0][0] = 1.0; matr[1][0] = f1; matr[2][0] = f2; matr[3][0] = f3; matr[4][0] = f4;

					matr[0][1] = d0; matr[1][1] = d1; matr[2][1] = d2; matr[3][1] = d3; matr[4][1] = d4;
					deter1 = determinant55(matr);
					matr[0][1] = f1; matr[1][1] = f2; matr[2][1] = f3; matr[3][1] = f4; matr[4][1] = f5;

					matr[0][2] = d0; matr[1][2] = d1; matr[2][2] = d2; matr[3][2] = d3; matr[4][2] = d4;
					deter2 = determinant55(matr);
					matr[0][2] = f2; matr[1][2] = f3; matr[2][2] = f4; matr[3][2] = f5; matr[4][2] = f6;

					matr[0][3] = d0; matr[1][3] = d1; matr[2][3] = d2; matr[3][3] = d3; matr[4][3] = d4;
					deter3 = determinant55(matr);
					matr[0][3] = f3; matr[1][3] = f4; matr[2][3] = f5; matr[3][3] = f6; matr[4][3] = f7;

					matr[0][4] = d0; matr[1][4] = d1; matr[2][4] = d2; matr[3][4] = d3; matr[4][4] = d4;
					deter4 = determinant55(matr);
					matr[0][4] = f4; matr[1][4] = f5; matr[2][4] = f6; matr[3][4] = f7; matr[4][4] = f8;

					c0 = deter0 / deter; c1 = deter1 / deter; c2 = deter2 / deter; c3 = deter3 / deter; c4 = deter4 / deter;



					matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
					matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
					matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
					matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
					matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

					Csolve = solve(d0, d1, d2, d3, d4, matr);
					c0 = Csolve[0]; c1 = Csolve[1]; c2 = Csolve[2]; c3 = Csolve[3]; c4 = Csolve[4];

					err_parab = 0.0;
					for (i = is_min; i <= is_max; i++)
					{
						tmp = Y1[i] - c0 - c1 * Xr[i] - c2 * Xr[i] * Xr[i] - c3 * Xr[i] * Xr[i] * Xr[i] - c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i];
						err_parab += tmp * tmp;
					}
					if (err_parab < err_parab_min) { err_parab_min = err_parab; j1_min = j; }
					j++;
				} /*end j*/

				is_min = j1_min - Nsim; is_max = j1_min + Nsim;
				f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; f7 = 0.0; f8 = 0.0;
				d0 = 0.0; d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
				for (i = is_min; i <= is_max; i++)
				{
					tmp2 = Xr[i] * Xr[i]; tmp3 = tmp2 * Xr[i]; tmp4 = tmp2 * tmp2;
					f1 += Xr[i]; f2 += tmp2; f3 += tmp3; f4 += tmp4;
					f5 += tmp4 * Xr[i]; f6 += tmp4 * tmp2; f7 += tmp3 * tmp4; f8 += tmp4 * tmp4;
					d0 += Y1[i]; d1 += Y1[i] * Xr[i]; d2 += Y1[i] * tmp2; d3 += Y1[i] * tmp3; d4 += Y1[i] * tmp4;
				}
				f1 /= (2 * Nsim + 1); f2 /= (2 * Nsim + 1); f3 /= (2 * Nsim + 1); f4 /= (2 * Nsim + 1);
				f5 /= (2 * Nsim + 1); f6 /= (2 * Nsim + 1); f7 /= (2 * Nsim + 1); f8 /= (2 * Nsim + 1);
				d0 /= (2 * Nsim + 1); d1 /= (2 * Nsim + 1); d2 /= (2 * Nsim + 1); d3 /= (2 * Nsim + 1); d4 /= (2 * Nsim + 1);

				matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
				matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
				matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
				matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
				matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

				deter = determinant55(matr);

				matr[0][0] = d0; matr[1][0] = d1; matr[2][0] = d2; matr[3][0] = d3; matr[4][0] = d4;
				deter0 = determinant55(matr);
				matr[0][0] = 1.0; matr[1][0] = f1; matr[2][0] = f2; matr[3][0] = f3; matr[4][0] = f4;

				matr[0][1] = d0; matr[1][1] = d1; matr[2][1] = d2; matr[3][1] = d3; matr[4][1] = d4;
				deter1 = determinant55(matr);
				matr[0][1] = f1; matr[1][1] = f2; matr[2][1] = f3; matr[3][1] = f4; matr[4][1] = f5;

				matr[0][2] = d0; matr[1][2] = d1; matr[2][2] = d2; matr[3][2] = d3; matr[4][2] = d4;
				deter2 = determinant55(matr);
				matr[0][2] = f2; matr[1][2] = f3; matr[2][2] = f4; matr[3][2] = f5; matr[4][2] = f6;

				matr[0][3] = d0; matr[1][3] = d1; matr[2][3] = d2; matr[3][3] = d3; matr[4][3] = d4;
				deter3 = determinant55(matr);
				matr[0][3] = f3; matr[1][3] = f4; matr[2][3] = f5; matr[3][3] = f6; matr[4][3] = f7;

				matr[0][4] = d0; matr[1][4] = d1; matr[2][4] = d2; matr[3][4] = d3; matr[4][4] = d4;
				deter4 = determinant55(matr);
				matr[0][4] = f4; matr[1][4] = f5; matr[2][4] = f6; matr[3][4] = f7; matr[4][4] = f8;

				c0 = deter0 / deter; c1 = deter1 / deter; c2 = deter2 / deter; c3 = deter3 / deter; c4 = deter4 / deter;

				/*
				error_Xr=1.0; check=0; Xr_tmp=XX[i1_max]-B0;
				while (fabs(error_Xr)>1.0e-12 && check<1000)
				{
				Xr_max=-(c1+3.0*c3*Xr_tmp*Xr_tmp+4.0*c4*Xr_tmp*Xr_tmp*Xr_tmp)/(2.0*c2);
				error_Xr=Xr_max-Xr_tmp; Xr_tmp=Xr_max; check++;
				}
				*/
				matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
				matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
				matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
				matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
				matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

				Csolve = solve(d0, d1, d2, d3, d4, matr);
				c0 = Csolve[0]; c1 = Csolve[1]; c2 = Csolve[2]; c3 = Csolve[3]; c4 = Csolve[4];

				Xr_max = root_cube1(0.75*c3 / c4, 0.5*c2 / c4, 0.25*c1 / c4, XX[i1_max] - B0);
				x1_max_s = Xr_max;
				Y1_max_s = c0 + c1 * Xr_max + c2 * Xr_max*Xr_max + c3 * Xr_max*Xr_max*Xr_max + c4 * Xr_max*Xr_max*Xr_max*Xr_max;



				/*	fp2 = fopen("tzall6asa55_parabol_max.dat", "w");
					fprintf(fp2, "\n i1_max=%d j1_min=%d err_parab_min=%le Xr_max=%le x1_max_s=%le Y1_max=%le",
						i1_max, j1_min, err_parab_min, Xr_max, x1_max_s, Y1_max_s);
					fprintf(fp2, "\n deter=%le deter0=%le deter1=%le deter2=%le deter3=%le deter4=%le", deter, deter0, deter1, deter2, deter3, deter4);
					fprintf(fp2, "\n c0=%le c1=%le c2=%le c3=%le c4=%le", c0, c1, c2, c3, c4);
					fprintf(fp2, "\n f1=%le f2=%le f3=%le f4=%le f5=%le f6=%le f7=%le f8=%le", f1, f2, f3, f4, f5, f6, f7, f8);
					fprintf(fp2, "\n d0=%le d1=%le d2=%le d3=%le d4=%le", d0, d1, d2, d3, d4);

					for (i = j1_min - Nsim; i <= j1_min + Nsim; i++)
					{
						fprintf(fp2, "\n%le\t%le\t%le", XX[i], c0 + c1 * Xr[i] + c2 * Xr[i] * Xr[i] + c3 * Xr[i] * Xr[i] * Xr[i] + c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i], Y1[i]);
					}
					fclose(fp2);


					fp2 = fopen("tzall6asa55_parabol_max_full.dat", "w");
					fprintf(fp2, "\n i1_max=%d j1_min=%d err_parab_min=%le Xr_max=%le  x1_max_s=%le Y1_max=%le",
						i1_max, j1_min, err_parab_min, Xr_max,  x1_max_s, Y1_max_s);
					fprintf(fp2, "\n deter=%le deter0=%le deter1=%le deter2=%le deter3=%le deter4=%le", deter, deter0, deter1, deter2, deter3, deter4);
					fprintf(fp2, "\n c0=%le c1=%le c2=%le c3=%le c4=%le", c0, c1, c2, c3, c4);
					fprintf(fp2, "\n f1=%le f2=%le f3=%le f4=%le f5=%le f6=%le f7=%le f8=%le", f1, f2, f3, f4, f5, f6, f7, f8);
					fprintf(fp2, "\n d0=%le d1=%le d2=%le d3=%le d4=%le", d0, d1, d2, d3, d4);

					for (i = 0; i < NN; i++)
					{
						fprintf(fp2, "\n%le\t%le\t%le", XX[i], c0 + c1 * Xr[i] + c2 * Xr[i] * Xr[i] + c3 * Xr[i] * Xr[i] * Xr[i] + c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i]);
					}
					fclose(fp2); */





				err_parab_min = 1.0e+100;
				j = i1_min - i1_scan; j1_min = i1_min;
				while (j < i1_min + i1_scan) {
					is_min = j - Nsim; is_max = j + Nsim;

					f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; f7 = 0.0; f8 = 0.0;
					d0 = 0.0; d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
					for (i = is_min; i <= is_max; i++)
					{
						tmp2 = Xr[i] * Xr[i]; tmp3 = tmp2 * Xr[i]; tmp4 = tmp2 * tmp2;
						f1 += Xr[i]; f2 += tmp2; f3 += tmp3; f4 += tmp4;
						f5 += tmp4 * Xr[i]; f6 += tmp4 * tmp2; f7 += tmp3 * tmp4; f8 += tmp4 * tmp4;
						d0 += Y1[i]; d1 += Y1[i] * Xr[i]; d2 += Y1[i] * tmp2; d3 += Y1[i] * tmp3; d4 += Y1[i] * tmp4;
					}
					f1 /= (2 * Nsim + 1); f2 /= (2 * Nsim + 1); f3 /= (2 * Nsim + 1); f4 /= (2 * Nsim + 1);
					f5 /= (2 * Nsim + 1); f6 /= (2 * Nsim + 1); f7 /= (2 * Nsim + 1); f8 /= (2 * Nsim + 1);
					d0 /= (2 * Nsim + 1); d1 /= (2 * Nsim + 1); d2 /= (2 * Nsim + 1); d3 /= (2 * Nsim + 1); d4 /= (2 * Nsim + 1);

					matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
					matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
					matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
					matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
					matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

					deter = determinant55(matr);

					matr[0][0] = d0; matr[1][0] = d1; matr[2][0] = d2; matr[3][0] = d3; matr[4][0] = d4;
					deter0 = determinant55(matr);
					matr[0][0] = 1.0; matr[1][0] = f1; matr[2][0] = f2; matr[3][0] = f3; matr[4][0] = f4;

					matr[0][1] = d0; matr[1][1] = d1; matr[2][1] = d2; matr[3][1] = d3; matr[4][1] = d4;
					deter1 = determinant55(matr);
					matr[0][1] = f1; matr[1][1] = f2; matr[2][1] = f3; matr[3][1] = f4; matr[4][1] = f5;

					matr[0][2] = d0; matr[1][2] = d1; matr[2][2] = d2; matr[3][2] = d3; matr[4][2] = d4;
					deter2 = determinant55(matr);
					matr[0][2] = f2; matr[1][2] = f3; matr[2][2] = f4; matr[3][2] = f5; matr[4][2] = f6;

					matr[0][3] = d0; matr[1][3] = d1; matr[2][3] = d2; matr[3][3] = d3; matr[4][3] = d4;
					deter3 = determinant55(matr);
					matr[0][3] = f3; matr[1][3] = f4; matr[2][3] = f5; matr[3][3] = f6; matr[4][3] = f7;

					matr[0][4] = d0; matr[1][4] = d1; matr[2][4] = d2; matr[3][4] = d3; matr[4][4] = d4;
					deter4 = determinant55(matr);
					matr[0][4] = f4; matr[1][4] = f5; matr[2][4] = f6; matr[3][4] = f7; matr[4][4] = f8;

					c0 = deter0 / deter; c1 = deter1 / deter; c2 = deter2 / deter; c3 = deter3 / deter; c4 = deter4 / deter;

					matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
					matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
					matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
					matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
					matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

					Csolve = solve(d0, d1, d2, d3, d4, matr);
					c0 = Csolve[0]; c1 = Csolve[1]; c2 = Csolve[2]; c3 = Csolve[3]; c4 = Csolve[4];

					err_parab = 0.0;
					for (i = is_min; i <= is_max; i++)
					{
						tmp = Y1[i] - c0 - c1 * Xr[i] - c2 * Xr[i] * Xr[i] - c3 * Xr[i] * Xr[i] * Xr[i] - c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i];
						err_parab += tmp * tmp;
					}
					if (err_parab < err_parab_min) { err_parab_min = err_parab; j1_min = j; }
					j++;
				}


				is_min = j1_min - Nsim; is_max = j1_min + Nsim;
				f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; f7 = 0.0; f8 = 0.0;
				d0 = 0.0; d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
				for (i = is_min; i <= is_max; i++)
				{
					tmp2 = Xr[i] * Xr[i]; tmp3 = tmp2 * Xr[i]; tmp4 = tmp2 * tmp2;
					f1 += Xr[i]; f2 += tmp2; f3 += tmp3; f4 += tmp4;
					f5 += tmp4 * Xr[i]; f6 += tmp4 * tmp2; f7 += tmp3 * tmp4; f8 += tmp4 * tmp4;
					d0 += Y1[i]; d1 += Y1[i] * Xr[i]; d2 += Y1[i] * tmp2; d3 += Y1[i] * tmp3; d4 += Y1[i] * tmp4;
				}
				f1 /= (2 * Nsim + 1); f2 /= (2 * Nsim + 1); f3 /= (2 * Nsim + 1); f4 /= (2 * Nsim + 1);
				f5 /= (2 * Nsim + 1); f6 /= (2 * Nsim + 1); f7 /= (2 * Nsim + 1); f8 /= (2 * Nsim + 1);
				d0 /= (2 * Nsim + 1); d1 /= (2 * Nsim + 1); d2 /= (2 * Nsim + 1); d3 /= (2 * Nsim + 1); d4 /= (2 * Nsim + 1);

				matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
				matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
				matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
				matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
				matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

				deter = determinant55(matr);

				matr[0][0] = d0; matr[1][0] = d1; matr[2][0] = d2; matr[3][0] = d3; matr[4][0] = d4;
				deter0 = determinant55(matr);
				matr[0][0] = 1.0; matr[1][0] = f1; matr[2][0] = f2; matr[3][0] = f3; matr[4][0] = f4;

				matr[0][1] = d0; matr[1][1] = d1; matr[2][1] = d2; matr[3][1] = d3; matr[4][1] = d4;
				deter1 = determinant55(matr);
				matr[0][1] = f1; matr[1][1] = f2; matr[2][1] = f3; matr[3][1] = f4; matr[4][1] = f5;

				matr[0][2] = d0; matr[1][2] = d1; matr[2][2] = d2; matr[3][2] = d3; matr[4][2] = d4;
				deter2 = determinant55(matr);
				matr[0][2] = f2; matr[1][2] = f3; matr[2][2] = f4; matr[3][2] = f5; matr[4][2] = f6;

				matr[0][3] = d0; matr[1][3] = d1; matr[2][3] = d2; matr[3][3] = d3; matr[4][3] = d4;
				deter3 = determinant55(matr);
				matr[0][3] = f3; matr[1][3] = f4; matr[2][3] = f5; matr[3][3] = f6; matr[4][3] = f7;

				matr[0][4] = d0; matr[1][4] = d1; matr[2][4] = d2; matr[3][4] = d3; matr[4][4] = d4;
				deter4 = determinant55(matr);
				matr[0][4] = f4; matr[1][4] = f5; matr[2][4] = f6; matr[3][4] = f7; matr[4][4] = f8;

				c0 = deter0 / deter; c1 = deter1 / deter; c2 = deter2 / deter; c3 = deter3 / deter; c4 = deter4 / deter;

				matr[0][0] = 1.0; matr[0][1] = f1; matr[0][2] = f2; matr[0][3] = f3; matr[0][4] = f4;
				matr[1][0] = f1; matr[1][1] = f2; matr[1][2] = f3; matr[1][3] = f4; matr[1][4] = f5;
				matr[2][0] = f2; matr[2][1] = f3; matr[2][2] = f4; matr[2][3] = f5; matr[2][4] = f6;
				matr[3][0] = f3; matr[3][1] = f4; matr[3][2] = f5; matr[3][3] = f6; matr[3][4] = f7;
				matr[4][0] = f4; matr[4][1] = f5; matr[4][2] = f6; matr[4][3] = f7; matr[4][4] = f8;

				Csolve = solve(d0, d1, d2, d3, d4, matr);
				c0 = Csolve[0]; c1 = Csolve[1]; c2 = Csolve[2]; c3 = Csolve[3]; c4 = Csolve[4];

				/*
				error_Xr=1.0; check=0; Xr_tmp=XX[i1_min]-B0;
				while (fabs(error_Xr)>1.0e-12 && check<1000)
				{
				Xr_min=-(c1+3.0*c3*Xr_tmp*Xr_tmp+4.0*c4*Xr_tmp*Xr_tmp*Xr_tmp)/(2.0*c2);
				error_Xr=Xr_min-Xr_tmp; Xr_tmp=Xr_min; check++;
				}
				*/

				Xr_min = root_cube1(0.75*c3 / c4, 0.5*c2 / c4, 0.25*c1 / c4, XX[i1_min] - B0);
				x1_min_s = Xr_min;
				Y1_min_s = c0 + c1 * Xr_min + c2 * Xr_min*Xr_min + c3 * Xr_min*Xr_min*Xr_min + c4 * Xr_min*Xr_min*Xr_min*Xr_min;
			if (t == 40) ///////////////////////////

				{
					for (i = j1_min - Nsim; i <= j1_min + Nsim; i++)
					{
							fprintf(fp1, "\n%le\t%le\t%le", XX[i], c0 + c1 * Xr[i] + c2 * Xr[i] * Xr[i] + c3 * Xr[i] * Xr[i] * Xr[i] + c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i], Y1[i]);
					}
					fclose(fp1);
				}  ////////////////////////////
				/*fp2 = fopen("tzall6asa55_parabol_min.dat", "w");
				fprintf(fp2, "\n i1_max=%d j1_min=%d err_parab_min=%le Xr_min=%le x1_min_s=%le Y1_max=%le",
					i1_min, j1_min, err_parab_min, Xr_min, x1_min_s, Y1_max_s);
				fprintf(fp2, "\n deter=%le deter0=%le deter1=%le deter2=%le deter3=%le deter4=%le", deter, deter0, deter1, deter2, deter3, deter4);
				fprintf(fp2, "\n c0=%le c1=%le c2=%le c3=%le c4=%le", c0, c1, c2, c3, c4);
				fprintf(fp2, "\n f1=%le f2=%le f3=%le f4=%le f5=%le f6=%le f7=%le f8=%le", f1, f2, f3, f4, f5, f6, f7, f8);
				fprintf(fp2, "\n d0=%le d1=%le d2=%le d3=%le d4=%le", d0, d1, d2, d3, d4);

				for (i = j1_min - Nsim; i <= j1_min + Nsim; i++)
				{
					fprintf(fp2, "\n%le\t%le\t%le", XX[i], c0 + c1 * Xr[i] + c2 * Xr[i] * Xr[i] + c3 * Xr[i] * Xr[i] * Xr[i] + c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i], Y1[i]);
				}
				fclose(fp2);


				fp2 = fopen("tzall6asa55_parabol_min_full.dat", "w");
				fprintf(fp2, "\n i1_max=%d j1_min=%d err_parab_min=%le Xr_min=%le x1_min_s=%le Y1_max=%le",
					i1_min, j1_min, err_parab_min, Xr_min,  x1_min_s, Y1_max_s);
				fprintf(fp2, "\n deter=%le deter0=%le deter1=%le deter2=%le deter3=%le deter4=%le", deter, deter0, deter1, deter2, deter3, deter4);
				fprintf(fp2, "\n c0=%le c1=%le c2=%le c3=%le c4=%le", c0, c1, c2, c3, c4);
				fprintf(fp2, "\n f1=%le f2=%le f3=%le f4=%le f5=%le f6=%le f7=%le f8=%le", f1, f2, f3, f4, f5, f6, f7, f8);
				fprintf(fp2, "\n d0=%le d1=%le d2=%le d3=%le d4=%le", d0, d1, d2, d3, d4);

				for (i = 0; i < NN; i++)
				{
					fprintf(fp2, "\n%le\t%le", XX[i], c0 + c1 * Xr[i] + c2 * Xr[i] * Xr[i] + c3 * Xr[i] * Xr[i] * Xr[i] + c4 * Xr[i] * Xr[i] * Xr[i] * Xr[i]);
				}
				fclose(fp2); */


				App = Y1_max_s - Y1_min_s; delta_Bpp = (x1_min_s - x1_max_s);

			}
			B0 = (Xr_min + Xr_max) / 2.0 + B0;
			if (t == 0)
			{
				int l = 0;
				int left = 0, right = 0;
				i = 0;
				if ((2 * B0 - XX[0]) < XX[NN-1])
				{
					while (XX[i] <= (2 * B0 - XX[0]))
					{
						l = i;
						i++;
					}
					left = 0; right = l;
				}
				else
				{		while (i<NN)
						{
							if ((XX[i] >= (B0 - (XX[NN - 1] - B0))) && (XX[i] <= XX[NN - 1]))
							{
								l++;
							}
							i++;
						}
					left = NN - 1 - l; right = NN-1;
				}
				double dY=0, dY1=0,Ysum1=0, Ysum2 = 0;
				double deltaY=0;
				for (i = 0; i < (right-left)/10; i++)
				{
					Ysum1 = Ysum1 + Y1[left+i];
					Ysum2 = Ysum2 + Y1[right-i];
				}
				dY = (Ysum1 - Ysum2) / (double)((right - left) / 10);
				dY1 = (Ysum1)/ (double)((right - left) / 10);
				deltaY = dY / 2.0 - dY1;
				//NN = l;
				for (i = 0; i < NN; i++)
				{
					Y1[i] = Y1[i] + deltaY;
					if (flag == 0)
					Y1ish[i] = Y1ish[i] + deltaY;
				}
				if (flag == 1)
				{
					for (i = 0; i < NN*3; i++)
					{
						Y1ish[i] = Y1ish[i] + deltaY;
					}
				}
			}

			I10 = 0.0;
			B = XX[0]; Y1_old = Y1[0];
			for (i = 1; i < NN; i++)
			{
				B = XX[i];
				if (B <= B0)
					I10 += 0.5*(Y1[i] + Y1_old)*(XX[i]-XX[i-1]);
				Y1_old = Y1[i];
			}
			//I10 = 898.0;
			//B0 = 3247.391;
			/*
			while (fabs(error_q) > 1.0e-12) {
				dq = -1.0 + 8.0*I10*pow(2.0*q / (q + 1.0), a - 1.0) / (App*delta_Bpp*(1.0 - pow(1.0 + gamma * gamma*(q - 1.0) / (q + 1.0), a)));
				error_q = q - dq; q = dq; a = -1.0 / (q - 1.0); check++;
			}
			q_error_min = q; error_q_min = error_q;
			*/
			double a1, q1 = 1.000001, q3 = 1.000001, fa;
			double a2, q2 = 2.999999, fb;
			int ii = 0;
			while (ii < 4000)

			{

				a1 = -1.0 / (q1 - 1.0);
				a2 = -1.0 / (q2 - 1.0);

				gamma = 2.0*(XX[0] - B0) / delta_Bpp;
				//gamma = 2.0*(X[0] - Xmax) / dxpp(qq[k]);
				fa = (0.125*App / I10 - pow(2.0*q1 / (q1 + 1.0), a1 - 1.0) / ((q1 + 1.0)*delta_Bpp*(1.0 - pow(1.0 + gamma * gamma*(q1 - 1.0) / (q1 + 1.0), a1))));
				fb = (0.125*App / I10 - pow(2.0*q2 / (q2 + 1.0), a2 - 1.0) / ((q2 + 1.0)*delta_Bpp*(1.0 - pow(1.0 + gamma * gamma*(q2 - 1.0) / (q2 + 1.0), a2))));
				if (fa*fb < 0)
				{
					q3 = q1;
					q1 = (q1 + q2) / 2.0;
				}
				else
				{
					q2 = q1;
					q1 = q3;
				}
				ii++;
			}
			q_error_min = q1;
			q = q1;

			a = -1.0 / (q - 1.0); b = pow(2.0, q - 1.0) - 1.0;
			x0_min = -1.0 / sqrt(b*(1.0 - 2.0*a));
			G_min = -delta_Bpp / (2.0*x0_min);
			x1 = (XX[0] - B0) / G_min;
			Ymax_min = I10 / (1.0 - pow(1.0 + b * x1*x1, a));



			double *Y1T = NULL;
			double error=0.0;
			if (flag == 1) NN = NN * 3;
			Y1T = _cal_1d(Y1T, NN+1);
			for (i = 0; i < NN ; i++)
			{
				x = (XXish[i] - B0) / G_min;
				Y1T[i] = (Ymax_min / G_min)*(2.0*a*b*x)*pow(1.0 + b * x*x, a - 1.0);
				error = error + (Y1ish[i] - Y1T[i])*(Y1ish[i] - Y1T[i]);

			   if (Nsim == 4)
				fprintf(fp2, "\n %le\t %le\t%le\t%le\t%le\t  ",XXish[i], Y1ish[i], Y1T[i],XX[i],Y1[i]);
			}
	
			if (t == 0) { errorlast = error; qopt = q; G_opt = G_min; Ymax_opt = Ymax_min; B0_opt = B0; }
			if (error < errorlast) { errorlast = error; qopt = q; G_opt = G_min; Nsimopt = Nsim; Ymax_opt = Ymax_min; B0_opt = B0;}
			free(Y1T);
			B0 = B0teor;
		errorQ[t] = error;
		QN[t] = q;
		GN[t] = G_min;
		NsimM[t + 1] = NsimM[t] + 1;
		t = t + 1;
		fprintf(fp3, "\n %le\t %le\t%le\t %d\t  ", q, pow(error, 0.5) / (double)NN, G_min, Nsim);
		if (flag == 1) NN = NN / 3;
	}
	double QNmin = 5.0,QNmax = 0.0;
	double GNmin = 100.0, GNmax = 0.0;
	for (i = 0; i < t; i++)
	{
		if (fabs(errorQ[i]- errorlast) < (errorlast*0.3))
		{
			if (QN[i] < QNmin) { QNmin = QN[i]; }
			if (GN[i] < GNmin) { GNmin = GN[i]; }
			if (QN[i] > QNmax) { QNmax = QN[i]; }
			if (GN[i] > GNmax) { GNmax = GN[i]; }
		}
	}
	if ((fabs(QNmax - qopt)) > (fabs(QNmin - qopt)))
	{
		dq = fabs(QNmax - qopt);
	}
	else
	{
		dq = fabs(QNmin - qopt);
	}
	if ((fabs(GNmax - G_opt)) > (fabs(GNmin - G_opt)))
	{
		dG = fabs(GNmax - G_opt);
	}
	else
	{
		dG = fabs(GNmin - G_opt);
	}
	if (flag == 1) NN = NN * 3;
	fprintf(fp4, "\n %le\t %le\t %le\t %le\t %le\t  %d\t  ", qopt, G_opt, dq, dG, pow(errorlast, 0.5) / (double)NN, Nsimopt);
	
	FILE *fp10 = NULL;
	fp10 = fopen("DIPLOM1.dat", "w");
	double *Y1T = NULL;
	double error = 0.0;
	Y1T = _cal_1d(Y1T, NN + 1);
	a = -1.0 / (qopt - 1.0); b = pow(2.0, qopt - 1.0) - 1.0;
	for (i = 0; i < NN; i++)
	{
		x = (XXish[i] - B0_opt) / G_opt;
		Y1T[i] = (Ymax_opt / G_opt)*(2.0*a*b*x)*pow(1.0 + b * x*x, a - 1.0);
		error = error + (Y1ish[i] - Y1T[i])*(Y1ish[i] - Y1T[i]);
		fprintf(fp10, "\n %le\t %le\t%le\t", XXish[i], Y1ish[i], Y1T[i]);
		if (i == NN - 1)
		{
		fprintf(fp10, "\n %le\t", pow(error, 0.5) / (double)NN );
		}

	}
	fclose(fp10);

	fclose(fp3);
	fclose(fp);
	fclose(fp4);


}