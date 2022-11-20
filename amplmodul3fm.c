#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <process.h>
#include <time.h>



FILE *fp2 = NULL;
int NN;
double tmp, tmp1, tmp2, tmp3, tmp4, tmp5;
double *a1=NULL, *a1t = NULL, *Hdm = NULL, f_max, f_min, dH, H0, Hd, G, G2, q, Hmax, Hmin;
double err, err_min, err_max, qt, Gt, qt_min, Gt_min, dHpp_t;
int Kscan, i_err_max, i_err_min;


#define pi  3.1415926535897932384626433832795
#define pi2  6.283185307179586476925286766559             
#define odnashestaja  0.16666666666666666666666666666667
#define enat 2.71828182845904523536028747135266249775724709369995 
#define ln2  0.69314718055994530941723212145818


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;};
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXIT 60
#define UNUSED (-1.11e30)



double* _cal_1d(double *a, int N)
{
	if ((a = (double*)calloc(N, sizeof(double))) == NULL)
	{
		printf("Memory is not sufficient");
		exit(-1);
	}
	return(a);
}

int* _cal_1i(int *a, int N)
{
	if ((a = (int*)calloc(N, sizeof(int))) == NULL)
	{
		printf("Memory is not sufficient");
		exit(-1);
	}
	return(a);
}


/* logarithm of gamma-function */
double gammLn(double x)
{
	double y, ser, *co;
	int j;
	static double cof[8] = {
	   2.5066282746310005,
	   1.0000000000190015,
	   76.18009172947146,
	  -86.50532032941677,
	   24.01409824083091,
	  -1.231739572450155,
	   0.1208650973866179e-2,
	  -0.5395239384953e-5, };

	/* calculate the series */
	ser = cof[1]; y = x; co = cof + 2;
	for (j = 2; j < 8; j++) {
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
	return(exp(gammLn(x)));
}

/* the beta function */
double Beta_LM(double x, double y)
{
	return(exp(gammLn(x) + gammLn(y) - gammLn(x + y)));
}




void qG_finder(double qt_start, double qt_fin, double Gt_start, double Gt_fin, double d_qt, double d_Gt)
{
	int i;
	double err0;

	err_min = 1.0e+100; err_max = -1.0e+100;
	qt = qt_start;
	while (qt < qt_fin)
	{
		Gt = Gt_start;
		while (Gt < Gt_fin)
		{
			tmp1 = pow(2.0, qt - 1.0) - 1.0;
			Hd = -Kscan * G; i = 0; f_max = -1.0e+10; f_min = 1.0e+10;
			while (Hd <= Kscan * G)
			{
				tmp2 = Hd / Gt;
				a1t[i] = -Hd * pow(1.0 + tmp1 * tmp2*tmp2, -qt / (qt - 1.0));
				if (a1t[i] > f_max) { f_max = a1t[i]; Hmax = Hd; }
				if (a1t[i] < f_min) { f_min = a1t[i]; Hmin = Hd; }
				Hd += dH; i++;
			}
			for (i = 0; i < NN; i++) { a1t[i] /= f_max; }

			err = 0.0;
			for (i = 0; i < NN; i++) { err += (a1t[i] - a1[i])*(a1t[i] - a1[i]); }
			if (err < err_min) { err_min = err; Gt_min = Gt; qt_min = qt; dHpp_t = Hmin - Hmax; }
			Gt += d_Gt;
		}
		qt += d_qt;
	}
	err = sqrt(err_min) / NN;

	tmp1 = pow(2.0, qt_min - 1.0) - 1.0;
	Hd = -Kscan * G; i = 0; f_max = -1.0e+10; f_min = 1.0e+10;
	while (Hd <= Kscan * G)
	{
		tmp2 = Hd / Gt_min;
		a1t[i] = -Hd * pow(1.0 + tmp1 * tmp2*tmp2, -qt_min / (qt_min - 1.0));
		if (a1t[i] > f_max) { f_max = a1t[i]; Hmax = Hd; }
		if (a1t[i] < f_min) { f_min = a1t[i]; Hmin = Hd; }
		Hd += dH; i++;
	}
	for (i = 0; i < NN; i++) { a1t[i] /= f_max; }

	err_max = -1.0e+10; err_min = 1.0e+10;
	for (i = 0; i < NN; i++)
	{
		err0 = (a1t[i] - a1[i])*(a1t[i] - a1[i]);
		/*err0/=(a1t[i]*a1t[i]);*/
		if (err0 > err_max) { err_max = err0; i_err_max = i; }
		if (err0 < err_min) { err_min = err0; i_err_min = i; }
	}

	err_max = sqrt(err_max);
	err_min = sqrt(err_min);

}



void main()
{
	int i, j;
	FILE *fp=NULL, *fp1 = NULL;
	double sum;
	double teta, dHpp, App, hm, Ha, d_teta, f, fmax, q, beta, beta0;
	double *sn= NULL;
	double q_s, q_f, G_s, G_f, d_q, d_G;
	int z = 0;
	Kscan = 18;

	sn = _cal_1d(sn, 18010);
	a1 = _cal_1d(a1, Kscan * 1010);
	a1t = _cal_1d(a1t, Kscan * 1010);
	Hdm = _cal_1d(Hdm, Kscan * 1010);



	H0 = 3250.0; G = 1.0 ; G2 = G * G; d_teta = 2.0*pi / 10000.0;
	i = 0; while (i <= 10000) { teta = i * d_teta - pi; sn[i] = sin(teta); i++; }

	dH = 2.0*Kscan*G / (Kscan*1000.0);

	q = 1.000001;
	while (q <= 2.6)
	{
		fp = fopen("amplmodul3fm_spectra.dat", "a");
		fprintf(fp, "\n q=%le G=%le\t Kscan=%d", q, G, Kscan);
		fclose(fp);

		fp = fopen("amplmodul3fm_norm_spectra.dat", "a");
		fprintf(fp, "\n q=%le G=%le\t Kscan=%d", q, G, Kscan);
		fclose(fp);

		fp1 = fopen("amplmodul3fm_hm,dHpp,App,qt,Gt,dHppt,err.dat", "a");
		fprintf(fp1, "\n q=%le G=%le\t Kscan=%d", q, G, Kscan);
		fclose(fp1);

		printf("\n q=%le", q);
		double dHpp_real;
		dHpp_real = 2.0 * G*pow(((q - 1) / ((q + 1)*(pow(2, (q - 1)) - 1))), 0.5);
		hm = 0.1;
		while (hm <= 4.0)
		{
			printf("\n hm=%le", hm);

			
			tmp1 = pow(2.0, q - 1.0) - 1.0;
			/*beta=beta0*exp(gammln(1.0/(q-1.0)-0.5))/exp(gammln(1.0/(q-1.0)));*/
			beta = Beta_LM(0.5, 1.0 / (q - 1.0) - 0.5);
			fmax = sqrt(tmp1) / (G*beta);
			//fmax = 1.0; //гюлемх рср еякх врн
			Hd = -Kscan * G; f_max = -1.0e+10; f_min = 1.0e+10; j = 0;
			while (Hd <= Kscan * G)
			{
				i = 0; sum = 0.0;
				while (i <= 10000)
				{
					tmp2 = (Hd + 0.5*hm*sn[i]) / G;
					f = sn[i] * pow(1.0 + tmp1 * tmp2*tmp2, -1.0 / (q - 1.0));
					sum += f;
					i++;
				}
				a1[j] = sum * fmax*d_teta; Hdm[j] = Hd;
				if (a1[j] > f_max) { f_max = a1[j]; Hmax = Hd; }
				if (a1[j] < f_min) { f_min = a1[j]; Hmin = Hd; }
				Hd += dH;  j++;
			}
			NN = j;
			dHpp = Hmin - Hmax; App = f_max - f_min;

			fp = fopen("amplmodul3fm_spectra.dat", "a");
			fprintf(fp, "\n hm=%le ", hm);
			for (i = 0; i < NN; i++) { fprintf(fp, "\n%le\t%le", Hdm[i], a1[i]); }
			fclose(fp);

			for (i = 0; i < NN; i++) { a1[i] /= f_max; }

			q_s = 0.1; q_f = 3.0; G_s = G / 10.0; G_f = 10.0*G; d_q = 0.1; d_G = 0.1;
			for (i = 0; i < 10; i++)
			{
				qG_finder(q_s, q_f, G_s, G_f, d_q, d_G);
				q_s = qt_min - 2.0*d_q; q_f = qt_min + 2.0*d_q; G_s = Gt_min - 2.0*d_G; G_f = Gt_min + 2.0*d_G; d_q /= 2.0; d_G /= 2.0;
			}

			fp1 = fopen("amplmodul3fm_hm,dHpp,App,qt,Gt,dHppt,err.dat", "a");
			fprintf(fp1, "\n%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le",
				hm, dHpp, App, qt_min, Gt_min, dHpp_t, err, err_min, err_max, Hdm[i_err_min] / G, Hdm[i_err_max] / G);
			fclose(fp1);

			fp = fopen("amplmodul3fm_norm_spectra.dat", "a");
			fprintf(fp, "\n hm=%le ", hm);
			for (i = 0; i < NN; i++) { fprintf(fp, "\n%le\t%le\t%le", Hdm[i], a1[i], a1t[i]); }
			fclose(fp);

			fp = fopen("APPOT(HM).dat", "a");
			fprintf(fp, "\n%le\t%le\t%le\t%le\t%le", hm, App, dHpp/dHpp_real, (qt-q)/q, qt);
			fclose(fp);
			FILE *dHpp_q1 = NULL, *dHpp_q15 = NULL, *dHpp_q2 = NULL, *dHpp_q25 = NULL ;
			if (z == 0)
			{
				dHpp_q1 = fopen("dHpp_q1.dat", "a");
				fprintf(dHpp_q1, "\n%le\t%le\t%le\t%le\t%le", hm, App, dHpp / dHpp_real, (qt_min - q) / q, qt_min);
				fclose(dHpp_q1);
			}
			if (z == 1)
			{
				dHpp_q1 = fopen("dHpp_q15.dat", "a");
				fprintf(dHpp_q1, "\n%le\t%le\t%le\t%le\t%le", hm, App, dHpp / dHpp_real, (qt_min - q) / q, qt_min);
				fclose(dHpp_q1);
			}
			if (z == 2)
			{
				dHpp_q1 = fopen("dHpp_q2.dat", "a");
				fprintf(dHpp_q1, "\n%le\t%le\t%le\t%le\t%le", hm, App, dHpp / dHpp_real, (qt_min - q) / q, qt_min);
				fclose(dHpp_q1);
			}
			if (z == 3)
			{
				dHpp_q1 = fopen("dHpp_q25.dat", "a");
				fprintf(dHpp_q1, "\n%le\t%le\t%le\t%le\t%le", hm, App, dHpp / dHpp_real, (qt_min - q) / q, qt_min);
				fclose(dHpp_q1);
			}

			hm += 0.1;

		}
		z += 1;
		q += 0.5;
	}

} /*end main*/
