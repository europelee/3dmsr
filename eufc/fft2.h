/*-------------------------------------------------------------------------
Perform a 2D FFT inplace given a complex 2D array
The direction dir, 1 for forward, -1 for reverse
The size of the array (nx,ny)
Return false if there are memory problems or
the dimensions are not powers of 2
*/
#ifndef FFT2_H
#define FFT2_H
#include"myconfig.h"
int FFT(int dir,int m,float *x,float *y);
int FFT2D(_complex **c,int nx,int ny,int dir)
{
	int i,j;
	int m,twopm;
	float *real,*imag;

	/* Transform the rows */
	real = (float *)malloc(nx * sizeof(float));
	imag = (float *)malloc(nx * sizeof(float));
	if (real == NULL || imag == NULL)
		return(FALSE);
	//if (!Powerof2(nx,&m,&twopm) || twopm != nx)
		//return(FALSE);
	if(nx==128)
		m = 7;
	if(nx==256)
		m = 8;
	for (j=0;j<ny;j++) {
		for (i=0;i<nx;i++) {
			real[i] = c[i][j].x;
			imag[i] = c[i][j].y;
		}
		FFT(dir,m,real,imag);
		for (i=0;i<nx;i++) {
			c[i][j].x = real[i];
			c[i][j].y = imag[i];
		}
	}
	free(real);
	free(imag);

	/* Transform the columns */
	real = (float *)malloc(ny * sizeof(float));
	imag = (float*)malloc(ny * sizeof(float));
	if (real == NULL || imag == NULL)
		return(FALSE);
	//if (!Powerof2(ny,&m,&twopm) || twopm != ny)
		//return(FALSE);
	if(ny==128)
		m = 7;
	if(ny==256)
		m = 8;
	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++) {
			real[j] = c[i][j].x;
			imag[j] = c[i][j].y;
		}
		FFT(dir,m,real,imag);
		for (j=0;j<ny;j++) {
			c[i][j].x = real[j];
			c[i][j].y = imag[j];
		}
	}
	free(real);
	free(imag);

	return(TRUE);
}

/*-------------------------------------------------------------------------
This computes an in-place complex-to-complex FFT
x and y are the real and imaginary arrays of 2^m points.
dir =  1 gives forward transform
dir = -1 gives reverse transform

Formula: forward
N-1
---
1   \          - j k 2 pi n / N
X(n) = ---   >   x(k) e                    = forward transform
N   /                                n=0..N-1
---
k=0

Formula: reverse
N-1
---
\          j k 2 pi n / N
X(n) =       >   x(k) e                    = forward transform
/                                n=0..N-1
---
k=0
*/
int FFT(int dir,int m,float *x,float *y)
{
	long nn,i,i1,j,k,i2,l,l1,l2;
	float c1,c2,tx,ty,t1,t2,u1,u2,z;

	/* Calculate the number of points */
	nn = 1;
	for (i=0;i<m;i++)
		nn *= 2;

	/* Do the bit reversal */
	i2 = nn >> 1;
	j = 0;
	for (i=0;i<nn-1;i++) {
		if (i < j) {
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<m;l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j=0;j<l1;j++) {
			for (i=j;i<nn;i+=l2) {
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = x[i] - t1;
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for forward transform */
	if (dir == 1) {
		for (i=0;i<nn;i++) {
			x[i] /= (float)nn;
			y[i] /= (float)nn;
		}
	}

	return(TRUE);
}
#endif
