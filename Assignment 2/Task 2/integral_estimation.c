#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "ranlxd.h"

void ranlxd(double r[ ],int n);

// f(s) = (e^-s) / (1 + (s/9))
double func(double s){
	double temp1 = pow(M_E, -s);
	double temp2 = 1 + (s / 9);
	return (temp1 / temp2);
}

// Limits of the integrals
const double a = 0;
const double b = 3;

const double ACTUAL_RESULT = 0.873109;

// The number of random numbers to be used in each test
#define NUM_SIZES 4
long long SIZES[NUM_SIZES] = { 100, 1000, 10000, 100000 };


int main(void){
	int i;
	for(i  = 0; i < NUM_SIZES; i++){
		int n = SIZES[i];
		double *rn = malloc(sizeof(double) * n);
		struct timeval tv;
		gettimeofday(&tv, NULL);
		rlxd_init(1, tv.tv_usec); 
		ranlxd(rn, n);
		double result = 0.0;
		int k;
		for(k = 0; k < n; k++){
			result = result + ((func(rn[k] * 3) * (b - a)) / n);
		}
		free(rn);
		printf("n=%d, result=%f, absolute error=%f\n", n, result, fabs(result - ACTUAL_RESULT));
	}
	return 0;
}

