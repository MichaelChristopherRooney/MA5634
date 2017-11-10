#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "ranlxd.h"

void ranlxd(double r[ ],int n);

// g(x) = x
double func_1(double x){
	return x;
}

// g(x) = x * x
double func_2(double x){
	return x * x;
}

// g(x) = sqrt(x);
double func_3(double x){
	return sqrt(x);
}

// g(x) = (2 / (sqrt(1 - (x*x)))
double func_4(double x){
	double temp;
	temp = x * x;
	temp = 1 - temp;
	temp = sqrt(temp);
	temp = 2 / temp;
	return temp;		
}

// Array of the functions, their descriptions and their actual values when integrated from 0 to 1.
// Used so we can iterate over them, print them and call them easily.
#define NUM_FUNCS 4
const char *FUNC_NAMES[NUM_FUNCS] = { "g(x) = x", "g(x) = x*x", "g(x) = sqrt(x)", "g(x) = (2/(sqrt(1-(x*x)))" };
double (*FUNCS[NUM_FUNCS]) (double x) = {func_1, func_2, func_3, func_4 };
double FUNC_ACTUAL_RESULT[NUM_FUNCS] = { 0.5, 0.33333333, 0.66666666, M_PI };


// Limits of the integrals
const double a = 0;
const double b = 1;

// The number of random numbers to be used in each test
#define NUM_SIZES 7
long long SIZES[NUM_SIZES] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };


int main(void){
	int i;
	for(i  = 0; i < NUM_SIZES; i++){
		int n = SIZES[i];
		double *rn = malloc(sizeof(double) * n);
		struct timeval tv;
		gettimeofday(&tv, NULL);
		rlxd_init(1, tv.tv_usec); 
		ranlxd(rn, n);
		int j;
		printf("================================================\n");		
		printf("n=%d\n", n);
		printf("FUNCTION\t\t\tRESULT\n");
		printf("================================================\n");
		for(j = 0; j < NUM_FUNCS; j++){
			double result = 0.0;
			int k;
			for(k = 0; k < n; k++){
				// Var[g(x)] = E[g(x)^2] - μ^2
				result = result + ((pow(FUNCS[j](rn[k]) * (b - a), 2)) - pow(FUNC_ACTUAL_RESULT[j], 2)) / n;
			}
			printf("%s\t\t\t%f\n", FUNC_NAMES[j], result);
		}
		printf("\n");
		free(rn);
	}
	return 0;
}

