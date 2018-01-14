#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "ranlxd.h"

// TODO:
// make a param object that is passed around instead of having so many params
// in function calls

// Declaration of a function defined in variance.c
void calculate_variances(double *results, double mean, int num_results, int bin_size);

// This is g(x), which is also the PDF
static double g(double x){
	double temp = x * x;
	temp = pow(M_E, -1 * temp);
	temp = temp / sqrt(M_PI);
	return temp;
}

static double cos_x(double x){
	return cos(x);
}

static double x_squared(double x){
	return x*x;
}

// Generates a RV between in range [0, 1] using RANLUX
// Then shifts it to be in range [x - delta, x + delta]
static double generate_x_proposal(double x, double delta){
	double rv;
	ranlxd(&rv, 1);
	return (x - delta) + (rv * (2 * delta));
}

// Returns 1 for accept, 0 for reject
static int accept_or_reject(double x, double x_proposal){
	double res = g(x_proposal) / g(x);
	if(res > 1){
		return 1;
	}
	double u;
	ranlxd(&u, 1);
	if(res < u){
		return 0;
	}
	return 1;
}

// Generates the markov chain.
// The chain is returned and f(x) is applied to it elsewhere.
static double *run_metropolis(double x_start, double delta, int num_iter, int *num_accepted, int *num_rejected){
	*num_accepted = 0;
	*num_rejected = 0;
	double x = x_start;
	double x_proposal;
	double *results = calloc(num_iter, sizeof(double));
	results[0] = x_start;
	for(int i = 1; i < num_iter; i++){ // note: starts at 1
		x_proposal = generate_x_proposal(x, delta);
		int accept = accept_or_reject(x, x_proposal);
		if(accept == 1){
			(*num_accepted)++;
			x = x_proposal;
		} else {
			(*num_rejected)++;
		}
		results[i] = x;
	}
	return results;
}

static void save_results(double *results, int size, char *filename){
	FILE *fp = fopen(filename, "w+");
	if(fp == NULL){
		printf("Error: unable to write to file %s\n", filename);
		return;
	}
	int i;
	for(i = 0; i < size; i++){
		fprintf(fp, "%d, %f\n", i, results[i]);
	}
	fclose(fp);
}

static void print_stats(char *f_desc, double x_start, double delta, int num_iter, int discard, double estimate, int num_accepted, int num_rejected){
	printf("=========================================\n");
	printf("Statistics for f(x) = %s\n", f_desc);
	printf("Starting x = %f\n", x_start);
	printf("Delta = %f\n", delta);
	printf("Number of iterations = %d\n", num_iter);
	printf("Discarding first %d results to account for thermalisation\n", discard);
	printf("Total accepted: %d\n", num_accepted);
	printf("Total rejected: %d\n", num_rejected);
	printf("Acceptance rate: %f%%\n", ((float) num_accepted / (float) num_iter) * 100.0f);
	printf("Estimate is: %f\n", estimate);
	printf("=========================================\n");
}

// Pass a null filename and the results won't be saved to disk
static double *estimate_integral(double (*f)(double), char *f_desc, double x_start, double delta, int num_iter, int discard, char *filename, int print_results, double *estimate){
	int num_accepted;
	int num_rejected;
	double *f_results = calloc(num_iter, sizeof(double));
	double *results = run_metropolis(x_start, delta, num_iter, &num_accepted, &num_rejected);
	*estimate = 0.0;
	int i;
	for(i = discard; i < num_iter; i++){
		f_results[i] = f(results[i]);
		*estimate += f_results[i];
	}
	*estimate = (*estimate) / num_iter;
	if(filename != NULL){
		save_results(results, num_iter, filename);
	}
	free(results);
	if(print_results == 1){
		print_stats(f_desc, x_start, delta, num_iter, discard, *estimate, num_accepted, num_rejected);
	}
	return f_results;
}

static void init_ranlux(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec); 
}

// This produces the data needed for the graphs in question (b).
static void produce_graph_data(){
	double start_x = 10.0;
	int num_iter = 1000;
	int discard = 100;
	double estimate;
	// cos(x) with delta = 1.5, 3, and 6
	estimate_integral(&cos_x, "cos(x)", start_x, 1.5, num_iter, discard, "data/cosx-1.txt", 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", start_x, 15, num_iter, discard, "data/cosx-2.txt", 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", start_x, 150, num_iter, discard, "data/cosx-3.txt", 1, &estimate);
	// x*x with delta = 1.5, 3, and 6
	estimate_integral(&x_squared, "x*x", start_x, 1.5, num_iter, discard, "data/xsqr-1.txt", 1, &estimate);
	estimate_integral(&x_squared, "x*x", start_x, 15, num_iter, discard, "data/xsqr-2.txt", 1, &estimate);
	estimate_integral(&x_squared, "x*x", start_x, 150, num_iter, discard, "data/xsqr-3.txt", 1, &estimate);
}

static void delta_vs_acceptance(){
	int num_iter = 1000000;
	int discard = 100;
	double estimate;
	estimate_integral(&cos_x, "cos(x)", 0, 0.5, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 1.0, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 2.0, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 5.0, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 10.0, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 25.0, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 50.0, num_iter, discard, NULL, 1, &estimate);	
	estimate_integral(&cos_x, "cos(x)", 0, 100.0, num_iter, discard, NULL, 1, &estimate);
	estimate_integral(&cos_x, "cos(x)", 0, 150.0, num_iter, discard, NULL, 1, &estimate);		
}

static void variance_calulcations(){
	double x_start = 0.0;
	double delta = 2.4;
	int num_iter = 10000000;
	int discard = 100;
	double estimate;
	double *cosx_results = estimate_integral(&cos_x, "cos(x)", x_start, delta, num_iter, discard, NULL, 1, &estimate);
	calculate_variances(cosx_results, estimate, num_iter, 10);
	calculate_variances(cosx_results, estimate, num_iter, 100);
	calculate_variances(cosx_results, estimate, num_iter, 1000);
	calculate_variances(cosx_results, estimate, num_iter, 10000);
	calculate_variances(cosx_results, estimate, num_iter, 100000);
	calculate_variances(cosx_results, estimate, num_iter, 1000000);
}

int main(void){
	init_ranlux();
	variance_calulcations();
	//delta_vs_acceptance();
	//produce_graph_data();
	return 0;
}

