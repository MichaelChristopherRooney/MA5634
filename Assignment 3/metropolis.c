#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "ranlxd.h"

// Declaration of a function defined in variance.c
void calculate_variances(double *results, int num_results, int bin_size);


// TODO: pass g(x) and f(x) as parameters

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

// TODO: accept a function for f(x) and g(x)
// TODO: seems to be working "ok" with small delta, but I need to verify it
static double *run_metropolis(double x_start, double delta, int num_iter){
	int num_accepted = 0;
	int num_rejected = 0;
	double x = x_start;
	double x_proposal;
	double *results = calloc(num_iter, sizeof(double));
	results[0] = x_start;
	for(int i = 1; i < num_iter; i++){ // note: starts at 1
		x_proposal = generate_x_proposal(x, delta);
		int accept = accept_or_reject(x, x_proposal);
		if(accept == 1){
			num_accepted++;
			x = x_proposal;
		} else {
			num_rejected++;
		}
		results[i] = x;
	}
	printf("Total accepted: %d\n", num_accepted);
	printf("Total rejected: %d\n", num_rejected);
	printf("Acceptance rate: %f%%\n", ((float) num_accepted / (float) num_iter) * 100.0f);
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
	printf("Markov chain saved to %s\n", filename);
}


// Pass a null filename and the results won't be saved to disk
// Pass bin size = 0 to avoid calculating variance
static double estimate_integral(double (*f)(double), char *f_desc, double x_start, double delta, int num_iter, int discard, char *filename, int bin_size){
	printf("=========================================\n");
	printf("Statistics for f(x) = %s\n", f_desc);
	printf("Starting x = %f\n", x_start);
	printf("Delta = %f\n", delta);
	printf("Number of iterations = %d\n", num_iter);
	printf("Discarding first %d results to account for thermalisation\n", discard);
	double *results = run_metropolis(x_start, delta, num_iter);
	double estimate = 0.0;
	int i;
	for(i = discard; i < num_iter; i++){
		estimate += f(results[i]);
	}
	estimate = estimate / num_iter;
	if(filename != NULL){
		save_results(results, num_iter, filename);
	}
	if(bin_size != 0){
		calculate_variances(results, num_iter, bin_size);
	} else {
		printf("Not calculating variance\n");
	}
	printf("Estimate is: %f\n", estimate);
	printf("=========================================\n");
	free(results);
	return estimate;
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
	// cos(x) with delta = 1.5, 3, and 6
	estimate_integral(&cos_x, "cos(x)", start_x, 1.5, num_iter, discard, "data/cosx-1.txt", 0);
	estimate_integral(&cos_x, "cos(x)", start_x, 15, num_iter, discard, "data/cosx-2.txt", 0);
	estimate_integral(&cos_x, "cos(x)", start_x, 150, num_iter, discard, "data/cosx-3.txt", 0);
	// x*x with delta = 1.5, 3, and 6
	estimate_integral(&x_squared, "x*x", start_x, 1.5, num_iter, discard, "data/xsqr-1.txt", 0);
	estimate_integral(&x_squared, "x*x", start_x, 15, num_iter, discard, "data/xsqr-2.txt", 0);
	estimate_integral(&x_squared, "x*x", start_x, 150, num_iter, discard, "data/xsqr-3.txt", 0);
}

static void delta_vs_acceptance(){
	int num_iter = 1000000;
	int discard = 100;
	estimate_integral(&cos_x, "cos(x)", 0, 0.5, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 1.0, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 2.0, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 5.0, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 10.0, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 25.0, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 50.0, num_iter, discard, NULL, 0);	
	estimate_integral(&cos_x, "cos(x)", 0, 100.0, num_iter, discard, NULL, 0);
	estimate_integral(&cos_x, "cos(x)", 0, 150.0, num_iter, discard, NULL, 0);		
}

int main(void){
	init_ranlux();
	float delta = 2.4;
	estimate_integral(&cos_x, "cos(x)", 0, delta, 1000000, 100, NULL, 100);
	//estimate_integral(&x_squared, "x*x", 0, delta, 10000000, 100, NULL, 0);
	//delta_vs_acceptance();
	//produce_graph_data();
	return 0;
}

