#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "ranlxd.h"

// TODO:
// make a param object that is passed around instead of having so many params
// in function calls

struct met_params {
	double x_start;
	double delta;
	int num_iter;
	int discard;
	char *filename;
	int num_accepted;
	int num_rejected;
	double *results;
	double *f_results;
	double estimate;
	char *f_desc;
	double (*f)(double);
};

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
static void run_metropolis(struct met_params *params){
	double x = params->x_start;
	double x_proposal;
	params->results[0] = params->x_start;
	for(int i = 1; i < params->num_iter; i++){ // note: starts at 1
		x_proposal = generate_x_proposal(x, params->delta);
		int accept = accept_or_reject(x, x_proposal);
		if(accept == 1){
			params->num_accepted++;
			x = x_proposal;
		} else {
			params->num_rejected++;
		}
		params->results[i] = x;
	}
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

static void print_stats(struct met_params *params){
	printf("=========================================\n");
	printf("Statistics for f(x) = %s\n", params->f_desc);
	printf("Starting x = %f\n", params->x_start);
	printf("Delta = %f\n", params->delta);
	printf("Number of iterations = %d\n", params->num_iter);
	printf("Discarding first %d results to account for thermalisation\n", params->discard);
	printf("Total accepted: %d\n", params->num_accepted);
	printf("Total rejected: %d\n", params->num_rejected);
	printf("Acceptance rate: %f%%\n", ((float) params->num_accepted / (float) params->num_iter) * 100.0f);
	printf("Estimate is: %f\n", params->estimate);
	printf("=========================================\n");
}

static struct met_params *init_params(double (*f)(double), char *f_desc, double x_start, double delta, int num_iter, int discard, char *filename){
	struct met_params *params = calloc(1, sizeof(struct met_params));
	params->x_start = x_start;
	params->delta = delta;
	params->num_iter = num_iter;
	params->discard = discard;
	params->filename = filename;
	params->num_accepted = 0;
	params->num_rejected = 0;
	params->f_results = calloc(num_iter, sizeof(double));
	params->results = calloc(num_iter, sizeof(double));
	params->estimate = 0.0;
	params->f = f;
	params->f_desc = f_desc;
	return params;
}

static void free_params(struct met_params *params){
	free(params->results);
	free(params->f_results);
	free(params);
}

// Pass a null filename and the results won't be saved to disk
static void estimate_integral(struct met_params *params){
	run_metropolis(params);
	int i;
	for(i = params->discard; i < params->num_iter; i++){
		params->f_results[i] = params->f(params->results[i]);
		params->estimate += params->f_results[i];
	}
	params->estimate = params->estimate / params->num_iter;
	if(params->filename != NULL){
		save_results(params->results, params->num_iter, params->filename);
	}
	print_stats(params);
}

static void init_ranlux(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec); 
}

// This produces the data needed for the graphs in question (b).
// TODO: use params here
static void produce_graph_data(){
	/*
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
	*/
}

// TODO: use params here
static void delta_vs_acceptance(){
	/*
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
	*/
}

static void variance_calulcations(){
	struct met_params *params = init_params(&cos_x, "cos(x)", 0.0, 2.4, 10000000, 100, NULL);
	estimate_integral(params);
	//calculate_variances(cosx_results, estimate, num_iter, 10);
	//calculate_variances(cosx_results, estimate, num_iter, 100);
	//calculate_variances(cosx_results, estimate, num_iter, 1000);
	//calculate_variances(cosx_results, estimate, num_iter, 10000);
	//calculate_variances(cosx_results, estimate, num_iter, 100000);
	//calculate_variances(cosx_results, estimate, num_iter, 1000000);
	free_params(params);
}

int main(void){
	init_ranlux();
	variance_calulcations();
	//delta_vs_acceptance();
	//produce_graph_data();
	return 0;
}

