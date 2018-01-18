#include "common.h"

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
	int i;
	for(i = 1; i < params->num_iter; i++){ // note: starts at 1
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

void print_met_stats(struct met_params *params){
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

// TODO: how to account for thermalisation here
void estimate_integral(struct met_params *params){
	run_metropolis(params);
	int i;
	for(i = 0; i < params->num_iter; i++){
		// pass x to f(x)
		double est = params->f(params->results[i]);
		params->f_results[i] = est;
		params->estimate += params->f_results[i];
		params->running_estimates[i] = params->estimate / (i + 1);
	}
	params->estimate = params->estimate / params->num_iter;
	if(params->filename != NULL){
		save_results(params->results, params->num_iter, params->filename);
	}
}



