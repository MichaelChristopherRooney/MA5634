#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "ranlxd.h"

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
	double *running_estimates;
	char *f_desc;
	double (*f)(double);
};

double g(double x);
double cos_x(double x);
double x_squared(double x);

struct variance_results {
	long long bin_size;
	long long num_bins;
	double bin_variance;
	double naive_variance;
	double mean;
};

void calculate_variances(double *results, double mean, int num_results, int bin_size);

void estimate_integral(struct met_params *params);
void print_met_stats(struct met_params *params);

void do_gaussian_rv();
