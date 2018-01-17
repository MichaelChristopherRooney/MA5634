#include "common.h"

// Uses method described here:
// https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform

static double grv_1;
static double grv_2;
static int iter_count = 0;

static double generate_normal_rv(double mean, double std_dev) {
	if(iter_count == 1){
		iter_count = 0;
		return (grv_2 * std_dev) + mean;
	}
	iter_count++;
	double u1, u2;
	ranlxd(&u1, 1);
	ranlxd(&u2, 1);
	grv_1 = sqrt(-2.0 * log(u1)) * (cos(2.0 * M_PI * u2));
	grv_2 = sqrt(-2.0 * log(u1)) * (sin(2.0 * M_PI * u2));
	return (grv_1 * std_dev) + mean;
}

void do_gaussian_rv(){
	int num_iter = 10000;
	double cosx_est = 0.0;
	double xsqr_est = 0.0;
	int i;
	for(i = 0; i < num_iter; i++){
		double grv = generate_normal_rv(0.0, 0.70785);
		cosx_est += cos_x(grv);
		xsqr_est += x_squared(grv);
	}
	cosx_est = cosx_est / num_iter;
	xsqr_est = xsqr_est / num_iter;
	printf("Cos(x) estimate: %f\n", cosx_est);
	printf("x*x estimate: %f\n", xsqr_est);
}
