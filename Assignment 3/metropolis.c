#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "ranlxd.h"

// TODO: pass g(x) and f(x) as parameters

// This is g(x), which is also the PDF
static double g(double x){
	double temp = x * x;
	temp = pow(M_E, -1 * temp);
	temp = temp / sqrt(M_PI);
	return temp;
}

// Simple f(x) for now
static double f(double x){
	return x + 100;
}

// Generates a RV between in range [0, 1] using RANLUX
// Then shifts it to be in range [x - delta, x + delta]
static double generate_x_proposal(double x, double delta){
	double rv;
	ranlxd(&rv, 1);
	double u_min = x - delta;
	double u_max = x + delta;
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
		return 1;
	}
	return 0;
}

// TODO: accept a function for f(x) and g(x)
// TODO: seems to be working "ok" with small delta, but I need to verify it
static double run_metropolis(double x_start, double delta, int num_iter){
	int num_accepted = 0;
	int num_rejected = 0;
	double x = x_start;
	double x_proposal;
	double accept_prob;
	double f_est = 0.0;
	int therm_limit = 500;
	int accepted_after_therm = 0;
	for(int i = 0; i < num_iter; i++){
		x_proposal = generate_x_proposal(x, delta);
		int accept = accept_or_reject(x, x_proposal);
		if(accept == 1){
			//printf("Accepted: %f\n", x_proposal);
			num_accepted++;
			x = x_proposal;
			// is this the right way to sample ?
			if(i > therm_limit){
				accepted_after_therm++;
				double y = f(x);
				printf("x = %f, y = %f\n", x, y);
				f_est = f_est + y;
			}
		} else {
			//printf("Rejected: %f\n", x_proposal);
			num_rejected++;
		}
	}
	f_est = f_est / accepted_after_therm;
	printf("Total accepted: %d\n", num_accepted);
	printf("Total rejected: %d\n", num_rejected);
	printf("f(x) est: %f\n", f_est);
	return 0.0;
}

void init_ranlux(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec); 
}

int main(void){
	init_ranlux();
	run_metropolis(0, 0.2, 10000);
	return 0;
}

