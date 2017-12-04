#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "ranlxd.h"

int generate_bernoulli_rv(double p){
	double rv;
	ranlxd(&rv, 1);
	if(rv < p){
		return 1;
	}
	return 0;
}

double get_p_prime(double p, double e_t_star){
	return (p*e_t_star) / ((p*e_t_star) + (1-p));
}

// Finds the value of e^t*
// Uses a re-arranged equation specified in Simulation pages 169/170
double get_e_t_star(double p, int n, int a){
	return ((p*a) - a) / ((p*a) - (n*p));
}

void init_ranlux(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec); 
}

double get_std_dev(double avg, double *vals, int n){
	double res = 0.0;
	int i;
	for(i = 0; i < n; i++){
		double temp = vals[i] - avg;
		temp = pow(temp, 2);
		res += temp;
	}
	res = res / (n);
	res = sqrt(res);
	return res;
}

#define NUM_RUNS 10000

// Note: a total of NUM_RUNS simulations will be performed using p, n and a
void run_simulation(double p, int n, int a){
	printf("==================================================\n");
	printf("Running %d simulations with p = %f, n = %d, a = %d\n", NUM_RUNS, p, n, a);
	double e_t_star = get_e_t_star(p, n, a);
	printf("Found e^t* = %f\n", e_t_star);
	double p_prime = get_p_prime(p, e_t_star);
	printf("Found p' = %f\n", p_prime);
	printf("Found M(t*) = %f^%d\n", (e_t_star * p) + (1-p), n);
	double avg_sum = 0;
	double *est = malloc(sizeof(double) * NUM_RUNS);
	double avg_est = 0.0;
	int k;
	for(k = 0; k < NUM_RUNS; k++){
		int sum = 0;
		int i;
		for(i = 0; i < n; i++){
			int brv = generate_bernoulli_rv(p_prime);
			sum += brv;
		}
		if(sum >= a){ // calculate the value of the estimator for this iteration and store it
			est[k] = pow((1.0/e_t_star), sum) * pow(((e_t_star * p) + (1-p)), n);
			//printf("Estimator: %f\n", is[k]);
		} else {
			est[k] = 0.0;
		}
		avg_est += est[k];
		avg_sum += sum;
	}
	avg_sum = avg_sum / NUM_RUNS;
	avg_est = avg_est / NUM_RUNS;
	double est_std_dev = get_std_dev(avg_est, est, n);
	double est_var = pow(est_std_dev, 2);
	printf("Average sum: %f\n", avg_sum);
	printf("Average estimator: %.17g\n", avg_est);
	printf("Estimator standard deviation: %.17g\n", est_std_dev);
	printf("Estimator variance: %.17g\n", est_var);
	printf("==================================================\n");
}

int main(void){
	init_ranlux();
	run_simulation(0.4, 20, 16);
	run_simulation(0.2, 20, 16);
	return 0;
}
