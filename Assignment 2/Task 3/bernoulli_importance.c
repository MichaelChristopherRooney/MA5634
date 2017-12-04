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
	int k;
	for(k = 0; k < NUM_RUNS; k++){
		int sum = 0;
		int i;
		for(i = 0; i < n; i++){
			int brv = generate_bernoulli_rv(p_prime);
			sum += brv;
		}
		avg_sum += sum;
	}
	avg_sum = avg_sum / NUM_RUNS;
	printf("Average sum: %f\n", avg_sum);
	printf("==================================================\n");
}

int main(void){
	init_ranlux();
	run_simulation(0.4, 20, 16);
	run_simulation(0.2, 20, 16);
	return 0;
}
