#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "ranlxd.h"

// Generates an exponentially distributed RV with a given rate using the method described in
// Simulation 3rd edition by Sheldon Ross, page 65
static double generate_exponential_rv(double rate){
	double U;
	ranlxd(&U, 1);
	double ret = (-1/rate);
	ret = ret * log(U);
	return ret;
}

// Generates a poisson distributed RV with a given mean using the method described in
// Simulation 3rd edition by Sheldon Ross, page 51
static int generate_poisson_rv(int mean){
	double U;
	ranlxd(&U, 1);
	int i = 0;
	double p = pow(M_E, mean * -1);
	double F = p;
	while(U >= F){
		p = (mean * p)/(i + 1);
		F = F + p;
		i++;
	}
	return i;
}

// Sets the seed
static void init_ranlux(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec);
}

static const int ACTUAL_AVG_NUM_WITHDRAWALS = 50;
static const int ACTUAL_AVG_AMOUNT_WITHDRAWN = 800;
static const int SIMULATION_ITERATIONS = 1000000;
static const double SIMULATION_TARGET = 50000.00;

int main(){
	init_ranlux();
	double sim_avg_num_withdrawals = 0.0;
	double sim_avg_amount_withdrawn = 0.0;
	int count;
	int i;
	// Each iteration is one month
	for(i = 0; i < SIMULATION_ITERATIONS; i++){
		// Calculate the number of withdrawals and the average amount withdrawn for this month (iteration).
		int iter_withdrawals = generate_poisson_rv(ACTUAL_AVG_NUM_WITHDRAWALS);
		double iter_amount = generate_exponential_rv(1.0 / ACTUAL_AVG_AMOUNT_WITHDRAWN);
		double total = iter_withdrawals * iter_amount;
		// Add to the averages over the entire simulation.
		sim_avg_num_withdrawals = sim_avg_num_withdrawals + (iter_withdrawals);		
		sim_avg_amount_withdrawn = sim_avg_amount_withdrawn + iter_amount;
		if(total > SIMULATION_TARGET){
			count++;
		}
	}
	sim_avg_num_withdrawals = sim_avg_num_withdrawals / SIMULATION_ITERATIONS;
	sim_avg_amount_withdrawn = sim_avg_amount_withdrawn / SIMULATION_ITERATIONS;
	printf("Simulated %d months.\n", SIMULATION_ITERATIONS);
	printf("Average number of withdrawals over the simulation: %f\n", sim_avg_num_withdrawals);
	printf("Average amount per withdrawal over the simulation: %f\n", sim_avg_amount_withdrawn);
	printf("Probability that total withdrawn in a month > 50'000: %f\n", (float) count / (float) i);
	return 0;
}
