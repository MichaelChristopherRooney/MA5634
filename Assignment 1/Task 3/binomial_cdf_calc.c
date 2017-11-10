// A simple program for calculating the expected CDF of a binomial RV
// Uses the equation for the expected CDF found here:
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda366i.htm 
#include <math.h>

const int n = 20;
const double p = 0.5;

// Fails for large values of x, but it is ok for this simple use case
long long calculate_factorial(long long x){
	long long result = 1;
	int k;
	for(k = 1; k <= x; k++){
		result = result * k;
	}
	return result;
}

int main(void){
	int i;
	double sum = 0.0;
	for(i = 0; i < n; i++){
		double temp = 0.0;
		// First calculate ( n! / (i!*(n-i)!)
		long long n_fact = calculate_factorial(n);
		long long i_fact = calculate_factorial(i);
		long long n_i_fact = calculate_factorial(n - i);
		long long fact_part = n_fact / (i_fact * n_i_fact);
		// Now calculate (p^i) * ((1-p)^(n - i))
		double p_pow_i = pow(p, i);
		double one_p_pow_n_i = pow(1-p, n-i);
		// Start multipling all the parts together
		temp = p_pow_i * one_p_pow_n_i;
		temp = temp * fact_part;
		sum += temp;
		printf("%d %f\n", i, sum);
		
	}
	return 0;
}

