#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <math.h>

#include "ranlxd.h"

void ranlxd(double r[ ],int n);

const double p = 0.5;
const int n = 20;
const int num_runs = 100000;
// Stores the number of occurences for each number at that number's index.
// So if 5 occured 100 times then results[5] will contain 100.
int *results;

int get_binomial_rv(){
	int count = 0;
	double *rn = malloc(sizeof(double) * n);
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec); 
	ranlxd(rn, n);
	int i;
	for(i = 0; i < n; i++){
		if(rn[i] < p){
			count++;
		}
	}
	return count;
}

void gather_results(){
	int i;
	//printf("Peforming %d runs with p=%f and n=%d\n", num_runs, p, n);
	for(i = 0; i < num_runs; i++){
		int bi_rv = get_binomial_rv();
		results[bi_rv] = results[bi_rv] + 1;
	}
}


void display_raw_results(){
	int i;
	printf("Results: (value, number of occurences)\n");
	for(i = 0; i < n; i++){
		printf("%d: %d\n", i, results[i]);
	}
}


void analyse_results(){
	int i;
	for(i = 0; i < n; i++){
		int sum = 0;
		int j;
		for(j = 0; j <= i; j++){
			sum += results[j];
		}
		double chance = (double) sum / (double) num_runs;
		//printf("P(X <= %d) = %f\n", i, chance);
		printf("%d %f\n", i, chance);
	}
}

int main(void){
	results = calloc(n, sizeof(int));
	gather_results();
	//display_raw_results();
	analyse_results();
	return 0;
}

