#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// TODO: this isn't correct at the moment

static double calulcate_bin_mean(double *results, int bin_start, int bin_size){
	double bin_mean = 0.0;
	int i;
	for(i = bin_start; i < bin_start + bin_size; i++){
		bin_mean += results[i];
	}
	bin_mean = bin_mean / bin_size;
	//printf("Bin mean is %f\n", bin_mean);
	return bin_mean;
}

static double *calculate_bin_means(double *results, int num_bins, int bin_size){
	double *bin_means = calloc(num_bins, sizeof(double));
	int i;
	for(i = 0; i < num_bins; i++){
		bin_means[i] = calulcate_bin_mean(results, i * bin_size, bin_size);
	}
	return bin_means;
}

static double calculate_binned_variance(double *results, int num_results, int bin_size, double mean){
	int num_bins = num_results / bin_size;
	double *bin_means = calculate_bin_means(results, num_bins, bin_size);
	double binned_variance = 0.0;
	int i;
	for(i = 0; i < num_bins; i++){
		binned_variance = binned_variance + pow(bin_means[i] - mean, 2);
	}
	binned_variance = binned_variance / (num_bins * (num_bins - 1));
	return binned_variance;
}

static double calculate_mean(double *results, int num_results){
	double mean = 0.0;
	int i;
	for(i = 0; i < num_results; i++){
		mean += results[i];
	}
	mean = mean / num_results;
	return mean;
}

static double calculate_naive_variance(double *results, int num_results, double mean){
	double variance = 0.0;
	int i;
	for(i = 0; i < num_results; i++){
		variance = variance + (pow(results[i], 2) - pow(mean, 2));
	}
	variance = variance / num_results;
	return variance;
}

static void print_stats(double mean, double naive_variance, double bin_variance, int bin_size){
	printf("=========================================\n");
	printf("Stats for provided results and bin size = %d\n", bin_size);
	printf("Mean is %f\n", mean);
	printf("Naive variance is %f\n", naive_variance);
	printf("Bin variance is %E\n", bin_variance);
	printf("Integrated autocorrelation time is %f\n", bin_variance / naive_variance);
	printf("=========================================\n");
}

void calculate_variances(double *results, double mean, int num_results, int bin_size){
	//double mean = calculate_mean(results, num_results);
	double naive_variance = calculate_naive_variance(results, num_results, mean);
	double bin_variance = calculate_binned_variance(results, num_results, bin_size, mean);
	print_stats(mean, naive_variance, bin_variance, bin_size);

}

