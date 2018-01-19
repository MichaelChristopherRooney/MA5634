#include "common.h"

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

static double *calculate_bin_means(double *results, struct variance_results *var_results){
	double *bin_means = calloc(var_results->num_bins, sizeof(double));
	int i;
	for(i = 0; i < var_results->num_bins; i++){
		bin_means[i] = calulcate_bin_mean(results, i * var_results->bin_size, var_results->bin_size);
		//printf("%d, %f\n", i, bin_means[i]);
	}
	return bin_means;
}

// First calculates sum of (bin_mean[i] - mean)^2
// Then divides this by Nb * (Nb - 1)
static double calculate_binned_variance(struct variance_results *var_results, double *results){
	double *bin_means = calculate_bin_means(results, var_results);
	double binned_variance = 0.0;
	int i;
	for(i = 0; i < var_results->num_bins; i++){
		double temp = pow(bin_means[i] - var_results->mean, 2);
		binned_variance = binned_variance + temp;
	}
	double div = ((double) var_results->num_bins * (double) (var_results->num_bins - 1));
	binned_variance = binned_variance / div;
	free(bin_means);
	return binned_variance;
}

static double calculate_naive_variance(struct variance_results *var_results, double *results, int num_results){
	double variance = 0.0;
	int i;
	for(i = 0; i < num_results; i++){
		variance = variance + pow(results[i] - var_results->mean, 2);
	}
	variance = variance / num_results;
	return variance;
}

void print_var_stats(struct variance_results *var_results){
	printf("=========================================\n");
	printf("Variance stats for provided results and bin size = %lld\n", var_results->bin_size);
	printf("Number of bins: %lld\n", var_results->num_bins);
	printf("Mean is %f\n", var_results->mean);
	printf("Naive variance is %f\n", var_results->naive_variance);
	printf("Bin variance is %E\n", var_results->bin_variance);
	printf("Integrated autocorrelation time is %E\n", var_results->naive_variance / var_results->bin_variance);
	printf("=========================================\n");
}

struct variance_results *calculate_variances(double *results, double mean, int num_results, int bin_size){
	struct variance_results *var_results = calloc(1, sizeof(struct variance_results));
	var_results->bin_size = bin_size;
	var_results->num_bins = num_results / var_results->bin_size;
	var_results->mean = mean;
	var_results->naive_variance = calculate_naive_variance(var_results, results, num_results);
	var_results->bin_variance = calculate_binned_variance(var_results, results);
	return var_results;

}

