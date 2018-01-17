#include "common.h"

// TODO: make functions save data to file
// TODO: rename functions to be more descriptive
// TODO: add comments
// TODO: make sure thermalisation is being accounted for

// This is g(x), which is also the PDF
double g(double x){
	double temp = x * x;
	temp = pow(M_E, -1 * temp);
	temp = temp / sqrt(M_PI);
	return temp;
}

double cos_x(double x){
	return cos(x);
}

double x_squared(double x){
	return x*x;
}

static void init_ranlux(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	rlxd_init(1, tv.tv_usec); 
}

static struct met_params *init_params(double (*f)(double), char *f_desc, double x_start, double delta, int num_iter, int discard, char *filename){
	struct met_params *params = calloc(1, sizeof(struct met_params));
	params->x_start = x_start;
	params->delta = delta;
	params->num_iter = num_iter;
	params->discard = discard;
	params->filename = filename;
	params->num_accepted = 0;
	params->num_rejected = 0;
	params->f_results = calloc(num_iter, sizeof(double));
	params->results = calloc(num_iter, sizeof(double));
	params->estimate = 0.0;
	params->f = f;
	params->f_desc = f_desc;
	return params;
}

static void free_params(struct met_params *params){
	free(params->results);
	free(params->f_results);
	free(params);
}

// This produces the data needed for the graphs in question (b).
static void create_history_data(){
	double (*fs[2])(double) = { &cos_x, &x_squared };
	char *f_str[2] = {"cos(x)", "x*x"};
	char *file_names[2] = {"cosx-", "xsqr-"};
	double deltas[3] = { 1.5, 15.0, 150.0 };
	char buf[128];
	int i, n, j;
	for(i = 0; i < 2; i++){
		for(n = 0; n < 3; n++){
			struct met_params *params = init_params(fs[i], f_str[i], 0.0, deltas[n], 1000000, 100, NULL);
			sprintf(buf, "data/%s%d.txt", file_names[i], n + 1);
			FILE *fp = fopen(buf, "w");
			estimate_integral(params);
			for(j = 0; j < params->num_iter; j++){
				fprintf(fp, "%d, %f\n", j, params->results[j]);
			}
			fclose(fp);
			free_params(params);
		}
	}
}

static void create_delta_vs_acceptance_data(){
	FILE *fp = fopen("data/delta_vs_acceptance_rate.txt", "w");
	double deltas[] = { 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0, 150.0 };
	int i;
	for(i = 0; i < sizeof(deltas) / sizeof(deltas[0]); i++){
		struct met_params *params = init_params(&cos_x, "cos(x)", 0.0, deltas[i], 1000000, 100, NULL);
		estimate_integral(params);
		fprintf(fp, "%f, %f\n", deltas[i], ((float) params->num_accepted / (float) params->num_iter) * 100.0f);
		free_params(params);
	}
	fclose(fp);
}

static void variance_calulcations(){
	struct met_params *params = init_params(&cos_x, "cos(x)", 0.0, 2.4, 10000000, 100, NULL);
	estimate_integral(params);
	int bin_sizes[] = { 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000 };
	int i;
	for(i = 0; i < sizeof(bin_sizes) / sizeof(bin_sizes[0]); i++){
		calculate_variances(params->f_results, params->estimate, params->num_iter, bin_sizes[i]);
	}
	free_params(params);
}

//TODO: command line arguments that let you pick what code to run
int main(void){
	init_ranlux();
	//variance_calulcations();
	//create_delta_vs_acceptance_data();
	create_history_data();
	do_gaussian_rv();
	return 0;
}
