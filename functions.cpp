#include "functions.h"

Fiber::Fiber(){
	x = 0.0;
	y = 0.0;
	z = 0.0;

	d1 = 105;
	d2 = 125;

	n_x = 0.0;
	n_y = 0.0;
	n_z = 1.0;

	N = 100;
}


Fiber::Fiber(double x, double y, double z, double d1, double d2, double n_x, double n_y, double n_z, unsigned int N) : x(x), y(y), z(z), d1(d1), d2(d2), n_x(n_x), n_y(n_y), n_z(n_z), N(N) {
}


void Fiber::generateSources(Generator &generator){


}

Generator::Generator(){
	seed = 0;
	count = 0;
	T = gsl_rng_mt19937; 
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);
}

Generator::Generator(size_t s){
	seed = s;
	count = 0;
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);
}

Generator::~Generator(){
	gsl_rng_free(r);
}

void Generator::set_seed(size_t s){
	gsl_rng_set(r, s);
	seed = s;
}

double Generator::uniform(){
	double ret = gsl_rng_uniform(r);
	count++;
	return ret;
}

double Generator::gaussian(){
	double ret = gsl_ran_ugaussian(r);
	count++;
	return ret;
}
