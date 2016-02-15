#include "functions.h"
#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;

mat bin_photons(vec &x, vec &y, double x0, size_t bins, size_t N){
	vec centers_x = linspace<vec>(-x0,x0,bins);
	vec centers_y = linspace<vec>(-x0,x0,bins);
	double delta = centers_x(1) - centers_x(0);
	mat C = zeros<mat>(bins,bins);
	double temp = 0.0;
	for(size_t i=0; i<N; i++){
		double i_t = (x(i)+x0)/delta;
		double j_t = (y(i)+x0)/delta;
		size_t i_int = (size_t) i_t;
		size_t j_int = (size_t) j_t;
		double i_dec = modf(i_t,&temp);
		double j_dec = modf(j_t,&temp);
		if(i_dec >= 0.5){
			i_int++;
		}
		if(j_dec >= 0.5){
			j_int++;
		}
		C(i_int,j_int) = C(i_int, j_int) + 1.0;
	}
	//C = C/((double) N);
	return C;
}

mat propagate_photons(Generator &gen, size_t seed, size_t N, double r1, double z0, double theta_0){

	gen.set_seed(seed);
	vec x = zeros<vec>(N);
	vec y = zeros<vec>(N);
	double cos_theta_0 = std::cos(theta_0);
	for(size_t i=0; i<N; i++){
		double xs = 2.0*gen.uniform()-1.0;
		double ys = 2.0*gen.uniform()-1.0;
		while((xs*xs+ys*ys) > 1.0){
			xs = 2.0*gen.uniform()-1.0;
			ys = 2.0*gen.uniform()-1.0;
		}
		xs = r1*xs;
		ys = r1*ys;
		double psi = 2.0*M_PI*gen.uniform();
//		double theta = gen.uniform()*theta_0;
//		double theta = std::acos(2.0*gen.uniform() - 1.0);
		double theta = std::acos((1.0-cos_theta_0)*gen.uniform() + cos_theta_0);
		double cos_phi = std::cos(psi);
		double sin_phi = std::sin(psi);
		double tan_theta = std::tan(theta);
		x(i) = xs + z0*cos_phi*tan_theta;
		y(i) = ys + z0*sin_phi*tan_theta;
	}
	mat positions = zeros<mat>(N,2);
	positions.col(0) = x;
	positions.col(1) = y;
	return positions;
}	



int main(){

	Generator G;
	double r1 = 0.0525;
	double z0 = 0.2;
	double theta_0 = 0.22;
	double x0 = r1 + z0*std::tan(theta_0);
	size_t N = 10000000;
	size_t bins = 101;
	size_t seed = 3;
	mat pos = propagate_photons(G, seed, N, r1, z0, theta_0);

	vec x = pos.col(0);
	vec y = pos.col(1);
	mat distr = bin_photons(x, y, x0, bins, N);
	pos = pos/r1;
	pos.save("positions.dat", raw_ascii);
	vec counts1 = distr.col((bins-1)/2);
	counts1 = counts1/(accu(counts1));
	vec x_centers = linspace<vec>(-x0, x0, bins);
	uvec counts_2 = arma::hist(y, x_centers);
	vec counts2 = conv_to<vec>::from(counts_2);
	counts2 = counts2/((double) N);
	counts2.save("y2.dat", raw_ascii);
	counts1.save("y.dat", raw_ascii);
	distr.save("distribution.dat", raw_ascii);




	/*
	double N_double = (double) N;
	size_t bins = 100;
	Generator gen(10);
	vec R = zeros<vec>(N);
	vec X = zeros<vec>(N);
	vec Y = zeros<vec>(N);
	vec Xs = zeros<vec>(N);
	vec Ys = zeros<vec>(N);
	for(size_t i=0; i<N; i++){
		double xs = 2.0*gen.uniform()-1.0;
		double ys = 2.0*gen.uniform()-1.0;
		while((xs*xs+ys*ys) > 1.0){
			xs = 2.0*gen.uniform()-1.0;
			ys = 2.0*gen.uniform()-1.0;
		}
		double psi = 2.0*M_PI*gen.uniform();
		Xs(i) = std::cos(psi);
		Ys(i) = std::sin(psi);
		xs = r1*xs;
		ys = r1*ys;
		double theta = gen.uniform()*theta_0;
		psi = 2.0*M_PI*gen.uniform();
		double cos_phi = std::cos(psi);
		double sin_phi = std::sin(psi);

		double tan_theta = std::tan(theta);
		X(i) = xs + z0*cos_phi*tan_theta;
		Y(i) = ys + z0*sin_phi*tan_theta;
		R(i) = std::sqrt(std::pow(xs+z0*cos_phi*tan_theta,2.0) + std::pow(ys+z0*sin_phi*tan_theta,2.0));
	}
	//std::cout << R << std::endl;
	std::cout << R.max() << std::endl;
	std::cout << "numbers generated " << gen.count << std::endl;

	vec edge = linspace<vec>(R.min(), R.max(), bins);
	double dr = edge(1)-edge(0);
	vec area = M_PI*(square(edge.subvec(1,bins-1))-square(edge.subvec(0,bins-2)));
	vec radius = edge.subvec(0,bins-2) + 0.5*dr;
	uvec h0 = histc(R,edge);
	vec h0_new = conv_to<vec>::from(h0);
	vec counts = h0_new.subvec(0,bins-2);
	counts = counts/area;
	counts = counts/accu(counts);
	//counts.save("histogram.dat", raw_ascii);
	mat stats = zeros<mat>(bins-1,2);
	stats.col(0) = radius;
	stats.col(1) = counts;
	stats.save("histogram.dat", raw_ascii);

	mat XY = zeros<mat>(N,2);
	XY.col(0) = X;
	XY.col(1) = Y;
	XY.save("density.dat", raw_ascii);

	mat XYs = zeros<mat>(N,2);
	XY.col(0) = Xs;
	XY.col(1) = Ys;
	XY.save("densityS.dat", raw_ascii);

	double grid_max = r1 + z0*std::tan(theta_0);
	cout << grid_max << endl;

	//Binning
	size_t g_size = 301;
	vec x_g = linspace<vec>(-x0,x0,g_size);
	vec y_g = linspace<vec>(-x0,x0,g_size);
	double delta = x_g(1)-x_g(0);
	umat C = zeros<umat>(g_size,g_size);
	for(size_t i=0; i<N; i++){
		double i_t = (X(i)+x0)/delta;
		double j_t = (Y(i)+x0)/delta;
		size_t i_int = (size_t) i_t;
		size_t j_int = (size_t) j_t;
		double dec = 0.0;
		double i_dec = modf(i_t,&dec);
		double j_dec = modf(j_t,&dec);
		if(i_dec >= 0.5){
			i_int++;
		}
		if(j_dec >= 0.5){
			j_int++;
		}
		C(i_int,j_int)++;
	}
	uvec counts_x = hist(X,x_g);
	counts_x.save("countsX",raw_ascii);
	cout << arma::accu(C) << endl;
	C.save("3d.dat",raw_ascii);

	uvec Cx = C.col((g_size-1)/2);

	vec c = conv_to<vec>::from(Cx);
	c = c/N_double;
	mat distr(g_size,2);
	cout << size(y_g) << " " << size(c) << endl;
	distr.col(0) = y_g;
	distr.col(1) = c;
	distr.save("counts.dat",raw_ascii);
	*/











	return 0;
}
