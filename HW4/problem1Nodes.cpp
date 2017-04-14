#ifndef PROBLEM1_CPP
#define PROBLEM1_CPP

#include "problem1Nodes.h"

// Constants
double T_inf_1 = 1700; //K
double T_inf_2 = 400; //K
double h1 = 1000; // W/(m^2 K)
double h2 = 200; // W/(m^2 K)
double k = 25; // W/ (m K)
double delta_x = 1;
double delta_y = 1;
//double delta_x = 0.001; // m 
//double delta_y = 0.001; // m

double delta_alpha = delta_x/delta_y;
double delta_beta = delta_y/delta_x;


double node1(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*h1*delta_x*T_inf_1 + 0.5*k*(delta_x/delta_y)*T_b + 0.5*k*(delta_y/delta_x)*T_r;
	result = result / (0.5*h1*delta_x + 0.5*k*(delta_x/delta_y) + 0.5*k*(delta_y/delta_x));
	
	return result;
}

double node2(double T_t, double T_b, double T_l, double T_r) {
	double result = h1*delta_x*T_inf_1 + k*(delta_x/delta_y)*T_b + 0.5*k*(delta_y/delta_x)*T_l + 0.5*k*(delta_y/delta_x)*T_r;
	result = result / (h1*delta_x + k*(delta_x/delta_y) + k*(delta_y/delta_x));
	
	return result;
}

double node3(double T_t, double T_b, double T_l, double T_r) {
	return node2(T_t, T_b, T_l, T_r);
}

double node4(double T_t, double T_b, double T_l, double T_r) {
	return node2(T_t, T_b, T_l, T_r);
}

double node5(double T_t, double T_b, double T_l, double T_r) {
	return node2(T_t, T_b, T_l, T_r);
}

double node6(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*h1*delta_x*T_inf_1 + 0.5*k*(delta_x/delta_y)*T_b + 0.5*k*(delta_y/delta_x)*T_l;
	result = result / (0.5*h1*delta_x + 0.5*k*(delta_x/delta_y) + 0.5*k*(delta_y/delta_x));
	
	return result;
}

double node7(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*(delta_x/delta_y)*T_t + 0.5*(delta_x/delta_y)*T_b + (delta_y/delta_x)*T_r;
	result = result / ((delta_x/delta_y) + (delta_y/delta_x));
	
	return result;
}

double node8(double T_t, double T_b, double T_l, double T_r) {
	double result = (delta_y/delta_x)*T_l + (delta_y/delta_x)*T_r + (delta_x/delta_y)*T_t + (delta_x/delta_y)*T_b;
	result = 0.5 * result / ((delta_y/delta_x) + (delta_x/delta_y));
	
	return result;
}

double node9(double T_t, double T_b, double T_l, double T_r) {
	return node8(T_t, T_b, T_l, T_r);
}

double node10(double T_t, double T_b, double T_l, double T_r) {
	return node8(T_t, T_b, T_l, T_r);
}

double node11(double T_t, double T_b, double T_l, double T_r) {
	return node8(T_t, T_b, T_l, T_r);
}

double node12(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*(delta_x/delta_y)*T_t + 0.5*(delta_x/delta_y)*T_b + (delta_y/delta_x)*T_l;
	result = result / ((delta_x/delta_y) + (delta_y/delta_x));
	
	return result;
}

double node13(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*(delta_x/delta_y)*T_t + 0.5*(delta_x/delta_y)*T_b + (delta_y/delta_x)*T_r;
	result = result / ((delta_x/delta_y) + (delta_y/delta_x));
	
	return result;
}

double node14(double T_t, double T_b, double T_l, double T_r) {
	return node8(T_t, T_b, T_l, T_r);
}

double node15(double T_t, double T_b, double T_l, double T_r) {
	double result = k*delta_alpha*T_t + 0.5*k*delta_alpha*T_b + 0.5*h2*delta_x*T_inf_2 + k*delta_beta*T_l + 0.5*k*delta_beta*T_r + 0.5*h2*delta_y*T_inf_2;
	result = result/ (k*delta_alpha + 0.5*k*delta_alpha + 0.5*h2*delta_x + k*delta_beta + 0.5*k*delta_beta + 0.5*h2*delta_y);
	
	return result;
}

double node16(double T_t, double T_b, double T_l, double T_r) {
	double result = k*delta_alpha*T_t + h2*delta_x*T_inf_2 + 0.5*k*delta_beta*T_l + 0.5*k*delta_beta*T_r;
	result = result / (k*delta_alpha + h2*delta_x + 0.5*k*delta_beta + 0.5*k*delta_beta);
	return result;
}

double node17(double T_t, double T_b, double T_l, double T_r) {
	return node16(T_t, T_b, T_l, T_r);
}

double node18(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*k*delta_alpha*T_t + 0.5*h2*delta_x*T_inf_2 + 0.5*k*delta_beta*T_l;
	result = result / (0.5*k*delta_alpha + 0.5*h2*delta_x + 0.5*k*delta_beta);
	return result;
}

double node19(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*delta_alpha*T_t + 0.5*delta_beta*T_r;
	result = result / (0.5*delta_alpha + 0.5*delta_beta);
	return result;
}

double node20(double T_t, double T_b, double T_l, double T_r) {
	double result = delta_x*T_t + 0.5*delta_beta*T_l + 0.5*delta_beta*T_r;
	result = result / (delta_alpha + delta_beta);
	return result;
}

double node21(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*k*delta_alpha*T_t + 0.5*k*delta_beta*T_l + 0.5*h2*delta_y*T_inf_2;
	result = result / (0.5*k*delta_alpha + 0.5*k*delta_beta + 0.5*h2*delta_y);
	return result;
}
#endif
