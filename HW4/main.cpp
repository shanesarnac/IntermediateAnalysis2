
#include <iostream> 
#include <math.h>
#include <vector>
#include <cstring>

using namespace std;

#define PI 3.14159265359

// Constants
double T_inf_1 = 1700; //K
double T_inf_2 = 400; //K
double h1 = 1000; // W/(m^2 K)
double h2 = 200; // W/(m^2 K)
double k = 25; // W/ (m K)
double delta_x = 0.001; // m 
double delta_y = 0.001; // m

double delta_alpha = delta_x/delta_y;
double delta_beta = delta_y/delta_x;

double node1(double T_t, double T_b, double T_l, double T_r) {
	double result = 0.5*h1*delta_x*T_inf_1 + 0.5*k*(delta_x/delta_y)*T_b + 0.5*k*(delta_y/delta_x)*T_r;
	result = result / (0.5*h1*delta_x + 0.5*k*(delta_x/delta_y) + 0.5*k*(delta_y/delta_x));
	
	return result;
}

double node2Through5(double T_t, double T_b, double T_l, double T_r) {
	double result = h1*delta_x*T_inf_1 + k*(delta_x/delta_y)*T_b + 0.5*k*(delta_y/delta_x)*T_l + 0.5*k*(delta_y/delta_x)*T_r;
	result = result / (h1*delta_x + k*(delta_x/delta_y) + k*(delta_y/delta_x));
	
	return result;
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
	double result = ;
	result = result/ (k*delta_alpha + 0.5*k*delta_alpha + 0.5*h2*delta_x + k*delta_beta + 0.5*k*delta_beta + 0.5*h2*delta_y);
	
	return result;
}


// Problem 1 
void problem1() {

	
	int L_x = 6;
	int L_y = 5;
	int row = L_y - 1;
	int column = 0;
	
	double dist[L_y][L_x];
	
	// Initialize all points to zero
	for (int i = 0; i < L_y; i++) {
		for (int j = 0; j < L_x; j++) {
			dist[i][j] = 0;
		}
	}
	
	// Node 1
	dist[row][column] = node1(0, dist[row-1][column], 0, 0);
	column++;
	
	// Nodes 2-5
	for (int i = column; i < 5; i++) {
		dist[row][column] = node2Through5(0, dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
	} 
	
	// Node 6
	dist[row][column] = node6(0, dist[row-1][column], dist[row][column-1], 0);
	column = 0;
	row--;
	
	// Node 7
	dist[row][column] = node7(dist[row+1][column], dist[row-1][column], 0, dist[row][column+1]);
	column++;
	
	// Node 8
	
	
	
	
	
	
}

// Problem 6
// Part a initial conditions
double initialConditions6a (double x) {
	return sin(2.0*PI*x);
}

// Part b initial conditions
double initialConditions6b (double x) {
	if ((0 <= x) && (x <= 0.5)) {
		return 2.0*x;
	}
	else if ((0.5 < x) && (x <= 1.0)) {
		return 2.0 - 2.0*x;
	}
	else {
		return 0.0;
	}
}

void problem6(double (*f)(double)) {
	double L = 1.0;
	double c = 1.0;
	double delta_x = 0.05;
	double delta_t = 0.025;
	double delta_ratio = delta_t/delta_x;
	int num_x = (int) L/delta_x + 1;
	int num_t = 200;
	double dist[num_t][num_x]; 
	
	int t = 0;
	double T = 0;
	double x = 0;
	// Initialize t = 0
	for (int i = 0; i < num_x; i++) {
		dist[t][i] = f(x);
		x += delta_x;
	}
	
	dist[t][0] = 0;
	dist[t][num_x-1] = 0;
	
	// Using phantom node to determine t = 1 layer
	for (int i = 1; i < num_x-1; i++) {
		dist[t+1][i] = 0.5*(c*pow(delta_ratio, 2.0)*(dist[t][i-1] -2.0*dist[t][i] + dist[t][i+1]) + 2.0*dist[t][i]);
	}
	t++;
	
	for (int i = t; i < num_t-1; i++) {
		dist[i][0] = 0;
		dist[i][num_x-1] = 0;
		for (int j = 1; j < num_x-1; j++) {
			dist[i+1][j] = c*pow(delta_ratio,2.0)*(dist[i][j-1] - 2.0*dist[i][j] + dist[i][j+1]) - dist[i-1][j] + 2.0*dist[i][j];
		}
	}
	
	// print matrix
	cout << "t, ";
	for (int i = 0; i < num_x-1; i++) {
		cout << "x" << i << ",";
	}
	cout << "x" << num_x-1 << endl;
	
	for (int i = 0; i < num_t; i++) {
		cout << T << ",";
		for (int j = 0; j < num_x-1; j++) {
			cout << dist[i][j] << ", ";
		}
		cout << dist[i][num_x-1] << endl;
		T += delta_t;
	}
	cout << endl;
	
}

// Part b initial conditions


int main() {
	
	cout.setf(ios::fixed,ios::floatfield);
    cout.precision(3);

	// Problem 1
	problem1();
	
	
	// Problem 6
	// part a
	/** problem6(initialConditions6a); **/
	
	// part b
	//problem6(initialConditions6b);
		
	return 0;
}
