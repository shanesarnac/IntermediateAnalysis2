
#include <iostream> 
#include <math.h>
#include <vector>
#include <cstring>

using namespace std;

#define PI 3.14159265359

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
	
	// Problem 6
	// part a
	problem6(initialConditions6a);
	
	// part b
	//problem6(initialConditions6b);
		
	return 0;
}
