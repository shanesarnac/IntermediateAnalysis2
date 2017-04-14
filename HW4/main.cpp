
#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream> 
#include <math.h>
#include <vector>
#include <cstring>

#include "problem1Nodes.h"

using namespace std;

#define PI 3.14159265359

// Problem 1 
void problem1() {
	int L_x = 6;
	int L_y = 4;
	int row = L_y - 1;
	int column = 0;
	
	double dist[L_y][L_x];
	
	// Initialize all points to zero
	for (int i = 0; i < L_y; i++) {
		for (int j = 0; j < L_x; j++) {
			dist[i][j] = 400;
		}
	}
	
	// Iterate until stable
	for (int i = 0; i < 40; i++) {
		row = L_y-1;
		column = 0;
		
		// Node 1
		dist[row][column] = node1(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 2
		dist[row][column] = node2(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 3
		dist[row][column] = node3(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 4
		dist[row][column] = node4(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 5
		dist[row][column] = node5(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 6
		dist[row][column] = node6(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column = 0;
		row--;
		
		// Node 7
		dist[row][column] = node7(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 8
		dist[row][column] = node8(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 9
		dist[row][column] = node9(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 10
		dist[row][column] = node10(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 11
		dist[row][column] = node11(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 12
		dist[row][column] = node12(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column = 0;
		row--;
		
		// Node 13
		dist[row][column] = node13(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 14
		dist[row][column] = node14(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 15
		dist[row][column] = node15(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 16
		dist[row][column] = node16(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 17
		dist[row][column] = node17(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 18
		dist[row][column] = node18(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column = 0;
		row--;

		// Node 19
		dist[row][column] = node19(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 20
		dist[row][column] = node20(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
		
		// Node 21
		dist[row][column] = node21(dist[row+1][column], dist[row-1][column], dist[row][column-1], dist[row][column+1]);
		column++;
	}
	
	for (int i = L_y-1; i >= 0; i--) {
		for (int j = 0; j < L_x-1; j++) {
			cout << dist[i][j] << ",";
		}
		cout << dist[i][L_x-1] << endl;
	}
	

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
#endif
