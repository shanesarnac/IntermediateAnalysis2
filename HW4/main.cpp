
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
double intitialCondition6b (double x) {
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

void problem6(double ((*f)(double)) {
	
	
}

// Part b initial conditions


int main() {
	
	cout.setf(ios::fixed,ios::floatfield);
    cout.precision(3);

	// Problem 1
	
	// Problem 6
	
		
	return 0;
}
