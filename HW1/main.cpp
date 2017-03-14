
#include <iostream> 
#include <math.h>
#include "ode_solver.h"

using namespace std;

double y_prime(double x, double y) {
	return -y;
}

double solutionPos(double x) {
	return exp(x);
}

double solutionNeg(double x) {
	return exp(-x);
}

void printSolutionPos(double x0, double xmax, double h) {
	for (double i = x0; i < xmax; i += h) {
		cout << "y(" << i << ") = " << solutionPos(i) << endl;
	}
}

void printSolutionNeg(double x0, double xmax, double h) {
	for (double i = x0; i < xmax; i += h) {
		cout << "y(" << i << ") = " << solutionNeg(i) << endl;
	}
}

void printImplicitEulersHardcoded(double h, double A, double max, double x_0, double y_0) {
	cout << "Implicit Eulers" << endl;
	double y_i, x_i = x_0;
	
	for (double i = 0; i < max; i++) {
		y_i = y_0 * pow((1.0/(1.0 - A*h)), i); 
		cout << "y(" << x_i << ") = " << y_i << endl;
		x_i = x_i + h;
	}
}

void printLeapFrogHardcodedNeg(double h, double A, double max, double x_0, double y_0) {
	cout << "Leap Frog " << endl;
	double y_i, x_i = x_0;
	
	double c1 = 1.0 - (-1.0 + sqrt(2) - exp(-1.0))/(2.0*sqrt(2.0));
	double c2 = (-1.0 + sqrt(2) - exp(-1.0))/(2.0*sqrt(2.0));
	
	
	for (double i = 0; i < max; i++) {
		y_i = c1*pow(A*h + sqrt(pow(A*h, 2.0) + 1.0), i) + c2*pow(A*h - sqrt(pow(A*h, 2.0) + 1.0), i);
		cout << "y(" << x_i << ") = " << y_i << endl;
		x_i = x_i + h;
	}
}

void printLeapFrogHardcodedPos(double h, double A, double max, double x_0, double y_0) {
	cout << "Leap Frog " << endl;
	double y_i, x_i = x_0;
	
	double c2 = (1.0 + sqrt(2.0) - exp(1.0))/(2.0*sqrt(2.0));
	
	double c1 = 1.0 - c2;
	
	for (double i = 0; i < max; i++) {
		y_i = c1*pow(A*h + sqrt(pow(A*h, 2.0) + 1.0), i) + c2*pow(A*h - sqrt(pow(A*h, 2.0) + 1.0), i);
		cout << "y(" << x_i << ") = " << y_i << endl;
		x_i = x_i + h;
	}
}

void printRK2Hardcoded(double h, double A, double max, double x_0, double y_0) {
	cout << "RK2 " << endl;
	
	double y_i, x_i = x_0; 
	
	for (double i = 0; i < max; i++) {
		y_i = y_0*pow(1.0 + A*h + (0.5)*pow(A*h, 2.0), i);
		cout << "y(" << x_i << ") = " << y_i << endl;
		x_i = x_i + h;
	}
}

void printMilnesHardcoded(double h, double A, double max, double x_0, double y_0) {
	cout << "Milnes" << endl;
	
	double y_i, x_i = x_0;
	double c1 = 1.0;
	double c2 = 1.0;
	
	for (double i = 0; i < max; i++) {
		y_i = c1*pow((2.0*A*h/3.0 + sqrt(3.0*pow(A*h/3.0, 2.0) + 1)*(1.0/(1.0 - A*h/3.0))), i);
		y_i += c2*pow((2.0*A*h/3.0 - sqrt(3.0*pow(A*h/3.0, 2.0) + 1)*(1.0/(1.0 - A*h/3.0))), i);
		cout << "y(" << x_i << ") = " << y_i << endl;
		x_i = x_i + h;
	}
}

int main() {
	
	double h = 6.0;
	double A = 1.0;
	//printSolutionNeg(0.0, 15, h); cout << endl;
	printSolutionPos(0.0, 15, h); 
	//printLeapFrogHardcodedNeg(h, A, 20, 0.0, 1.0);
	//printLeapFrogHardcodedPos(h, A, 20, 0.0, 1.0);

	//improved_eulers(h, y_prime, 0.0, 15, 1.0);	
	//printImplicitEulersHardcoded(h, A, 60, 0.0, 1.0);
	
	//printRK2Hardcoded(h, A, 20, 0.0, 1.0);
	printMilnesHardcoded(h, A, 20, 0.0, 1.0);
	
	
	return 0;
}
