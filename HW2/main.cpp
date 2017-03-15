
#include <iostream> 
#include <math.h>
#include "ode_solver.h"

using namespace std;

/* Problem 1 Begin */
double realSolution(double t) {
	return 3.0 - 0.998*exp(-1000.0*t) - 2.002*exp(-t);
}
void printSolution(double h, double t0, double tmax) {
	cout << "Real Solution" << endl;
	cout << "t, y" << endl;
	for (double t = t0; t <= tmax; t += h) {
		cout << t << "," << realSolution(t) << endl;
	}
	cout << endl;
}

void implicitEulers(double h, double t0, double y0, double tmax) {
	double y_k, y_kp1, t_kp1, coeff;
	int k = 0;
	
	cout << "Implicit Eulers" << endl;
	cout << "k, t_k, y_k, y_real" << endl;
	cout << k++ << "," << t0 << "," << y0 << "," << realSolution(t0) << endl;
	
	y_k = y0;
	t_kp1 = t0 + h;
	coeff = 1.0 / (1.0 + 1000.0*h);
	
	for (double t_k = t0; t_k < tmax; t_k += h) {
		y_kp1 = coeff * (y_k + 3000.0*h - 2000*h*exp(-t_kp1));
		cout << k++ << "," << t_kp1 << "," << y_kp1 << "," <<  realSolution(t_kp1) << endl;
		y_k = y_kp1;
		t_kp1 += h;
	}
	cout << endl;
}
/* Problem 1 End */

int main() {
	double h = 0.001;
	double t0 = 0;
	double y0 = 0;
	double tmax = 5.0;
	
	implicitEulers(h, t0, y0, tmax);
	
	//printSolution(h, t0, tmax);
	
	return 0;
}
