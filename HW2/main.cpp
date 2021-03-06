
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

void gearsThreePoint(double h, double t0, double y0, double t1, double y1, double t2, double y2, double tmax) {
	double y_km2, y_km1, y_k, y_kp1, t_kp1, coeff;
	
	cout << "Gears Three Point" << endl;
	cout << "t_k, y_k, y_real" << endl;
	cout << t2 << "," << y2 << "," << realSolution(t2) << endl;
	// How do we determine y1 and y2?
	
	y_km2 = y0;
	y_km1 = y1;
	y_k = y2;
	t_kp1 = t2 + h;
	coeff = 1.0 / (1.0 + 6000.0*h/11.0);
	
	for (double t_k = t2; t_k < tmax; t_k += h) {
		y_kp1 = (1.0/11.0)* coeff * (2.0*y_km2 - 9.0*y_km1 + 18.0*y_k + 6.0*h*(3000.0 - 2000*exp(-t_kp1)));
		cout << t_kp1 << "," << y_kp1 << "," <<  realSolution(t_kp1) << endl;
		y_km2 = y_km1;
		y_km1 = y_k;
		y_k = y_kp1;
		t_kp1 += h;
	}
	cout << endl;
}

/* Problem 1 End */

/* Problem 2 Begin */
double PI = 3.1415926535897;
double L = 0.1;
double d = 0.005;
double k = 400;
double Tb = 100;
double Tinf = 25;
double conv_h = 100;
double Area = PI*pow(0.5*d, 2);
double P = PI*d;
double m = conv_h*P/(k*Area);

double realSolutionProblem2(double X) {
	//double X = x/L;
	return (cosh(m*L*(1.0-X)) + (conv_h/(m*k))*sinh(m*L*(1.0-X))) / (cosh(m*L) + (conv_h/(m*k))*sinh(m*L));
}

void printSolutionProblem2() {
	cout << "x, y_real" << endl;
	for (double x = 0.0; x <= 0.1; x += 0.01) {
		cout << x << "," << realSolutionProblem2(x) << endl;
	}
	cout << endl;
}

void dividedDifferences() {
	double theta[] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double deltaX = 0.1;
	double delta = 2.0 + pow(m*L*deltaX, 2.0);
	double x = 0.0;
	
	for (int i = 0; i < 100; i++) {
		theta[0] = 1;
		for (int j = 1; j < 10; j++) {
			theta[j] = (theta[j-1] + theta[j+1])/delta;
		}
		theta[10] = 2*theta[9]/(delta + 2*(conv_h*L*deltaX/k));
		
		//cout << "iteration " << i << endl;

	}
	
	cout << "x, theta_est, theta_real" << endl;
	x = 0.0;
	
	for (int k = 0; k < 11; k++) {
		cout << x << "," << theta[k] << "," << realSolutionProblem2(x) << endl;
		x += deltaX;
	}
	cout << endl;
}

/* Problem 2 End */

int main() {
	/*double h = 0.01;
	double t0 = 0;
	double y0 = 0;
	double tmax = 5.0;
	
	double tm1 = t0 - 0.0001;
	double tm2 = tm1 - 0.0001; 
	double ym1 = 1.0-exp(-1000.0*tm1);
	double ym2 = 1.0-exp(-1000.0*tm2);
	
	//implicitEulers(h, t0, y0, tmax);
	gearsThreePoint(h, tm2, ym2, tm1, ym1, t0, y0, tmax);
	
	//printSolution(h, t0, tmax);*/
	
	//printSolutionProblem2();
	
	dividedDifferences();
		
	return 0;
}
