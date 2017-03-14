#ifndef ODE_SOLVER_CPP
#define ODE_SOLVER_CPP

#include "ode_solver.h"

double rk4(double h, double (*f)(double,double), double x_0, double x_max, double y_0) {
	cout << "RK4" << endl;
	double acc = 0.0;
	double y_ip1 = y_0;
	double y_i;
	cout << "y_{i + 1}(" << x_0 << ") = " << y_ip1 << endl;
	for (double x_i = x_0; x_i < x_max; x_i += h) {
		y_i = y_ip1;
		double k1 = h*f(x_i, y_i);
		double k2 = h*f(x_i + h/2.0, y_i + k1/2.0);
		double k3 = h*f(x_i + h/2.0, y_i + k2/2.0);
		double k4 = h*f(x_i + h, y_i + k3);
		y_ip1 = y_i + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
		cout << "k1 = " << k1 << endl;
		cout << "k2 = " << k2 << endl;
		cout << "k3 = " << k3 << endl;
		cout << "k4 = " << k4 << endl;
		cout << "y_{i + 1}(" << x_i + h << ") = " << y_ip1 << endl << endl;
	}
	
	return acc;
}

double rk4_coupled(double h, double (*f1)(double,double,double), double (*f2)(double, double, double), double t_max, double t_0, double x_0, double y_0) {
	cout << "RK4 Coupled" << endl;
	double x_ip1 = x_0;
	double y_ip1 = y_0;
	double x_i, y_i;
	
	cout << "x(" << t_0 << ") = " << x_0 << endl;
	//cout << "y(" << t_0 << ") = " << y_0 << endl;
	
	for(double t_i = t_0; t_i < t_max; t_i += h) {
		x_i = x_ip1;
		y_i = y_ip1;
		
		double k0 = h*f1(t_i, x_i, y_i);
		double l0 = h*f2(t_i, x_i, y_i);
		double k1 = h*f1(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0); 
		double l1 = h*f2(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0); 
		double k2 = h*f1(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1); 
		double l2 = h*f2(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1); 
		double k3 = h*f1(t_i + h, x_i + k2, y_i + l2); 
		double l3 = h*f2(t_i + h, x_i + k2, y_i + l2); 
		
		//cout << "k0 = " << k0 << endl; 
		//cout << "k1 = " << k1 << endl; 
		//cout << "k2 = " << k2 << endl; 
		//cout << "k3 = " << k3 << endl; 
		
		x_ip1 = x_i + (1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3);
		y_ip1 = y_i + (1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3);
		
		cout << "x(" << t_i + h << ") = " << x_ip1 << endl;
		//cout << "y(" << t_i + h << ") = " << y_ip1 << endl;
	} 
	cout << endl;
	
	return 0.0;
}

double adams_bashforth_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1){
	cout << "Adams-Bashforth Two-Step" << endl;
	double acc = 0.0;
	double y_i = y_1;
	double y_im1 = y_0;
	double y_ip1;
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_0 + h << ") = " << y_1 << endl;
	
	for (double x_i = x_0 + h; x_i + h <= x_max; x_i += h) {
		/*cout << "x_i = " << x_i << endl;
		cout << "y_i = " << y_i << endl;
		cout << "y_{i  -1} = " << y_im1 << endl;*/
		
		y_ip1 = y_i + (h/2.0)*(3*f(x_i, y_i) - f(x_i - h, y_im1));
		y_im1 = y_i;
		y_i = y_ip1;

		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 2.0) << endl;
	}
	cout << endl;
	
	return acc;
}

// Note: Only works for y' = y + x (had to solve for y_ip1)
double adams_moulton_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1) {
	cout << "Adams-Moulton Two Step" << endl;
	double acc = 0.0;
	double x_im1 = x_0;
	double y_i = y_1;
	double y_im1 = y_0;
	double y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_0 + h << ") = " << y_1 << endl;
	
	for (double x_i = x_0 + h; x_i + h <= x_max; x_i += h) {
		y_ip1 = (y_i + (h/12.0) * (5.0*(x_i + h) + 8.0*f(x_i, y_i) - f(x_im1, y_im1)))*pow(1.0 - (5.0*h)/12.0, -1.0);
		y_im1 = y_i;
		x_im1 = x_i;
		y_i = y_ip1;
		
		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 2.0) << endl;
	}
	cout << endl;
	
	return acc;
}

double improved_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0) {
	cout << "Improved Euler's Method" << endl;
	double y_i = y_0;
	double y_bar, y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	
	for (double x_i = x_0; x_i < x_max; x_i += h) {
		y_bar = y_i + h*f(x_i, y_i);
		y_ip1 = y_i + (h/2.0)*(f(x_i,y_i) + f(x_i + h, y_bar));
		y_i = y_ip1;
		//cout << "y bar = " << y_bar << endl;
		cout << "y(" << x_i + h<< ") = " << y_ip1 << endl;
	}
	cout << endl;
	
	return 0.0;
}

/*
double implicit_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0) {
	cout << "Implicit Euler's Method" << endl;
	double y_i = y_0;
	double y_bar, y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	
	for (double x_i = x_0; x_i < x_max; i_i += h) {
		y_ip1 = y_i + h*f(x_i + h, y_i + h);
		y_i = y_ip1;
		cout << "y(" << x_i + h<< ") = " << y_ip1 << endl;
	}
	
	return 0.0;
}*/

double adams_bashforth_four_step_predictor(double h, double (*f)(double, double), double x_i, double y_i, double x_im1, double y_im1, double x_im2, double y_im2, double x_im3, double y_im3) {
	return y_i + (h/24.0)*(55.0*f(x_i, y_i) - 59.0*f(x_im1, y_im1) + 37.0*f(x_im2, y_im2) - 9.0*f(x_im3, y_im3));
}

double adams_moulton_three_step_corrector(double h, double (*f)(double, double), double x_max, double x_0, double y_0, double x_1, double y_1, double x_2, double y_2, double x_3, double y_3) {
	cout << "Adams-Moulton Three Step" << endl;
	double acc = 0.0;
	double x_im3 = x_0;
	double x_im2 = x_1;
	double x_im1 = x_2;
	double y_im3 = y_0;
	double y_im2 = y_1;
	double y_im1 = y_2;
	double y_i = y_3;
	double y_ip1, y_bar;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_1 << ") = " << y_1 << endl;
	cout << "y(" << x_2 << ") = " << y_2 << endl;
	
	for (double x_i = x_3; x_i + h <= x_max; x_i += h) {
		y_bar = adams_bashforth_four_step_predictor(h, f, x_i, y_i, x_im1, y_im1, x_im2, y_im2, x_im3, y_im3);
		y_ip1 = y_i + (h/24.0)*(9.0*f(x_i + h, y_bar) + 19.0*f(x_i, y_i) -5.0*f(x_im1, y_im1) + f(x_im2, y_im2));
		y_im3 = y_im2;
		y_im2 = y_im1;
		y_im1 = y_i;
		y_i = y_ip1;
		
		x_im3 = x_im2;
		x_im2 = x_im1;
		x_im1 = x_i;
		
		
		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 5.0) << endl;
	}
	cout << endl;
	
	return acc;
	
}



#endif
