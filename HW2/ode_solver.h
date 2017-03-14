#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <iostream>
#include <math.h>

using namespace std;

double rk4(double h, double (*f)(double,double), double x_0, double x_max, double y_0);

double rk4_coupled(double h, double (*f1)(double,double,double), double (*f2)(double, double, double), double t_max, double t_0, double x_0, double y_0);

double adams_bashforth_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1);

double adams_bashforth_four_step_predictor(double h, double (*f)(double, double), double x_i, double y_i, double x_im1, double y_im1, 
																				double x_im2, double y_im2, double x_im3, double y_im3);

double adams_moulton_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1);

double adams_moulton_three_step_corrector(double h, double (*f)(double, double), double x_max, double x_0, double y_0, 
																double x_1, double y_1, double x_2, double y_2, double x_3, double y_3);



double improved_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0);

//double implicit_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0);




#endif 
