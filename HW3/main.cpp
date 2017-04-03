
#include <iostream> 
#include <math.h>
#include <vector>
#include <cstring>

using namespace std;

void printArray(int size, double* array) {
	for (int i = 0; i < size-1; i++) {
		cout << array[i] << ",";
	}
	cout << array[size-1] << endl;
}

void printArray2D(int size_x, int size_y, double** array) {
	for (int i = 0; i < size_y; i++) {
		for (int j = 0; j < size_x-1; j++) {
			cout << array[i][j] << ",";
		}
		cout << array[i][size_x-1] << endl;
	}
}

/* Problem 1 Begin */


/* Problem 1 End */


/* Problem 2 Begin*/

double L = 2.0; 
double k = 0.13;
double c = 0.11;
double rho = 7.8;
double alpha = k/(rho*c);
double t_s = 10;
double T_s = 100;
double T_r = 0;
double A = (alpha*t_s)/pow(L, 2.0);

double initial_condition(double x) {
	if ((x >= 0) && (x < 1.0/L)) {
		return L*x;
	}
	if ((x >= 1.0/L) && (x <= 2.0/L)) {
		return 2.0 - L*x;
	}
	else {
		return 0;
	}
}

void dividedDifferences(double delta_x, double delta_t, double (*f)(double)) {
	int n_nodes_x = 1.0/delta_x + 1;
	int n_nodes_t = 1.0/delta_t;
	double distribution[n_nodes_t][n_nodes_x];
	double x_i = 0;
	
	//cout << "nodes: " << n_nodes << endl;
	
	// fill array with initial values
	for (int i = 0; i < n_nodes_x; i++) {
		//cout << "x_i = " << x_i << endl;
		distribution[0][i] = f(x_i);
		x_i += delta_x;
		//cout << initial_vector[i] << endl;
	}
	
	//double * t_im1 = initial_vector;
	//printArray(n_nodes, t_im1);
	
	for (int i = 1; i < n_nodes_t; i++) {
		
		//double t_i[] = new double[n_nodes];
		for (int j = 0; j < n_nodes_x; j++) {
			distribution[i][j] = 0;
		}
		
		for (int j = 1; j < n_nodes_x-1; j++) {
			distribution[i][j] = 0.5*(distribution[i][j-1] + distribution[i][j+1] + ((distribution[i-1][j]*pow(delta_x, 2.0)/(2*A*delta_t))));
		}
		

		//printArray2D(n_nodes_x, n_nodes_t, distribution);
		//printArray(n_nodes, t_im1);
	}
	
	for (int i = 0; i < n_nodes_t; i++) {
		for (int j = 0; j < n_nodes_x-1; j++) {
			cout << distribution[i][j] << ",";
		}
		cout << distribution[i][n_nodes_x-1] << endl;
	}
	
	
	
	
}

/* Problem 2 End */

/* Problem 3 Begin */


/* Problem 3 End */


/* Problem 4 Begin */


/* Problem 4 End */

int main() {
	double delta_x = 1.0/8.0;
	double delta_t = 1.0/10.0;
	
	cout << "A = " << A << endl;
	dividedDifferences(delta_x, delta_t, initial_condition);
		
	return 0;
}
