
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

double L_y = 5;
double L_x = 9;
double h = 0.073;
double k1 = 0.16;
double delta_z = 0.5;
double q = 0.6;
double T_inf = 25;
double T_ref = 20;
double delta_x = 1;
double alpha1 = -delta_z*q*pow(L_x, 2.0)/(T_inf - T_ref);

void finiteDifferences() {
	int x_max = L_x + 1.0;
	int y_max = L_y + 1.0;
	double distribution[y_max][x_max];
	double delta_X = delta_x/L_x;

	
	
	// initialize distribution with zeroes
	for (int i = 0; i < y_max; i++) {
		for (int j = 0; j < x_max; j++) {
			distribution[i][j] = 0;
		}
	}
	
	// Bottom Edge
	for (int i = 1; i < x_max-1; i++) {
		distribution[0][i] = 0.25*(distribution[0][i-1] + distribution[0][i+1] +2*distribution[1][i] - 30*L_x*(delta_X)/(T_inf - T_ref) - alpha1*pow(delta_X, 2.0));
	}
	
	// Interior Points
	for (int i = 1; i < y_max-1; i++) {
		for (int j = 1; j < x_max-1; j++) {
			distribution[i][j] = 0.25*(distribution[i][j-1] + distribution[i][j+1] + distribution[i+1][j] + distribution[i-1][j] - alpha1*pow(delta_X, 2.0)); 
		}
	}
	
	// Top Edge
	for (int i = 1; i < x_max-1; i++) {
		double a = 1.0 + h*L_x*delta_X/(2.0*k1);
		double b = distribution[y_max-1][i-1] + distribution[y_max-1][i+1] + 2*h*L_x*delta_X/k1 + 2.0*distribution[y_max-2][i] - alpha1*pow(delta_X,2.0);
		distribution[y_max-1][i] = 0.25*pow(a, -1.0)*b;
	}
	
	// Print Result (Non-dimensionalized)
	/*cout << endl;
	cout << "Problem 1: Non-Dimensionalized" << endl;
	for (int i = y_max-1; i >= 0; i--) {
		for (int j = 0; j < x_max-1; j++) {
			cout << distribution[i][j] << ",";
		}
		cout << distribution[i][x_max-1] << endl;
	}
	cout << endl;*/
	
	//cout << endl;
	cout << "Problem 1: Re-Dimensionalized" << endl;
	for (int i = y_max-1; i >= 0; i--) {
		for (int j = 0; j < x_max-1; j++) {
			cout << distribution[i][j]*(T_inf - T_ref) + T_ref << ",";
		}
		cout << distribution[i][x_max-1]*(T_inf - T_ref) + T_ref << endl;
	}
	cout << endl;
}


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

void generalMethod(double delta_x, double delta_t, double (*f)(double), double w) {
	int n_nodes_x = 1.0/delta_x + 1;
	int n_nodes_t = 1.0/delta_t;
	double distribution[n_nodes_t][n_nodes_x];
	double redimensionalized[n_nodes_t][n_nodes_x];
	double x_i = 0;
	cout << "A = " << A << endl;
	
	// fill array with initial values
	for (int i = 0; i < n_nodes_x; i++) {
		distribution[0][i] = f(x_i);
		redimensionalized[0][i] = distribution[0][i]*T_s;
		x_i += delta_x;
	}
	
	double r = A*delta_t/pow(delta_x, 2.0);
	cout << "r = " << r << endl;
	for (int i = 0; i < n_nodes_t-1; i++) {
		// set current row to all zeroes
		for (int j = 0; j < n_nodes_x; j++) {
			distribution[i+1][j] = 0;
			redimensionalized[i+1][j] = distribution[i][j]*T_s;
		}
		
		distribution[i+1][n_nodes_x-1] = 0;
		redimensionalized[i+1][n_nodes_x-1] = distribution[i+1][n_nodes_x-1]*T_s;
		
		for (int j = 1; j < n_nodes_x-1; j++) {
			distribution[i+1][j] = w*r*distribution[i+1][j-1] + w*r*distribution[i+1][j+1];
			distribution[i+1][j] += (1-w)*r*(distribution[i][j-1] - 2.0*distribution[i][j] + distribution[i][j+1]);
			distribution[i+1][j] += distribution[i][j];
			distribution[i+1][j] = distribution[i+1][j]/(1.0+2.0*w*r);
			//distribution[i+1][j] = r*distribution[i][j-1] + r*distribution[i][j+1] + (1.0-2.0*r)*distribution[i][j];
			//distribution[i][j] = 0.5*(distribution[i][j-1] + distribution[i][j+1] + ((distribution[i-1][j]*pow(delta_x, 2.0)/(2.0*A*delta_t))));
			redimensionalized[i+1][j] = distribution[i+1][j]*T_s;
		}
	}
	
	double time_elapsed = 0;
	int middle = n_nodes_x/2;
	cout << "middle = " << middle << endl;
	// Non-dimensionalized
	//cout << "Non-dimensionalized" << endl;
	//cout << "t, T" << endl;
	
	//for (int i = 0; i < n_nodes_t; i++) {
		////for (int j = 0; j < n_nodes_x-1; j++) {
			////cout << distribution[i][j] << ",";
		////}
		////cout << distribution[i][n_nodes_x-1] << endl;
		//cout << time_elapsed << ", " << distribution[i][middle] << endl;
		//time_elapsed += delta_t;
	//}
	cout << endl;
	
	
	// Redimensionalized
	cout << "Re-dimensionalized" << endl;
	cout << "t, T" << endl;
	for (int i = 0; i < 10/*n_nodes_t*/; i++) {
		//for (int j = 0; j < n_nodes_x-1; j++) {
			//cout << redimensionalized[i][j] << ",";
		//}
		//cout << redimensionalized[i][n_nodes_x-1] << endl;
		cout << time_elapsed*t_s << ", " << redimensionalized[i][middle] << endl;
		time_elapsed += delta_t;
	}
	cout << endl;
}


/* Problem 2 End */

//////////////////////////////////////

/* Problem 3 Begin */
int grid_width = 7;
int grid_height = 6;

void problem3() {
	cout << "Problem 3!!!" << endl;
	double dist[grid_height][grid_width];
	double B,D,E,F,H;
	double up,down,left,right;
	
	// Fill Distribution with Zeroes
	for (int i = 0; i < grid_height; i++) {
		for (int j = 0; j < grid_width; j++) {
			dist[i][j] = 0;
		}
	}
	
	// Set known values to 1
	dist[5][1] = 1;
	dist[4][4] = 1;
	dist[3][5] = 1;
	dist[0][6] = 1;
	
	for (int i = 1; i < grid_height; i++) {
		for (int j = 0; j < grid_width; j++) {
			up = dist[i+1][j];
			down = dist[i-1][j];
			left = dist[i][j-1];
			right = dist[i][j+1];
			B = 25;
			D = 25;
			E = -100;
			F = 25;
			H = 25;
			
			// Stencil for 17
			if ((i == 1) && (j == 5)) {
				B = 25;
				D = 26.33;
				E = -105.6175;
				F = 29.2875;
				H = 25;
				right = 1;
			}
			// Stencil for 12
			if ((i == 2) && (j == 5)) {
				B = 25;
				D = 31.5925;
				E = -135.8225;
				F = 54.2275;
				H = 25;
				right = 1;
			}
			// Stencil for 3
			if ((i == 4) && (j == 3)) {
				B = 54.2275;
				D = 25;
				E = -135.8225;
				F = 25;
				H = 31.5925;
				up = 1;
			}
			// Stencil for 2 and 2's reflection
			if (((i == 4) && (j == 2)) || ((i == 4) && (j == 0))) {
				B = 29.2875;
				D = 25;
				E = -105.6175;
				F = 25;
				H = 26.33;
				up = 1;
			}
			
			dist[i][j] = -(B*up + D*left + F*right + H*down)/E;
		}
	}
	
	// Print results
	for (int i = grid_height -2; i > 0; i--) {
		cout << "i = " << i << ": ";
		for (int j = 1; j < grid_width-1; j++) {
			cout << dist[i][j] << ",";
		}
		cout << dist[i][grid_width-1] << endl;
	}
}

/* Problem 3 End */

///////////////////////////////////////

/* Problem 4 Begin */


/* Problem 4 End */

int main() {
	
	cout.setf(ios::fixed,ios::floatfield);
    cout.precision(3);
    // Problem 1
    
    /**finiteDifferences();**/
    
    // Problem 2
    /**
	double delta_x = 1.0/8.0;
	//double delta_t = 1.0/48.0; // makes r approx 0.5
	double delta_t = 1.0/24; // makes r approx 1.0
	//double w = 0.0;
	//double w = 2.0/3.0;
	//double w = 0.878;
	double w = 1.0;
	
	generalMethod(delta_x, delta_t, initial_condition, w);
	**/
	
	// Problem 3
	problem3();
	
	
	// Problem 4
		
	return 0;
}
