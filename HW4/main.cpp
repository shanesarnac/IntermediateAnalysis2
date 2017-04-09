
#include <iostream> 
#include <math.h>
#include <vector>
#include <cstring>

using namespace std;



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
	/**problem3(); **/
	
	// Problem 4
	problem4();
		
	return 0;
}
