
import math

def function1(x,y):
	return pow(x,2.0) + pow(y, 2.0)
	
def function2(x,y):
	return 5*pow(x,4.0) + 3*pow(x, 2.0) + 2*pow(y, 4.0) + 5*pow(y, 2.0)

def printFunction(f, h, xmin, xmax, ymin, ymax):
	x = xmin
	y = ymin
	printout = "x,y,f(x,y)"
	print(printout)
	while (x <= xmax):
		while (y <= ymax):
			printout = "" + str(x) + "," + str(y) + "," + str(f(x,y))
			print(printout)
			y += h
		x += h
		y = ymin
	

def main():
	xmin = 0.0
	xmax = 2.0
	ymin = 0.0
	ymax = 2.0
	h = 0.1
	#printFunction(function1, h, xmin, xmax, ymin, ymax)
	printFunction(function2, h, xmin, xmax, ymin, ymax)
	
	
main()
