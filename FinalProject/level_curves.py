import sys
import math
import matplotlib.pyplot as plt

class LevelCurvesSolver:
	def __init__(self, h, x, y, function_data):
		self.levelCurveX = [x]
		self.levelCurveY = [y]
		self.minX = min(function_data.x)
		self.maxX = max(function_data.x) + function_data.deltaX
		self.minY = min(function_data.y) 
		self.maxY = max(function_data.y) + function_data.deltaY
		self.findLevelCurve(h, x, y, function_data)
		
		max_index = 1500
		
		for i in range(max_index):
			x2 = self.levelCurveX[-1]
			y2 = self.levelCurveY[-1]
			self.findLevelCurve(h, x2, y2, function_data)
		self.findLevelCurve(h, x, y, function_data)
		for i in range(max_index):
			x1 = self.levelCurveX[-2]
			y1 = self.levelCurveY[-2]
			self.findLevelCurve(h, x1, y1, function_data)
			
		self.plotLevelCurve()
		
	def findLevelCurve(self, h, x, y, function_data):
		factory = StrategyFactory(function_data)
		strategy = factory.getStrategy(x, y)
		
		print(strategy.getStrategyName())
		
		(x1, y1, x2, y2) = strategy.findValue(h,x,y)
		
		if function_data.isPointInData(x1,y1):
			self.levelCurveX.append(x1)
			self.levelCurveY.append(y1)
		if function_data.isPointInData(x2,y2):
			self.levelCurveX.append(x2)
			self.levelCurveY.append(y2)
	
	def plotLevelCurve(self):
		plt.plot(self.levelCurveX, self.levelCurveY, "ro")
		plt.xlabel("x")
		plt.ylabel("y")
		plt.title("Level Curve")
		plt.axis([self.minX, self.maxX, self.minY, self.maxY])
		plt.show()
		

# A factory class for generating the proper strategy 
class StrategyFactory:
	def __init__(self, function_data):
		self.data = function_data
		
	def determineEdgeStrategy(self, x, y):
		if float(y) == self.data.y[0]:
				if float(x) in self.data.x:
					return BottomEdgeVertexStrategy(self.data)
				else:
					return BottomEdgeMiddlePointStrategy(self.data)
		if float(y) == self.data.y[-1]:
			if float(x) in self.data.x:
				return TopEdgeVertexStrategy(self.data)
			else:
				return TopEdgeMiddlePointStrategy(self.data)
		if float(x) == self.data.x[0]:
			if float(y) in self.data.y:
				return LeftEdgeVertexStrategy(self.data)
			else:
				return LeftEdgeMiddlePointStrategy(self.data)
		if float(x) == self.data.x[-1]:
			if float(y) in self.data.y:
				return RightEdgeVertexStrategy(self.data)
			else:
				return RightEdgeMiddlePointStrategy(self.data)
					
	def getStrategy(self, x, y):
		if not self.data.isPointInData(x,y):
			print("Error: Provided data point not in data set")
			exit()
		if self.data.isPointOnEdge(x,y):
			return self.determineEdgeStrategy(x, y)
				
		if x in self.data.x and y in self.data.y:
			return VertexStrategy(self.data)
		
		if x in self.data.x and y not in self.data.y:
			return VerticalVertexMiddlePointStrategy(self.data)
			
		if x not in self.data.x and y in self.data.y:
			return HorizontalVertexMiddlePointStrategy(self.data)
		
		return MiddlePointStrategy(self.data)

# Abstract strategy class used as the parent for the Strategy Design Pattern
class Strategy:
	def __init__(self, function_data):
		self.strategyName = "Strategy"
		self.function_values = function_data
		
	def getStrategyName(self):
		return self.strategyName
		
	def determineVertecies(self, x, y):	
		raise NotImplementedError("This is an abstract class")
		
	def findValue(self, x, y):
		raise NotImplementedError("This is an abstract class")	
	


# The given point is in-between vertices (case A)
class MiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Middle Point Strategy"
	
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex
		left_vertex -=1 
		
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex 
		bottom_vertex -= 1
		
		top_left = self.function_values.data[left_vertex][top_vertex]
		bottom_left = self.function_values.data[left_vertex][bottom_vertex]
		top_right = self.function_values.data[right_vertex][top_vertex]
		bottom_right = self.function_values.data[right_vertex][bottom_vertex]
		
		return (top_left, bottom_left, top_right, bottom_right)
		
	
	def findValue(self, h, x, y):
		(top_left, bottom_left, top_right, bottom_right) = self.determineVertecies(x,y)
		dfdx = (top_right - top_left + bottom_right - bottom_left)/(4.0*self.function_values.deltaX)
		dfdy = (top_left - bottom_left + top_right - bottom_right)/(4.0*self.function_values.deltaY)
		
		print("dfdx = " + str(dfdx))
		print("dfdy = " + str(dfdy))
		
		x_new1 = x + h*(dfdy)
		y_new1 = y + h*(-dfdx)
		x_new2 = x + h*(-dfdy)
		y_new2 = y + h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is a vertex in the data grid (case B)
class VertexStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Vertex Strategy"
	
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex+1
		left_vertex -= 1
		
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex + 1
		bottom_vertex -= 1
		
		top = self.function_values.data[left_vertex+1][top_vertex]
		bottom = self.function_values.data[left_vertex+1][bottom_vertex]
		left = self.function_values.data[left_vertex][top_vertex-1]
		right = self.function_values.data[right_vertex][top_vertex-1]
		
		return (top, bottom, left, right)
	
	def findValue(self, h, x, y):
		(top, bottom, left, right) = self.determineVertecies(x, y)
		
		dfdx = (right - left)/(2.0*self.function_values.deltaX)
		dfdy = (top - bottom)/(2.0*self.function_values.deltaY)
		
		x_new1 = x + h*(dfdy)
		y_new1 = y + h*(-dfdx)
		x_new2 = x + h*(-dfdy)
		y_new2 = y + h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is on a vertex on the bottom edge (case C)
class BottomEdgeVertexStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Bottom Edge Vertex Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex + 1
		left_vertex -= 1
		
		top_vertex = 1
		center_vertex = 0
		
		top = self.function_values.data[left_vertex+1][top_vertex]
		center = self.function_values.data[left_vertex+1][center_vertex]
		left = self.function_values.data[left_vertex][center_vertex]
		right = self.function_values.data[right_vertex][center_vertex]
		
		return (top, center, left, right)
		
	def findValue(self, h, x, y):
		(top, center, left, right) = self.determineVertecies(x,y)
		
		dfdx = (right - left)/(2*self.function_values.deltaX)
		dfdy = (top - center)/(self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is on a vertex on the left edge (case D)
class LeftEdgeVertexStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Left Edge Vertex Strategy"
		
	def determineVertecies(self, x, y):
		center_vertex = 0
		right_vertex = center_vertex + 1
		
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex + 1
		bottom_vertex -= 1
		
		top = self.function_values.data[center_vertex][top_vertex]
		bottom = self.function_values.data[center_vertex][bottom_vertex]
		center = self.function_values.data[center_vertex][bottom_vertex+1]
		right = self.function_values.data[right_vertex][bottom_vertex+1]
		
		return (top, bottom, center, right)
		
	def findValue(self, h, x, y):
		(top, bottom, center, right) = self.determineVertecies(x,y)
		
		dfdx = (right - center)/(self.function_values.deltaX)
		dfdy = (top - bottom)/(2.0*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)
	
# The given point is on a vertex on the left edge (case E)	
class RightEdgeVertexStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Right Edge Vertex Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = -2
		center_vertex = -1
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex + 1
		bottom_vertex -= 1
		
		top = self.function_values.data[center_vertex][top_vertex]
		bottom = self.function_values.data[center_vertex][bottom_vertex]
		center = self.function_values.data[center_vertex][bottom_vertex+1]
		left = self.function_values.data[left_vertex][bottom_vertex+1]
		
		return (top, bottom, center, left)
		
	def findValue(self, h, x, y):
		(top, bottom, center, left) = self.determineVertecies(x,y)
		
		dfdx = (center - left)/(self.function_values.deltaX)
		dfdy = (top - bottom)/(2.0*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is on a vertex on the left edge (case F)
class TopEdgeVertexStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Top Edge Vertex Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex + 1
		left_vertex -= 1
		
		center_vertex = -1
		bottom_vertex = -2
		
		left = self.function_values.data[left_vertex][center_vertex]
		right = self.function_values.data[right_vertex][center_vertex]
		center = self.function_values.data[left_vertex+1][center_vertex]
		bottom = self.function_values.data[left_vertex+1][bottom_vertex]
		
		return (left, right, center, bottom)
		
	def findValue(self, h, x, y):
		(left, right ,center, bottom) = self.determineVertecies(x,y)
		
		dfdx = (right - left)/(2.0*self.function_values.deltaX)
		dfdy = (center - bottom)/(self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is in-between vertices on the bottom edge (case G)
class BottomEdgeMiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Bottom Edge Middle Point Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex
		left_vertex -= 1
		
		bottom_vertex = 0
		top_vertex = 1
		
		top_left = self.function_values.data[left_vertex][top_vertex]
		left = self.function_values.data[left_vertex][bottom_vertex]
		top_right = self.function_values.data[right_vertex][top_vertex]
		right = self.function_values.data[right_vertex][bottom_vertex]
		
		return (top_left, left, top_right, right)
	
	def findValue(self, h, x, y):
		(top_left, left, top_right, right) = self.determineVertecies(x,y)
		
		dfdx = (left - right)/(2*self.function_values.deltaX)
		dfdy = (top_left - left + top_right - right)/(2*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is in-between vertices on the left edge (case H)
class LeftEdgeMiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Left Edge Middle Point Strategy"
		
	def determineVertecies(self, x, y):
		center_vertex = 0
		right_vertex = 1
		
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex
		bottom_vertex -= 1
		
		top = self.function_values.data[center_vertex][top_vertex]
		bottom = self.function_values.data[center_vertex][bottom_vertex]
		top_right = self.function_values.data[right_vertex][top_vertex]
		bottom_right = self.function_values.data[right_vertex][bottom_vertex]
		
		return (top, bottom, top_right, bottom_right)
	
	def findValue(self, h, x, y):
		(top, bottom, top_right, bottom_right) = self.determineVertecies(x,y)
		
		dfdx = (top_right - top + bottom_right - bottom)/(2.0*self.function_values.deltaX)
		dfdy = (top - bottom)/(2*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is in-between vertices on the right edge (case I)
class RightEdgeMiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Right Edge Middle Point Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = -2
		right_vertex = -1
		
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex
		bottom_vertex -= 1
		
		top_left = self.function_values.data[left_vertex][top_vertex]
		top = self.function_values.data[right_vertex][top_vertex]
		bottom_left = self.function_values.data[left_vertex][bottom_vertex]
		bottom = self.function_values.data[right_vertex][bottom_vertex]
		
		return (top_left, top, bottom_left, bottom)
	
	def findValue(self, h, x, y):
		(top_left, top, bottom_left, bottom) = self.determineVertecies(x,y)
		
		dfdx = (top - top_left + bottom - bottom_left)/(2.0*self.function_values.deltaX)
		dfdy = (top - bottom)/(2.0*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is in-between vertices on the top edge (case J)
class TopEdgeMiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Top Edge Middle Point Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex
		left_vertex -= 1
		
		top_vertex = -1
		bottom_vertex = -2
		
		left = self.function_values.data[left_vertex][top_vertex]
		bottom_left = self.function_values.data[left_vertex][bottom_vertex]
		right = self.function_values.data[right_vertex][top_vertex]
		bottom_right = self.function_values.data[right_vertex][bottom_vertex]
		
		return (left, bottom_left, right, bottom_right)
	
	def findValue(self, h, x, y):
		(left, bottom_left, right, bottom_right) = self.determineVertecies(x,y)
		
		dfdx = (left - right)/(2.0*self.function_values.deltaX)
		dfdy = (left - bottom_left + right - bottom_right)/(2.0*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is in-between vertices on a horizontal edge 
# (not in middle of 4 vertices) (case K)
class HorizontalVertexMiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Horizontal Vertex Middle Point Strategy"
		
	def determineVertecies(self, x, y):
		left_vertex = 0
		while self.function_values.x[left_vertex] < x:
			left_vertex += 1
		right_vertex = left_vertex
		left_vertex -= 1
		
		center_vertex = 0
		while self.function_values.y[center_vertex] != y:
			center_vertex += 1
		bottom_vertex = center_vertex - 1
		top_vertex = center_vertex + 1
		
		top_left = self.function_values.data[left_vertex][top_vertex]
		left = self.function_values.data[left_vertex][center_vertex]
		bottom_left = self.function_values.data[left_vertex][bottom_vertex]
		top_right = self.function_values.data[right_vertex][top_vertex]
		right = self.function_values.data[right_vertex][center_vertex]
		bottom_right = self.function_values.data[right_vertex][bottom_vertex]
		
		return (top_left, left, bottom_left, top_right, right, bottom_right)
	
	def findValue(self, h, x, y):
		(top_left, left, bottom_left, top_right, right, bottom_right) = self.determineVertecies(x,y)
		
		dfdx = (right - left)/(2.0*self.function_values.deltaX)
		dfdy = (top_left - bottom_left + top_right - bottom_right)/(4.0*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)

# The given point is in-between vertices on a horizontal edge 
# (not in middle of 4 vertices) (case L)
class VerticalVertexMiddlePointStrategy(Strategy):
	def __init__(self, function_data):
		Strategy.__init__(self, function_data)
		self.strategyName = "Vertical Vertex Middle Point Strategy"
	
	def determineVertecies(self, x, y):
		center_vertex = 0
		while self.function_values.x[center_vertex] != x:
			center_vertex += 1
		left_vertex = center_vertex - 1
		right_vertex = center_vertex + 1
		
		bottom_vertex = 0
		while self.function_values.y[bottom_vertex] < y:
			bottom_vertex += 1
		top_vertex = bottom_vertex
		bottom_vertex -= 1
		
		top_left = self.function_values.data[left_vertex][top_vertex]
		bottom_left = self.function_values.data[left_vertex][bottom_vertex]
		top = self.function_values.data[center_vertex][top_vertex]
		bottom = self.function_values.data[center_vertex][bottom_vertex]
		top_right = self.function_values.data[right_vertex][top_vertex]
		bottom_right = self.function_values.data[right_vertex][bottom_vertex]
		
		return (top_left, bottom_left, top, bottom, top_right, bottom_right)
		
	def findValue(self, h, x, y):
		(top_left, bottom_left, top, bottom, top_right, bottom_right) = self.determineVertecies(x,y)
		
		dfdx = (top_right - top_left + bottom_right - bottom_left)/(4.0*self.function_values.deltaX)
		dfdy = (top - bottom)/(2.0*self.function_values.deltaY)
		
		x_new1 = x + 0.5*h*(dfdy)
		y_new1 = y + 0.5*h*(-dfdx)
		x_new2 = x + 0.5*h*(-dfdy)
		y_new2 = y + 0.5*h*(dfdx)
		
		return (x_new1, y_new1, x_new2, y_new2)
		
# Read in the data from the given file
class FunctionData:
	def readDataFromFile(self, filename):
		data_file = open(filename, "r")
		data_string = data_file.readlines()
		self.data_length = len(data_string)-1
		
		self.x = []
		self.y = []
		self.f = []
		self.data = [[0]* (int) (math.sqrt(self.data_length)) for _ in xrange((int) (math.sqrt(self.data_length)))]
		
		for i in range(1, len(data_string)):
			temp_string = data_string[i]
			temp_string = temp_string.strip()
			temp_string = temp_string.split(",")
			
			if float(temp_string[0]) not in self.x:
				self.x.append(float(temp_string[0]))
			if float(temp_string[1]) not in self.y:
				self.y.append(float(temp_string[1]))
			self.f.append(float(temp_string[2]))
			
		self.deltaX = self.x[1] - self.x[0]
		self.deltaY = self.y[1] - self.y[0]
			
		for i in range(0, len(self.data)):
			for j in range(0, len(self.data)):
				self.data[i][j] = self.f[i*len(self.data) + j]
		
		data_file.close()
	
	def isPointInData(self, x, y):
		if float(x) > max(self.x) or float(x) < min(self.x):
			return False
		if float(y) > max(self.y) or float(y) < min(self.y):
			return False
		return True
	
	def isPointOnEdge(self, x, y):
		if float(x) == self.x[0] or float(x) == self.x[-1]:
			return True
		if float(y) == self.y[0] or float(y) == self.y[-1]:
			return True
		return False

def testStrategyAssignment(function_data):
	factory = StrategyFactory(function_data)
	strategy1 = factory.getStrategy(1.51,1.51)
	strategy2 = factory.getStrategy(1.5, 1.5)
	strategy3 = factory.getStrategy(1.5, 0.0)
	strategy4 = factory.getStrategy(0.0, 1.5)
	strategy5 = factory.getStrategy(1.9, 1.5)
	strategy6 = factory.getStrategy(1.5, 1.9)
	strategy7 = factory.getStrategy(1.51, 0.0)
	strategy8 = factory.getStrategy(0.0, 1.51)
	strategy9 = factory.getStrategy(1.9, 1.51)
	strategy10 = factory.getStrategy(1.51, 1.9)
	strategy11 = factory.getStrategy(1.51, 1.5)
	strategy12 = factory.getStrategy(1.5, 1.51)
	
	if strategy1.getStrategyName() != "Middle Point Strategy":
		print("(1.51, 1.51) should have a Middle Point Strategy")
		
	if strategy2.getStrategyName() != "Vertex Strategy":
		print("(1.5, 1.5) should have a Vertex Strategy")
	
	if strategy3.getStrategyName() != "Bottom Edge Vertex Strategy":
		print("(1.5, 0) should have a Bottom Edge Vertex Strategy")
		
	if strategy4.getStrategyName() != "Left Edge Vertex Strategy":
		print("(0, 1.5) should have a Left Edge Vertex Strategy")
	
	if strategy5.getStrategyName() != "Right Edge Vertex Strategy":
		print("(1.9, 1.5) should have a Right Edge Vertex Strategy")
	
	if strategy6.getStrategyName() != "Top Edge Vertex Strategy":
		print("(1.5, 1.9) should have a Top Edge Vertex Strategy")
	
	if strategy7.getStrategyName() != "Bottom Edge Middle Point Strategy":
		print("(1.51, 0) should have a Bottom Edge Middle Point Strategy")
		
	if strategy8.getStrategyName() != "Left Edge Middle Point Strategy":
		print("(0, 1.51) should have a Left Edge Middle Point Strategy")
		
	if strategy9.getStrategyName() != "Right Edge Middle Point Strategy":
		print("(1.9, 1.51) should have a Right Edge Middle Point Strategy")
		
	if strategy10.getStrategyName() != "Top Edge Middle Point Strategy":
		print("(1.51, 1.9) should have a Top Edge Middle Point Strategy")
		
	if strategy11.getStrategyName() != "Horizontal Vertex Middle Point Strategy":
		print("(1.51, 1.5) should have a Horizontal Vertex Middle Point Strategy")
	
	if strategy12.getStrategyName() != "Vertical Vertex Middle Point Strategy":
		print("(1.5, 1.51) should have a Vertical Vertex Middle Point Strategy")


def main():
	x = 0.5
	y = 0.75
	h = 0.005
	#print("(x,y) = (" + str(x) + "," + str(y) + ")")
	filename = sys.argv[1]

	function_data = FunctionData()
	function_data.readDataFromFile(sys.argv[1])
	
	testStrategyAssignment(function_data)
	
	level_curve_solver = LevelCurvesSolver(h, x, y, function_data)
	#level_curve_solver.plotLevelCurve()
main()
