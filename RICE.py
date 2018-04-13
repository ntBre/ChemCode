import math

class RICE:

	def __init__(self, Ka, C_1_initial, C_2_initial = 0, C_3_initial = 0, v1 = 1, v2 = 1, v3 = 1):
		self.Ka = Ka
		self.C_1_initial = C_1_initial
		self.C_2_initial = C_2_initial
		self.C_3_initial = C_3_initial
		self.v1 = v1
		self.v2 = v2
		self.v3 = v3
		
	def get_Ka(self):
		return self.Ka
	
	def get_C_1_initial(self):
		return self.C_1_initial

	def get_C_2_initial(self):
		return self.C_2_initial

	def get_C_3_initial(self):
		return self.C_3_initial

	def get_table(self):
		return \
"######################################\n\
# R #   R  -->  #    P1    #  +  P2  #\n\
######################################\n\
# I #{:10.2f} #{:8.2f}  #{:8.2f} #\n\
######################################\n\
# C #    -{:d}x    #   +{:d}x    #   +{:d}x   #\n\
######################################\n\
# E #{:7.2f}-{:d}x #{:6.2f}+{:d}x #{:5.2f}+{:d}x #\n\
######################################".format(self.C_1_initial, self.C_2_initial, self.C_3_initial,self.v1,self.v2, self.v3, self.C_1_initial, self.v1, self.C_2_initial, self.v2, self.C_3_initial, self.v3)

	def get_equation(self):
		return str(self.Ka) + " = ({:.2f}+{:d}x)({:.2f}+{:d}x)/({:.2f}-{:d}x)".format(self.C_2_initial, self.v2, self.C_3_initial, self.v3, self.C_1_initial, self.v1)

	def get_x(self):
		x_not = 1
		def function(x):
			f = self.C_2_initial*self.C_3_initial + self.C_2_initial*self.v3*x + self.v2*self.C_3_initial*x + self.v2*self.v3*x**2-self.Ka*self.C_1_initial + self.Ka*self.v1*x
			return f
		def fprime(x):
			f1 = self.C_2_initial*self.v3 + self.C_3_initial*self.v2 + 2*self.v2*self.v3*x+self.Ka*self.v1
			return f1
		x_one = 1 - (function(1)/fprime(1))
		while function(x_one) > 0.000000000000001:
			x_one = x_not - (function(x_not)/fprime(x_not))
			x_not = x_one
		return x_one
	
	def get_pX(self):
		return -math.log10(self.get_x())

	def set_Ka(self, Ka):
		self.Ka = Ka

	def set_C_1_initial(self, C_1_initial):
		self.C_1_initial = C_1_initial

	def set_C_2_initial(self, C_2_initial):
		self.C_2_initial = C_2_initial

	def set_C_3_initial(self, C_3_initial):
		self.C_3_initial = C_3_initial


#test = RICE(0.000000000556, 0.14)
#print test.get_table()
#print test.get_equation()
#print test.get_x()
#print test.get_pX()

def main():
	params = input("Enter your Ka, and the initial concentration of your reactant separated by commas. Optional arguments include initial concentrations of two products, and the stoichiometric coefficients of all three compounds, in the order Ka, C1, C2, C3, v1, v2, v3: ")
	call = RICE(*params)
	print "RICE table:"
	print call.get_table()
	print "Evaluated equation: ", call.get_equation()
	print "X = ", call.get_x()
	print "pX = ", call.get_pX()
main()
