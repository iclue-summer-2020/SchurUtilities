from ortools.sat.python import cp_model
import copy
import sympy
from functools import reduce

class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
	def __init__(self, variables, young_diagram):
		cp_model.CpSolverSolutionCallback.__init__(self)
		self.variables = variables
		self.solution_count = 0
		self.young_diagram = young_diagram
		self.solutions = []

	def on_solution_callback(self):
		self.solution_count += 1
		current_solution = copy.deepcopy(self.young_diagram.diagram)
		solution = {}
		for v in self.variables:
			solution[v] = self.Value(v)
			coord = [int(i) for i in str(v).split(",")]
			current_solution[coord[0]][coord[1]] = int(self.Value(v))
		current_young_diagram = YoungDiagram(self.young_diagram.partition)
		current_young_diagram.label(current_solution)
		self.solutions.append(current_young_diagram)
	def solution_count(self):
		return self.solution_count


class YoungDiagram:
	def __init__(self, partition):
		self.partition = partition
		self.diagram = [[None]*i for i in partition]
	def label(self, label):
		if len(label) != len(self.diagram): raise ValueError("dimension mismatch!")
		for i in range(len(label)):
			if len(label[i]) != len(self.diagram[i]): raise ValueError("dimension mismatch!")
			for j in range(len(label[i])):
				self.diagram[i][j] = label[i][j]
	def __str__(self):
		diagram_string = [[str(j) if j else "x" for j in i] for i in self.diagram]
		return "\n".join([",".join(i) for i in diagram_string])

def get_all_semistandard_young_tableaux(young_diagram, max_label):
	model = cp_model.CpModel()
	variables = copy.deepcopy(young_diagram.diagram)
	for i in range(len(variables)):
		for j in range(len(variables[i])):
			variables[i][j] = model.NewIntVar(1, max_label, "{0},{1}".format(i,j))
	all_variables = sum(variables, [])
	padded_variables = copy.deepcopy(variables)
	length = max([len(i) for i in variables])
	for i in range(len(padded_variables)):
		while len(padded_variables[i]) != length:
			padded_variables[i].append(None)
	# Condition 1: Weakly increasing along rows
	for i in range(len(variables)):
		for j in range(len(variables[i])-1):
			model.Add(variables[i][j] <= variables[i][j+1])
	# Condition 2: Strictly increasing along columns
	transposed_variables = list(zip(*padded_variables))
	for i in range(length):
		column = transposed_variables[i]
		for j in range(len(column)-1):
			if column[j+1] == None: break
			model.Add(column[j] < column[j+1])
	solver = cp_model.CpSolver()
	solution_printer = VarArraySolutionPrinter(all_variables, young_diagram)
	status = solver.SearchForAllSolutions(model, solution_printer)
	return solution_printer.solutions

def count_statements(partition):
	out = 0
	variables = YoungDiagram(partition).diagram
	all_variables = sum(variables, [])
	padded_variables = copy.deepcopy(variables)
	length = max([len(i) for i in variables])
	for i in range(len(padded_variables)):
		while len(padded_variables[i]) != length:
			padded_variables[i].append(None)
	# Condition 1: Weakly increasing along rows
	for i in range(len(variables)):
		for j in range(len(variables[i])-1):
			out += 1
	# Condition 2: Strictly increasing along columns
	transposed_variables = list(zip(*padded_variables))
	for i in range(length):
		column = transposed_variables[i]
		for j in range(len(column)-1):
			if column[j+1] == None: break
			out += 1
	return out
def Schur(partition, N, variables = None):
	if variables == None: variables = [None] + list(sympy.symbols("x1:"+str(N+1)))
	out = 0
	for young_tableaux in get_all_semistandard_young_tableaux(YoungDiagram(partition), N):
		labels = sum(young_tableaux.diagram, [])
		out += reduce(lambda x,y: x*y, [variables[i]**(labels.count(i)) for i in range(1, N+1)])
	return out

def Schur_weyl(partition, N, variables = None):
	if variables == None: variables = [None] + list(sympy.symbols("x1:"+str(N+1)))
	numerator_matrix = [[None for i in range(N)] for j in range(N)]
	for i in range(1,N+1):
		for j in range(1,N+1):
			numerator_matrix[i-1][j-1] = variables[j] ** (partition[i-1] + N - i)
	denominator = 1
	for j in range(1,N+1):
		for i in range(1,j):
			denominator *= (variables[i] - variables[j])
	numerator = sympy.Matrix(numerator_matrix).det(method="LU")
	return numerator/denominator
	
def Schur_weyl_sympy(partition, N, variables = None):
	return sympy.factor(Schur_weyl(partition, N, variables))

def poly_equals(poly1, poly2):
	return sympy.expand(poly1) == sympy.expand(poly2)

def sanity_test():
	variables = [None] + list(sympy.symbols("x1 x2 x3"))
	_, x1, x2, x3 = variables
	assert Schur((2,1), 3, variables = variables) == x1*x2**2 + x2*x1**2 + x1*x3**2 + x3*x1**2 + x2*x3**2 + x3*x2**2 + 2*x1*x2*x3 # From Yong paper
	assert Schur((2,1,1), 3, variables = variables) == sympy.expand(x1*x2*x3*(x1+x2+x3)) # From Wikipedia article
	assert Schur((2,2,0), 3, variables = variables) == (x1*x2)**2 + (x1*x3)**2 + (x2*x3)**2 + x1**2*x2*x3 + x1*x2**2*x3 + x1*x2*x3**2 # From Wikipedia article
	assert poly_equals(Schur((2,1,1), 3), Schur_weyl((2,1,1), 3))
	assert poly_equals(Schur((4,2,1,1), 4), Schur_weyl((4,2,1,1),4))
	return "Passed"




