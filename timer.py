from schur import *
import time
import matplotlib.pyplot as plt

# Function to generate all partitions of n
# from http://jeromekelleher.net/category/combinatorics.html        
def partitions(n):
	a = [0 for i in range(n+1)]
	k = 1
	y = n-1
	while k != 0:
		x = a[k-1] + 1
		k -= 1
		while 2*x <= y:
			a[k] = x
			y -= x
			k += 1
		l = k+1
		while x <= y:
			a[k] = x
			a[l] = y
			yield a[:k+2]
			x += 1
			y -= 1
		a[k] = x+y
		y = x+y-1
		yield a[:k+1]

def time_function(function, args, trials):
	times = []
	for i in range(trials):
		start = time.time()
		_ = function(*args)
		end = time.time()
		times.append(end-start)
	return sum(times)/len(times)

def plot_times():
	print("Plotting times for normal calculation")
	times_cp, length_cp, statements_cp, n_cp = [], [], [], []
	for N in range(1,13):
		for partition in partitions(N):
			print(partition)
			time_ = time_function(Schur, [partition[::-1], len(partition)], 10 if N < 10 else 1)
			times_cp.append(time_)
			length_cp.append(len(partition))
			statements_cp.append(count_statements(partition))
			n_cp.append(N)
	print("Plotting times for using Weyl formula")
	times_det, length_det, n_det = [], [], []
	for N in range(1,20):
		for partition in partitions(N):
			print(partition)
			time_ = time_function(Schur_weyl, [partition[::-1], len(partition)], 10 if N < 10 else 1)
			times_det.append(time_)
			length_det.append(len(partition))
			n_det.append(N)
	print("Plotting times with simplification")
	times_sympy, length_sympy, n_sympy = [], [], []
	for N in range(1,6):
		for partition in partitions(N):
			print(partition)
			time_ = time_function(Schur_weyl_sympy, [partition[::-1], len(partition)], 1)
			times_sympy.append(time_)
			length_sympy.append(len(partition))
			n_sympy.append(N)
	plt.scatter(length_cp, times_cp, c = "r", label = "Constraint programming")
	plt.scatter(length_det, times_det, c = "g", label = "Weyl's formula")
	plt.scatter(length_sympy, times_sympy, c = "b", label = "Weyl's formula + sympy")
	plt.title("Length of partition vs time (s)")
	plt.legend()
	plt.savefig("length_plot.png")
	plt.clf()
	plt.scatter(n_cp, times_cp, c = 'r', label = "Constraint programming")
	plt.scatter(n_det, times_det, c = 'g', label = "Weyl's formula")
	plt.scatter(n_sympy, times_sympy, c = 'b', label = "Weyl's formula + sympy")
	plt.title("Sum of partition vs time (s)")
	plt.legend()
	plt.savefig("n_plot.png")
	plt.clf()
