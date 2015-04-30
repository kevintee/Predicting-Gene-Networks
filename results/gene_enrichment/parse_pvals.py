import numpy as np
import matplotlib.pyplot as plt
import pdb
import math

# plots pvals on histogram
def main():
	# change fname for other file of p-value
	f = open('merlin_p_vals.txt', 'rb')
	p_vals = f.read()
	f.close()

	splitted = p_vals.split('\n')
	numbers = []
	for item in splitted:
		ooo = item.split('E')
		print ooo
		if len(ooo) == 2:
			front, back = ooo[0], ooo[1]
			number = float(front) * 10**(int(back))
			print number
			print number == 0.0
			if number == 0.0:
				numbers.append(0.0)
			else:
				numbers.append(math.log(number, 10))

	plt.hist(numbers, bins=[-50, -20, -10, -5, -4, -3, -2, -1, 0])
	plt.title("MERLIN: P-values")
	plt.xlabel("P-Value")
	plt.ylabel("Frequency")
	plt.show()
	pdb.set_trace()
if __name__ == '__main__':
	main()