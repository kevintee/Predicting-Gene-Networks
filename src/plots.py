# Add plots here
import pylab
import matplotlib.pyplot as plt
import numpy as np

def plot_p_vals(p_vals, title):
    plt.xscale('log')
    plt.title(title)
    plt.ylabel('Number of Transcription Factors')
    plt.xlabel('p-value')
    plt.hist(p_vals, bins=np.logspace(-3.0, 1.0, 20))
    plt.show()


def plot_p_vals_both(pval1, pval2, title):
	ax = pylab.subplot(111)
	ax.set_xlim([-5, 1])
	ax.set_ylim([-5, 1])
	ax.set_title(title)
	ax.set_xlabel('p_value sbm')
	ax.set_ylabel('p_value merlin')
	ax.scatter(pval1, pval2)
	ax.plot([-5, 1], [-5, 1])


	plt.show()