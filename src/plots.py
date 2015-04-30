# Add plots here

import matplotlib.pyplot as plt
import numpy as np

def plot_p_vals(p_vals):
    plt.xscale('log')
    plt.hist(p_vals, bins=np.logspace(-3.0, 1.0, 20))
    plt.show()
