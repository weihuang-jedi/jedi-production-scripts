import numpy as np
import matplotlib.pyplot as plt

import tkinter
import matplotlib
matplotlib.use('TkAgg')

np.random.seed(1)
data = np.random.random(100) * 100
bins = np.linspace(0, 100, 10)

histogram = np.zeros_like(bins)

bin_indexes = np.searchsorted(bins, data)
np.add.at(histogram, bin_indexes, 1)

print('bin_indexes = ', bin_indexes)
print('histogram = ', histogram)

plt.bar(bins, histogram)
plt.show()

plt.hist(histogram, bins=10, color='blue')
plt.show()

