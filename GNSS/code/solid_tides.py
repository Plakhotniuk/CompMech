import matplotlib.pyplot as plt
import numpy as np


plt.figure()

data = np.loadtxt(f'../data/SolidTides24Days.txt')

times = data[:, 4] / 3600  # time in h
tides = data[:, 3]

plotData1 = np.sqrt(data[:, 1]**2 + data[:, 2]**2)
plt.plot(times, tides, label=r'measurements')


plt.grid()
plt.xlabel('time, h')
plt.ylabel('Tides, m')
plt.legend()
plt.savefig('../plots/solid_tides24days.png', dpi=1000)
plt.show()