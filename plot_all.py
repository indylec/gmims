import numpy as np
import healpy as hp
import theil_sen as ts
import matplotlib.pyplot as plt

data=np.genfromtxt('offsets_slopes.txt')


x=np.arange(-10,10,0.001)

for i in range(192):
    plt.plot(x,data[i][1]*x+data[i][2])


plt.show()
