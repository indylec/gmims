import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

coeffs=ascii.read('coeffs.txt')

a=coeffs['a']
b=coeffs['b']
a_inv=coeffs['a_inv']
b_inv=coeffs['b_inv']

intercept_1=-b/a
intercept_2=-b_inv/a_inv

plt.subplot(211)
plt.plot(b,'k-+',label='Direct intercept')
plt.plot(intercept_2,'r--^',label='Intercept from inverted plot')
plt.title('GMIMS vs Haslam')
#plt.xlabel('HEALPix pixel number')
plt.ylabel('y-axis intercept (K)')
plt.legend()

plt.subplot(212)
plt.plot(b_inv,'k-+',label='Direct intercept')
plt.plot(intercept_1,'r--^',label='Intercept from inverted plot')
plt.title('Haslam vs GMIMS')
plt.xlabel('HEALPix pixel number')
plt.ylabel('y-axis intercept (K)')
plt.legend()

plt.show()
