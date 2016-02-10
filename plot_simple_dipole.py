import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

g0_coeffs_tab=ascii.read('haslam_adjusted-g0g_coeffs_single_dipole.txt', guess=False, delimiter=' ')
g1_coeffs_tab=ascii.read('haslam_adjusted-g1g_coeffs_single_dipole.txt', guess=False, delimiter=' ')
g7_coeffs_tab=ascii.read('haslam_adjusted-g7g_coeffs_single_dipole.txt', guess=False, delimiter=' ')
g8_coeffs_tab=ascii.read('haslam_adjusted-g8g_coeffs_single_dipole.txt', guess=False, delimiter=' ')

offsets=[g0_coeffs_tab['x0'],g1_coeffs_tab['x0'],g7_coeffs_tab['x0'],g8_coeffs_tab['x0']]
err=[g0_coeffs_tab['x0err'],g1_coeffs_tab['x0err'],g7_coeffs_tab['x0err'],g8_coeffs_tab['x0err']]

freqs=[1300855232.,1349402512.,1640686192.,1689233472.]

plt.plot(freqs,offsets,'go',linestyle='')
plt.errorbar(freqs,offsets,yerr=err,ls='')
plt.xlim([1.2E9, 1.8E9])
plt.ylim([-0.15,0.05])
plt.xlabel('GMIMS frequency (Hz)')
plt.ylabel('Offset wrt Haslam (K)')
plt.show()
