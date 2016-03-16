import numpy as np
import math
import scipy
import matplotlib
import random
import matplotlib.pyplot as plt
import numpy as np
import random
import copy

class confusion:

	def sigmac(self, freq=1.0e9, B=35.0):
 		theta=180.0*3600.0*3.0e8/(freq*B*1000.0*np.pi)
		return 1.2e-6 * np.power(freq / 3.014e9, -0.7) * np.power(theta/8.0, 10.0/3.0)

freqs=np.array([5.0e7, 1.0e8, 1.5e8])		
plt.clf()
baselines=np.power(10.0, np.arange(-1.0, 2.0, 0.1))
for freq in freqs:
	conf=confusion().sigmac(freq=freq, B=baselines)
	plt.loglog(baselines, conf)
plt.xlabel('Baseline length (km)')
plt.ylabel('Confusion noise (Jy)')
plt.title('Confusion noise at %s (MHz)' % (freqs*1e-6))
plt.savefig('confusion.pdf')
