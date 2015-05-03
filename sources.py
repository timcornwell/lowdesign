import numpy as numpy
import math
import scipy
import matplotlib
import random
import matplotlib.pyplot as plt
import numpy as np
import random
import copy
from telopt import *


freqs=np.array([5.0e7, 1.0e8, 1.5e8])		

plt.clf()
baselines=np.power(10.0, np.arange(-2.0, 3.0, 0.1))
for freq in freqs:
	conf=sources().confusion(freq=freq, B=baselines)
	plt.loglog(baselines, conf)
plt.xlabel('Baseline length (km)')
plt.ylabel('Noise (Jy)')
plt.title('Noise at 50 (r), 100 (g), 150 (b) (MHz) (solid: confusion, dotted: thermal)')
plt.axvline(0.035, color='purple', ls='--')
plt.text(0.035, 1.5e-9, 'Station diameter')
plt.axvline(80.0, color='green', ls='--')
plt.text(80.0, 1.5e-9, 'Array diameter')
plt.axhline(sources().tnoise(5e7), color='red', ls=':')
plt.axhline(sources().tnoise(1e8), color='green', ls=':')

plt.axhline(sources().tnoise(5e7, 300.0), color='red', ls='--')
plt.axhline(sources().tnoise(1e8, 300.0), color='green', ls='--')
plt.savefig('confusion.pdf')

plt.clf()
fluxes=np.power(10, np.arange(-6.0, 1.0, 0.1))
for freq in freqs:
	integratedflux=sources().integratedflux(fluxes, freq=freq)
	omega=numpy.power(3.0e8/(2*35.0*freq), 2)
	plt.loglog(fluxes, omega*integratedflux)
plt.xlabel('Flux (Jy)')
# plt.axes().set_xlim([1e-6, 1e2])
# plt.axes().set_ylim([1, 1e3])
plt.axvline(sources().tnoise(5e7, 300.0)*numpy.sqrt(1024.0), color='red', ls='--')
plt.axvline(sources().tnoise(1e8, 300.0)*numpy.sqrt(1024.0), color='green', ls='--')
plt.axvline(sources().tnoise(5e7), color='red', ls=':')
plt.axvline(sources().tnoise(1e8), color='green', ls=':')

plt.ylabel(r'Integrated Flux per station beam')
plt.title('Integrated flux per station beam at 50 (r), 100 (g), 150 (b) (MHz)')
plt.savefig('integratedflux.pdf')

plt.clf()
for freq in freqs:
	omega=numpy.power(3.0e8/(2*35.0*freq), 2)
	numbers=omega*sources().numbers(fluxes, freq=freq)
	plt.loglog(fluxes, numbers)
plt.axes().set_ylim([1.0, 1e6])
# plt.axes().set_xlim([1e-6, 1e2])
plt.axvline(sources().tnoise(5e7, 300.0)*numpy.sqrt(1024.0), color='red', ls='--')
plt.axvline(sources().tnoise(1e8, 300.0)*numpy.sqrt(1024.0), color='green', ls='--')
plt.axvline(sources().tnoise(5e7), color='red', ls=':')
plt.axvline(sources().tnoise(1e8), color='green', ls=':')

plt.ylabel('Numbers of sources per station beam')
plt.xlabel('Flux (Jy)')
plt.title('Numbers of sources per station beam at 50 (r), 100 (g), 150 (b) (MHz)')
plt.savefig('numbers.pdf')

plt.clf()
for freq in freqs:
	omega=numpy.power(3.0e8/(2*35.0*freq), 2)
	numbers=omega*sources().numbers(fluxes, freq=freq)
	integratedflux=omega*sources().integratedflux(fluxes, freq=freq)
	plt.loglog(numbers, integratedflux)
plt.axhline(sources().tnoise(5e7, 300.0)*numpy.sqrt(1024.0), color='red', ls='--')
plt.axhline(sources().tnoise(1e8, 300.0)*numpy.sqrt(1024.0), color='green', ls='--')
# plt.axes().set_xlim([1.0, 1e8])
# plt.axes().set_ylim([1, 1e2])
plt.xlabel('Number of sources per station beam')
plt.ylabel('Flux per station beam')
plt.title('Flux vs Numbers of sources per station beam at 50 (r), 100 (g), 150 (b) (MHz)')
plt.savefig('numbersintegrated.pdf')