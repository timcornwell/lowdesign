from telopt import *
import numpy
		
		
hiono=300.0
waves=[6.0, 3.0, 2.0]
base=numpy.arange(0.01, 100.0, 0.01)
ti=TelIono()

plt.title(r'Ionosphere (r:50MHz, g:100MHz, b:150MHz)')
plt.xlabel('Distance (km)')
plt.ylabel('Phase')
plt.text(1.05*numpy.sqrt(hiono*6.0/1000.0), 1.1e-3, 'Fresnel zone')
plt.axes().set_xlim([1e-2, 100.0])
plt.axes().set_ylim([1e-3, 1e2])
ind=0
color=['r', 'g', 'b']
for wave in waves:
	phase=ti.ionosphere(base)*wave
	plt.loglog(base, phase, color=color[ind])
	plt.axvline(numpy.sqrt(hiono*wave/1000.0), ls='--', color=color[ind])
	ind=ind+1

plt.axvline(0.035, color='purple', ls='--')
plt.text(0.035, 1.1e-3, 'Station beam')
	
plt.savefig('ionosphere.pdf')

plt.clf()
base=numpy.arange(0.01, 100.0, 0.01)
veliono=500.0/3600.0
timeiono=base/veliono
plt.loglog(base, timeiono, color='black')
plt.title(r'Ionospheric motion')
plt.xlabel(r'Distance (km)')
plt.ylabel('Transit time (s)')

plt.axvline(0.035, color='red', ls='--')
plt.text(0.9*0.035, 1.1e-2, 'Station', horizontalalignment='right')

plt.axvline(6, color='green', ls='--')
plt.text(0.9*6.0, 1.1e-2, 'Core', horizontalalignment='right')

plt.axvline(80, color='blue', ls='--')
plt.text(0.9*80, 1.1e-2, 'Array', horizontalalignment='right')
	
plt.savefig('timescales.pdf')
