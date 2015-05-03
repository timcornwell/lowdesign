from telopt import *
import numpy
		
		
hiono=300.0
waves=[6.0, 3.0, 2.0]
base=numpy.arange(0.01, 10.0, 0.01)
ti=TelIono()

phase=ti.ionosphere(base)
plt.text(1.1*base[0], 1.1*phase[0], 'Ionosphere')
plt.loglog(base, phase, color='black')
plt.title(r'Ionosphere')
plt.xlabel(r'$km^{-1}$')
plt.ylabel('Phase')
plt.axes().set_xlim([1e-2, 10.0])
plt.axes().set_ylim([1e-3, 1e1])
offset=1.1e-3
ind=0
color=['r', 'g', 'b']
for wave in waves:
	plt.axvline(1.0/numpy.sqrt(hiono*wave/1000.0), ls='--', color=color[ind])
	plt.text(1.05/numpy.sqrt(hiono*wave/1000.0), offset, '%.0fm' % wave, color=color[ind])
	offset=1.3*offset
	ind=ind+1

plt.axvline(1.0/35.0, color='purple', ls='--')
plt.text(1.1/35.0, 1.1e-3, 'Station beam')
	
plt.savefig('ionosphere.pdf')
