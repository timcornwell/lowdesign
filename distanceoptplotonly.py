from telopt import *
import math
import random
import copy

#
# Radius=1/2 FWZ ~ FWHM
#
#random.seed(781450893)
fresh=False
##########################################################################################
# Height of ionosphere
hiono=300

# Number of zernike terms numbered using Noll notation
nnoll=200

# Number of sources that can be used for calibration by 35m diameter station
nsource=40

# Observing wavelength
wave=3.0

# Diameter of station
diameter=35.0

# Baseline length (radius)
rbase=40.0

# Radius of pierce points in the ionosphere
rpierce=2.0*hiono*wave/(diameter)

# Number of source trials
nsourcetrials=10

# Number of array trials
ntrials=1000000

# Number of stations
nstations = 45
##########################################################################################

slowbd=numpy.zeros(nnoll)

# Find the SVD for the Baseline Design. We use this as the reference.	
lowbd=TelArray()
lowbd.readCSV('LOW_BD', rcore=2.0, recenter=True, l1def='LOW_BD.csv')
slowbd=numpy.zeros(nnoll)
srb=numpy.zeros(nnoll)
ts=TelSources()
tp=TelPiercings()
for sourcetrial in range(nsourcetrials):
	ts.construct(nsources=nsource, radius=2.0*wave/35.0)
	print "Number of sources ", nsource, "Trial ", sourcetrial
	tp.construct(ts,lowbd,hiono=hiono,rmin=1.8)
	slowbd=slowbd+(1.0/float(nsourcetrials))*tp.assess(nnoll=nnoll, rmax=rpierce, doplot=False)
	if sourcetrial == 0:
		tp.plot(rmax=rpierce)
		
print "Singular value[0] %.3f" % slowbd[0]

rep=10
nconverged=0
while (rep > 0):

	sname='LOW_RB%d_source%d_wave%.1f_diam%.1f_rep%d' % (nsource, nstations, wave, diameter, rep)

# Recalculate all dependent values
	print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
	print "Number of stations = %d, baseline radius = %.3f (km), pierce radius= %.3f (km)" % (nstations, rbase, rpierce)
	ts=TelSources()
	ts.construct(nsources=nsource, radius=wave/diameter)

	lowrandbool=TelArray()
	lowrandbool.randomBoolardy(sname, nstations=nstations, nhalo=nstations, rhalo=rbase, diameter=diameter)
	trialconfig=copy.deepcopy(lowrandbool)
	best=copy.deepcopy(lowrandbool)

	bestDistance=best.excessDistance(False)
	T=bestDistance
	print "Initial distance metric = %.3f (km), T=%.3f (km)" % (bestDistance, T)
	tp=TelPiercings()
	tp.construct(ts,best,hiono=hiono,rmin=1.8)
	trial=1
	while (trial < ntrials) and (T > 0.01):
		trialconfig.stations=best.stations.copy()
		if bestDistance<0.001:
			print "Stalled - resetting"
			trialconfig=TelArray()
			trialconfig.randomBoolardy(sname, nstations=nstations, nhalo=nstations, rhalo=rbase)
		else:
			trialconfig.shakehalo(best.distance()/2.0)
		distance=trialconfig.excessDistance(False)
		prob=math.exp(-(bestDistance-distance)/T)
		rnd=random.uniform(0.0,1.0)
		if (distance>bestDistance):
			T=T*0.999		
			print "Trial %d: Selecting better config : %.3f / %.3f: T=%.3f" % (trial, bestDistance, distance, T)
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
		elif (distance<bestDistance) and (rnd<prob):
			print "Trial %d: Accepting worse config  : %.3f \ %.3f" % (trial, bestDistance, distance)
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
		trial=trial+1
		
	if (T < 0.01):
		nconverged=nconverged+1
		for sourcetrial in range(nsourcetrials):
			ts.construct(nsources=nsource, radius=2.0*wave/diameter)
			print "Number of sources ", nsource, "Trial ", sourcetrial
			tp.construct(ts,best,hiono=hiono,rmin=1.8)
			srb=srb+(1.0/float(nsourcetrials))*tp.assess(nnoll=nnoll, rmax=rpierce, doplot=False)
			if sourcetrial == 0:
				tp.plot(rmax=rpierce)
		

	# Scale the x axis to be a pseudo wave number
	print "Singular value[0] %.3f" % srb[0]
	print "Distance %.3f (km) MST %.3f (km)" % (best.distance(), best.mst(False))
	best.saveCSV('%s.csv' % sname)
	best.writeKML('%s.kml' % sname)
	best.plot(rbase, '%s_ARRAY.pdf' % sname)
	best.mst(True, '%s_MST.pdf' % sname)
	print " "
	rep=rep-1
	
	
plt.clf()
plt.rc('axes', color_cycle=['r', 'g', 'b', 'cyan', 'purple'])
scales=numpy.sqrt(numpy.arange(nnoll)/(rpierce*rpierce))
plt.loglog(scales, slowbd/slowbd[0], color='black')

print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
scales=numpy.sqrt(numpy.arange(nnoll)/(rpierce*rpierce))
plt.loglog(scales, srb/slowbd[0])

base=numpy.arange(0.01, 10.0, 0.01)
phase=numpy.power(14.0*base,-1.8/2.0)
plt.text(1.1*base[0], 1.1*phase[0], 'Ionosphere')
plt.loglog(base, phase, color='black')
plt.title(r'SVD 45 BDv1, $\lambda$%.1f (m), D%.1f (m)' % (wave, diameter))
plt.xlabel(r'$km^{-1}$')
plt.ylabel('Singular value/Phase')
plt.axes().set_xlim([1e-2, 10.0])
plt.axes().set_ylim([1e-3, 1e1])
plt.axvline(1.0/numpy.sqrt(hiono*wave/1000.0), color='r', ls='--')
plt.text(1.1/numpy.sqrt(hiono*wave/1000.0), 1.1e-3, 'Fresnel scale')
plt.axvline(1.0/rpierce, color='g', ls='--')
plt.text(1.1/rpierce, 1.1e-3, 'Station beam')
plt.savefig('SVD_OPT_wave%sm_%.1fm.pdf' % (wave, diameter))
