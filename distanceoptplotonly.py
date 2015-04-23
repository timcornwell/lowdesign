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
nsource=600

# Observing wavelength
wave=3.0

# Diameter of station
diameter=35.0

# Baseline length (radius)
rbase=hiono*wave/diameter

# Radius of pierce points in the ionosphere
rpierce=rbase

# Number of source trials
nsourcetrials=1

# Number of array trials
ntrials=1000000
##########################################################################################

slowbd=numpy.zeros(nnoll)

# Find the SVD for the Baseline Design. We use this as the reference.	
lowbd=TelArray()
lowbd.readCSV('LOWBD', rcore=2.0, recenter=True, l1def='LOW_BD.csv')
slowbd=numpy.zeros(nnoll)
ts=TelSources()
tp=TelPiercings()
for sourcetrial in range(nsourcetrials):
	ts.construct(nsources=nsource, radius=wave/35.0)
	print "Number of sources ", nsource, "Trial ", sourcetrial
	tp.construct(ts,lowbd,hiono=hiono,rmin=1.8)
	slowbd=slowbd+(1.0/float(nsourcetrials))*tp.assess(nnoll=nnoll, rmax=rpierce, doplot=False)
	if sourcetrial == 0:
		tp.plot(rmax=rpierce)
		
print "Zero singular value %.3f" % slowbd[0]

# Loop over all trial values for the number of stations
srb={}
for stations in [11, 21, 31, 41, 51]:
	sname='LOW_RANDOMBOOLARDY_WT%d_wave%.1f_diam%.1f' % (stations, wave, diameter)

# Recalculate all dependent values
# 	diameter=35.0*numpy.sqrt(46.0/float(stations))
	print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
	rbase=hiono*wave/diameter
	rpierce=rbase
	print "Number of stations = %d, baseline radius = %.3f (km), pierce radius= %.3f (km)" % (stations, rbase, rpierce)
	ts=TelSources()
	ts.construct(nsources=nsource, radius=wave/diameter)

	lowrandbool=TelArray()
	lowrandbool.randomBoolardy(sname, nstations=stations, nhalo=stations, rhalo=rbase, diameter=diameter)
	trialconfig=copy.deepcopy(lowrandbool)
	best=copy.deepcopy(lowrandbool)

	bestDistance=best.excessDistance(False)
	T=0.1*bestDistance
	print "Initial distance metric = %.3f (km)" % bestDistance
	tp=TelPiercings()
	tp.construct(ts,best,hiono=hiono,rmin=1.8)

	trial=1
	while trial < ntrials:
		trialconfig.stations=best.stations.copy()
		if fresh:
			trialconfig=TelArray()
			trialconfig.randomBoolardy(sname, nstations=stations, nhalo=stations, rhalo=rbase)
		else:
			trialconfig.shakehalo(best.distance()/2.0)
		distance=trialconfig.excessDistance(False)
		prob=math.exp(-(bestDistance-distance)/T)
		rnd=random.uniform(0.0,1.0)
		if (distance>bestDistance):
			print "Trial %d: Selecting better config : %.3f -> %.3f" % (trial, bestDistance, distance)
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
			T=T*0.99		
		elif (distance<bestDistance) and (rnd<prob):
			print "Trial %d: Accepting worse config  : %.3f -> %.3f" % (trial, bestDistance, distance)
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
		trial=trial+1
	srb['%d'%stations]=numpy.zeros(nnoll)
	for sourcetrial in range(nsourcetrials):
		ts.construct(nsources=nsource, radius=wave/diameter)
		print "Number of sources ", nsource, "Trial ", sourcetrial
		tp.construct(ts,best,hiono=hiono,rmin=1.8)
		srb['%d'%stations]=srb['%d'%stations]+(1.0/float(nsourcetrials))*tp.assess(nnoll=nnoll, rmax=rpierce, doplot=False)
		if sourcetrial == 0:
			tp.plot(rmax=rpierce)
		

	# Scale the x axis to be a pseudo wave number
	print "Zero Singular value %.3f" % srb['%d'%stations][0]
	print "Distance %.3f (km) MST %.3f (km)" % (best.distance(), best.mst(False))
	best.saveCSV('%s.csv' % sname)
	best.writeKML('%s.kml' % sname)
	best.plot(rbase, '%s_ARRAY.pdf' % sname)
	best.mst(True, '%s_MST.pdf' % sname)
	print " "
	
	
plt.clf()
plt.rc('axes', color_cycle=['r', 'g', 'b', 'cyan', 'purple'])
print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
rbase=hiono*wave/35.0
rpierce=rbase
scales=numpy.sqrt(numpy.arange(nnoll)/(rpierce*rpierce))
plt.loglog(scales, slowbd/slowbd[0], color='black')

for stations in [11, 21, 31, 41, 51]:
# 	diameter=35.0*numpy.sqrt(46.0/float(stations))
	print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
	rbase=hiono*wave/diameter
	rpierce=rbase
	scales=numpy.sqrt(numpy.arange(nnoll)/(rpierce*rpierce))
	plt.loglog(scales, srb['%d'%stations]/slowbd[0])

base=numpy.arange(0.01, 10.0, 0.01)
phase=numpy.power(14.0*base,-1.8/2.0)
plt.loglog(base, phase, color='black')
plt.title('SVD 11, 21, 31, 41, 51, BDv1, %.1f, %.1f' % (wave, diameter))
plt.xlabel(r'$km^{-1}$')
plt.ylabel('Singular value')
plt.axes().set_xlim([1e-2, 10.0])
plt.axes().set_ylim([1e-3, 1e1])
plt.axvline(1.0/numpy.sqrt(hiono*wave/1000.0), color='r', ls='--')
plt.axvline(1.0/rpierce, color='g', ls='--')
plt.savefig('SVD_OPT_wave%sm_%.1fm.pdf' % (wave, diameter))
