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
nsource=6

# Observing wavelength
wave=6.0

# Diameter of station
diameter=35.0

# Baseline length (radius)
rbase=hiono*wave/diameter

# Radius of pierce points in the ionosphere
rpierce=2.0*rbase

# Number of source trials
nsourcetrials=10

# Number of array trials
ntrials=1000000
##########################################################################################

slowbd=numpy.zeros(nnoll)

plt.clf()

# Find the SVD for the Baseline Design. We use this as the reference.	
lowbd=TelArray()
lowbd.readCSV('LOWBD', rcore=2.0, recenter=True, l1def='LOW_BD.csv')
slowbd=numpy.zeros(nnoll)
ts=TelSources()
tp=TelPiercings()
for sourcetrial in range(nsourcetrials):
	ts.construct(nsources=nsource, radius=wave/diameter)
	print "Number of sources ", nsource, "Trial ", sourcetrial
	tp.construct(ts,lowbd,hiono=300,rmin=1.8)
	slowbd=slowbd+(1.0/float(nsourcetrials))*tp.assess(nnoll=nnoll, rmax=rpierce, doplot=False)
# Scale the x axis to be a pseudo wave number
scales=numpy.sqrt(numpy.arange(nnoll)/(rpierce*rpierce))
plt.plot(scales, numpy.sqrt(slowbd/slowbd[0]), color='black')
print "Zero sqrt(singular value) %s" % numpy.sqrt(slowbd[0])

# Loop over all trial values for the number of stations
color={'11':'r', '21': 'g', '31':'b', '41':'yellow', '51':'cyan'}
for stations in [11, 21, 31, 41, 51]:

# Recalculate all dependent values
	diameter=35.0*numpy.sqrt(46.0/float(stations))
	ansource=int(nsource*(35.0/diameter)*(35.0/diameter))
	print "Station diameter = %s, number of sources = %s" % (diameter, ansource)
	rbase=hiono*wave/diameter
	rpierce=2.0*rbase
	print "Number of stations = %d, baseline radius = %s, pierce radius= %s" % (stations, rbase, rpierce)
	ts=TelSources()
	ts.construct(nsources=ansource, radius=wave/diameter)

	lowrandbool=TelArray()
	lowrandbool.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rbase, diameter=diameter)
	trialconfig=copy.deepcopy(lowrandbool)
	best=copy.deepcopy(lowrandbool)

# 	trialconfig=TelArray()
# 	trialconfig.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rbase, diameter=diameter)
# 
# 	best=TelArray()
# 	best.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rbase, diameter=diameter)
# 
	bestDistance=best.distance()
	T=0.1*bestDistance
	print "Initial distance metric = ", bestDistance
	tp=TelPiercings()
	tp.construct(ts,best,hiono=hiono,rmin=1.8)

	for trial in range(ntrials):
		trialconfig.stations=best.stations.copy()
		if fresh:
			trialconfig=TelArray()
			trialconfig.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rbase)
		else:
			trialconfig.shakehalo(bestDistance/2.0)
		distance=trialconfig.distance()
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
	srb=numpy.zeros(nnoll)
	for sourcetrial in range(nsourcetrials):
		ts.construct(nsources=nsource, radius=wave/diameter)
		print "Number of sources ", nsource, "Trial ", sourcetrial
		tp.construct(ts,best,hiono=300,rmin=1.8)
		srb=srb+(1.0/float(nsourcetrials))*tp.assess(nnoll=nnoll, rmax=rpierce, doplot=False)
	# Scale the x axis to be a pseudo wave number
	scales=numpy.sqrt(numpy.arange(nnoll)/(rpierce*rpierce))
	plt.plot(scales, numpy.sqrt(srb/slowbd[0]), color=color['%d'%stations])
	print "Zero sqrt(singular value) %s" % numpy.sqrt(srb[0])
	best.saveCSV('LOW_RANDOMBOOLARDY%d.csv' % stations)
	best.writeKML('LOW_RANDOMBOOLARDY%d.kml' % stations	)
	print " "
plt.title('Singular values 11, 21, 31, 41, 51, BDv1, %dm' % (wave))
plt.xlabel(r'$km^{-1}$')
plt.ylabel(r'$\sqrt{{\rm Singular value}}$')
plt.axes().set_xlim([0.0, 0.3])
plt.savefig('SVD_OPT_wave%sm.pdf' % wave)
