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
wave=6.0

# Diameter of station (m)
diameter=35.0

# Baseline length (radius)
rbase=40.0

# Radius of pierce points in the ionosphere
rpierce=rbase+hiono*wave/(diameter)

# Radius of core pierce points in the ionosphere
rpiercecore=hiono*wave/(diameter)

# Number of source trials
nsourcetrials=1

# Number of array trials
ntrials=1000000

# Number of stations
nstations = 25

# Number of core stations
ncore=0
##########################################################################################
# Recalculate all dependent values
print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
print "Number of stations = %d, baseline radius = %.3f (km), pierce radius= %.3f (km)" % (nstations, rbase, rpierce)

sname='LOW_RB%d_source%d_wave%.1f_diam%.1f_base%.1f' % (nstations, nsource, wave, diameter, rbase)
lowrandbool=TelArray()
lowrandbool.randomBoolardy(sname, nstations=nstations, rcore=3.0, nhalo=nstations, rhalo=rbase, diameter=diameter)
trialconfig=copy.deepcopy(lowrandbool)
best=copy.deepcopy(lowrandbool)

bestDistance=best.excessDistance(False)
T=2.0*bestDistance
print "Initial distance metric = %.3f (km), T=%.3f (km)" % (bestDistance, T)
trial=1
while (trial < ntrials) and (T > 0.001):
	trialconfig.stations=best.stations.copy()
	if bestDistance<0.001:
		print "Stalled - resetting"
		trialconfig=TelArray()
		trialconfig.randomBoolardy(sname, nstations=nstations, nhalo=nstations, rhalo=rbase)
	else:
		trialconfig.shakehalo(best.distance()/4.0)
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
	
best.shakehalo(best.distance()/2.0, one=False)

print "Distance %.3f (km) MST %.3f (km)" % (best.distance(), best.mst(False))
best.saveCSV('%s.csv' % sname)
best.writeWGS84('%s_WGS84.csv' % sname)
best.writeKML('%s.kml' % sname)
best.plot(rbase, '%s_ARRAY.pdf' % sname)
best.mst(True, '%s_MST.pdf' % sname)
uv=TelUV()
uv.construct(best)
uv.plot('%s_UV.pdf' % sname)
print " "	
	
print "Wavelength %.2f (m), Station diameter = %.1f (m), number of sources = %d" % (wave, diameter, nsource)
print "Number of stations = %d, baseline radius = %.3f (km), pierce radius= %.3f (km)" % (nstations, rbase, rpierce)

ts=TelSources()
tp=TelPiercings()
srb=numpy.zeros(nnoll)
for sourcetrial in range(nsourcetrials):
	ts.construct(nsources=nsource, radius=wave/diameter)
	print "Number of sources ", nsource, "Trial ", sourcetrial
	tp.construct(ts,best,hiono=hiono,rmin=0.0)
	srb=srb+tp.assess(nnoll=nnoll, rmax=rpierce, doplot=True)
	if sourcetrial == 0:
		tp.plot(rmax=rpierce, rcore=rpiercecore)
srb=srb/(float(nsourcetrials))
print "RB Singular value[0] %.3f" % srb[0]

