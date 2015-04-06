from telopt import *
import math
import random

#
# Radius=1/2 FWZ ~ FWHM
#
#random.seed(781450893)
fresh=False

rhalo=40.0
T=0.1

#for stations in [11, 21, 31, 41, 46, 51]:
for stations in [11, 21]:

	ntrials=1000000
	lowrandbool=TelArray()
	lowrandbool.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)

	trialconfig=TelArray()
	trialconfig.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)

	best=TelArray()
	best.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)

	bestDistance=best.distance()
	print "Initial distance metric = ", bestDistance
	for trial in range(ntrials):
		trialconfig.stations=best.stations.copy()
		if fresh:
			trialconfig=TelArray()
			trialconfig.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)
		else:
			trialconfig.shakehalo(1.0)
		distance=trialconfig.distance()
		prob=math.exp(-(bestDistance-distance)/T)
		rnd=random.uniform(0.0,1.0)
		if (distance>bestDistance):
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
			print "Trial ", trial, " Found better config ", best.distance()
			T=T*0.99
		elif (rnd<prob):
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
			print "Trial ", trial, " Found worse config ", best.distance(), " but accepting anyway"
			print prob, rnd, T
		else:
			trialconfig.stations=best.stations.copy()
		
	best.plot(rmax=rhalo)
	best.save('LOW_RANDOMBOOLARDY%d.csv' % stations)
	best.writeKML('LOW_RANDOMBOOLARDY%d.kml' % stations	)
	print best.mst()
