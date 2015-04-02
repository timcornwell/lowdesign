from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
#random.seed(781450893)
fresh=False

rhalo=15.0

#for stations in [11, 21, 31, 41, 46, 51]:
for stations in [46]:

	ntrials=1000000
	lowrandbool=TelArray()
	lowrandbool.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)

	trialconfig=TelArray()
	trialconfig.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)

	best=TelArray()
	best.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)

	bestDistance=best.distance()
	print "Initial distance metric = ", bestDistance
	frame=0
	for trial in range(ntrials):
		trialconfig.stations=best.stations.copy()
		if fresh:
			trialconfig=TelArray()
			trialconfig.randomBoolardy('LOW_RANDOMBOOLARDY%d' % stations, nstations=stations, nhalo=stations, rhalo=rhalo)
		else:
			trialconfig.shakehalo(0.1)
		distance=trialconfig.distance()
		if (distance>bestDistance):
			bestDistance=distance
			best.stations=trialconfig.stations.copy()
			print "Trial ", trial, " Found better config ", best.distance()
			best.plot(rmax=rhalo, plotfile='Array_RANDOMBOOLARDY%d_%d.jpg' % (stations, frame))
			frame=frame+1
		else:
			trialconfig.stations=best.stations.copy()
		
	best.plot(rmax=rhalo)
	best.save('LOW_RANDOMBOOLARDY%d.csv' % stations)
	print best.mst()
