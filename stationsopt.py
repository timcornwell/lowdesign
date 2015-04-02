from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=1

ntrials=1
rmax=50
nstationsmax=45

shalfrandom=numpy.zeros(nstationsmax)
stations=numpy.zeros(nstationsmax)

tuv=TelUV()

for station in range(5,nstationsmax):
	stations[station]=station
	lowrand=TelArray()
	lowrand.randomBoolardy('LOW_RANDOMBOOLARDY_%d'%station, nstations=station, nhalo=station, rhalo=40.0)
	lowrand.plot()
	print station, lowrand.mst()
	tuv.construct(lowrand);tuv.plot()

	tp=TelPiercings()

	random.seed(781490893)
	ts.construct(nsources=nsources, radius=3.0/35.0)
	for trial in range(ntrials):
		tp.construct(ts,lowrand,hiono=300,rmin=2.0)
		if trial==0:
			tp.plot(rmax=rmax)
		shalfrandom[station]=max(shalfrandom[station],tp.assess(nnoll=100, rmax=rmax, doplot=(trial==0)))
	
plt.clf()
plt.plot(stations, shalfrandom, color='r')
plt.title('S Average for fixed source density, random array')
plt.xlabel('Number of stations')
plt.ylabel('S ave')
plt.savefig('savestations.pdf')
