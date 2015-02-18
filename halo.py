from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

shalfcircles=numpy.zeros(30)
shalflowl1=numpy.zeros(30)
shalflofar=numpy.zeros(30)

for nsources in range(1,30):
	ts.construct(nsources=nsources, radius=6.0/35.0)
	tp=TelPiercings()
	lowrand=TelArray()
	lowrand.circles('LOW_CIRCLES', nstations=1024, nhalo=60, rhalo=40.0)
	lowrand.plot()
	tp.construct(ts,lowrand,hiono=300,rmin=0)
	tp.plot(rmax=50)
	shalfcircles[nsources]=tp.assess(nnoll=100, rmax=50)
	
	lowl1=TelArray()
	lowl1.readLOWL1('LOW_L1')
	lowl1.plot()
	tp.construct(ts,lowl1,hiono=300,rmin=2.5)
	tp.plot(rmax=50)
	shalflowl1[nsources]=tp.assess(nnoll=100, rmax=50)
	
	lofar=TelArray()
	lofar.readLOFAR('LOFAR')
	lofar.plot()
	tp.construct(ts,lofar,hiono=300,rmin=2.5)
	tp.plot(rmax=50)
	shalflofar[nsources]=tp.assess(nnoll=100, rmax=50)
	
plt.clf()
plt.plot(shalfcircles, color='r')
plt.plot(shalflowl1, color='g')
plt.plot(shalflofar, color='b')
plt.title('Shalf for (r:circles, g:Low_L1, and b:LOFAR) and varying number of sources ')
plt.xlabel('Number of sources')
plt.ylabel('Shalf')
plt.savefig('shalf.pdf')
