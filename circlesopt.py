from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

shalfcircles=numpy.zeros(10)
shalflowl1=numpy.zeros(10)
shalflofar=numpy.zeros(10)

ntrials=100

tuv=TelUV()

lowrand=TelArray()
lowrand.circles('LOW_CIRCLES', nstations=1024, nhalo=60, rhalo=40.0)
lowrand.shakehalo(rshake=5.0)
lowrand.plot()
tuv.construct(lowrand);tuv.plot()

lowl1=TelArray()
lowl1.readLOWL1('LOW_L1')
tuv.construct(lowl1);tuv.plot()

lofar=TelArray()
lofar.readLOFAR('LOFAR')
lofar.plot()
tuv.construct(lofar);tuv.plot()

tp=TelPiercings()

for nsources in range(1,10):
	ts.construct(nsources=nsources, radius=6.0/35.0)
	random.seed(781490893)
	for trial in range(ntrials):
		tp.construct(ts,lowrand,hiono=300,rmin=0)
		if trial==0:
			tp.plot(rmax=50)
		shalfcircles[nsources]=shalfcircles[nsources]+(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=50)
	
		tp.construct(ts,lowl1,hiono=300,rmin=2.5)
		if trial==0:
			tp.plot(rmax=50)
		shalflowl1[nsources]=shalflowl1[nsources]+(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=50)
	
		tp.construct(ts,lofar,hiono=300,rmin=2.5)
		if trial==0:
			tp.plot(rmax=50)
		shalflofar[nsources]=shalflofar[nsources]+(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=50)
	
plt.clf()
plt.plot(shalfcircles, color='r')
plt.plot(shalflowl1, color='g')
plt.plot(shalflofar, color='b')
plt.title('Shalf (r:circles, g:Low_L1, and b:LOFAR)')
plt.xlabel('Number of sources')
plt.ylabel('Shalf')
plt.savefig('shalf.pdf')
