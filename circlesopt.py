from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

shalfcircles=numpy.zeros(10)
shalflowbd=numpy.zeros(10)
shalflofar=numpy.zeros(10)

ntrials=100

tuv=TelUV()

lowrand=TelArray()
lowrand.circles('LOW_CIRCLES', nstations=1024, nhalo=60, rhalo=40.0)
lowrand.shakehalo(rshake=5.0)
lowrand.plot()
tuv.construct(lowrand);tuv.plot()

lowbd=TelArray()
lowbd.readLOWBD('LOWBD')
lowbd.plot()
tuv.construct(lowbd);tuv.plot()

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
		shalfcircles[nsources]=shalfcircles[nsources]+(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=50, doplot=(trial==0))
	
		tp.construct(ts,lowbd,hiono=300,rmin=2.5)
		if trial==0:
			tp.plot(rmax=50)
		shalflowbd[nsources]=shalflowbd[nsources]+(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=50, doplot=(trial==0))
	
		tp.construct(ts,lofar,hiono=300,rmin=2.5)
		if trial==0:
			tp.plot(rmax=50)
		shalflofar[nsources]=shalflofar[nsources]+(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=50, doplot=(trial==0))
	
plt.clf()
plt.plot(shalfcircles, color='r')
plt.plot(shalflowbd, color='g')
plt.plot(shalflofar, color='b')
plt.title('Shalf (r:circles, g:BD, and b:LOFAR)')
plt.xlabel('Number of sources')
plt.ylabel('Shalf')
plt.savefig('shalf.pdf')
