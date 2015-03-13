from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=10

c15circles46=numpy.zeros(nsources)
c15hand=numpy.zeros(nsources)
c15lowbd=numpy.zeros(nsources)
c15lofar=numpy.zeros(nsources)

ntrials=10
rmax=50

tuv=TelUV()

lowrand46=TelArray()
lowrand46.circles('LOW_CIRCLES46', nstations=1024, nhalo=46, rhalo=40.0)
lowrand46.shakehalo(rshake=5.0)
lowrand46.plot()
lowrand46.save('LOW_CIRCLES46_OPT.csv')
tuv.construct(lowrand46);tuv.plot()

lowhand=TelArray()
lowhand.readKML('boolardyhandedit.kml')
lowhand.plot()
lowhand.save('LOW_HAND_OPT.csv')
tuv.construct(lowhand);tuv.plot()


lowbd=TelArray()
lowbd.readLOWBD('LOWBD')
lowbd.plot()
lowbd.save('LOWBD_OPT.csv')
tuv.construct(lowbd);tuv.plot()

lofar=TelArray()
lofar.readLOFAR('LOFAR')
lofar.plot()
lofar.save('LOFAR_OPT.csv')
tuv.construct(lofar);tuv.plot()

tp=TelPiercings()

random.seed(781490893)
for trial in range(ntrials):
	print "Trial ", trial
	for nsource in range(nsources):
		ts.construct(nsources=nsource+1, radius=3.0/35.0)

		tp.construct(ts,lowrand46,hiono=300,rmin=0)
		if trial==0:
			tp.plot(rmax=rmax)
		c15circles46[nsource]=max(c15circles46[nsource],tp.assess(nnoll=60, rmax=rmax, doplot=(trial==0)))
	
		tp.construct(ts,lowhand,hiono=300,rmin=0)
		if trial==0:
			tp.plot(rmax=rmax)
		c15hand[nsource]=max(c15hand[nsource],tp.assess(nnoll=60, rmax=rmax, doplot=(trial==0)))
	
		tp.construct(ts,lowbd,hiono=300,rmin=1.8)
		if trial==0:
			tp.plot(rmax=rmax)
		c15lowbd[nsource]=max(c15lowbd[nsource],tp.assess(nnoll=60, rmax=rmax, doplot=(trial==0)))
	
		tp.construct(ts,lofar,hiono=300,rmin=1.8)
		if trial==0:
			tp.plot(rmax=rmax)
		c15lofar[nsource]=max(c15lofar[nsource],tp.assess(nnoll=60, rmax=rmax, doplot=(trial==0)))
	
plt.clf()
plt.plot(numpy.arange(nsources)+1.0, c15circles46, color='r')
plt.plot(numpy.arange(nsources)+1.0, c15hand, color='black')
plt.plot(numpy.arange(nsources)+1.0, c15lowbd, color='g')
plt.plot(numpy.arange(nsources)+1.0, c15lofar, color='b')
plt.title('Average singular value[0:15] (r:circles, g:BD, b:LOFAR, black:Boolardy)')
plt.xlabel('Number of sources')
plt.ylabel('Average singular value')
plt.savefig('averagesingularvalue.pdf')
