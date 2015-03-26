from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=10
nnoll=60

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
lowrand46.plot(rmax=40.0)
lowrand46.save('LOW_CIRCLES46_OPT.csv')
tuv.construct(lowrand46);tuv.plot()

lowhand=TelArray()
lowhand.readKML('Boolardy', 'boolardyhandedit.kml')
lowhand.plot(rmax=40.0)
lowhand.save('LOW_HAND_OPT.csv')
tuv.construct(lowhand);tuv.plot()

lowbd=TelArray()
lowbd.readLOWBD('LOWBD')
lowbd.plot(rmax=40.0)
lowbd.save('LOWBD_OPT.csv')
tuv.construct(lowbd);tuv.plot()

lofar=TelArray()
lofar.readLOFAR('LOFAR')
lofar.plot(rmax=40.0)
lofar.save('LOFAR_OPT.csv')
tuv.construct(lofar);tuv.plot()

tp=TelPiercings()

random.seed(781490893)
for nsource in range(nsources):

	scircles46=numpy.zeros(nnoll)
	shand=numpy.zeros(nnoll)
	slowbd=numpy.zeros(nnoll)
	slofar=numpy.zeros(nnoll)
	
	for trial in range(ntrials):
		ts.construct(nsources=nsource+1, radius=3.0/35.0)

		print "Number of sources ", nsource+1, "Trial ", trial
		tp.construct(ts,lowrand46,hiono=300,rmin=0)
		if trial==0:
			tp.plot(rmax=rmax)
		scircles46=scircles46+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))
	
		tp.construct(ts,lowhand,hiono=300,rmin=0)
		if trial==0:
			tp.plot(rmax=rmax)
		shand=shand+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))
	
		tp.construct(ts,lowbd,hiono=300,rmin=1.8)
		if trial==0:
			tp.plot(rmax=rmax)
		slowbd=slowbd+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))
	
		tp.construct(ts,lofar,hiono=300,rmin=1.8)
		if trial==0:
			tp.plot(rmax=rmax)
		slofar=slofar+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))
	
	plt.clf()
	plt.semilogy(scircles46, color='r')
	plt.semilogy(shand, color='black')
	plt.semilogy(slowbd, color='g')
	plt.semilogy(slofar, color='b')
	plt.title('Singular values %d sources (r:circles, g:BD, b:LOFAR, black:Boolardy)' % (nsource+1))
	plt.xlabel('Singular value index')
	plt.ylabel('Singular value')
	plt.axes().set_ylim([1e-3,1e3])
	plt.savefig('Singularvalue_%d_sources.pdf' % (nsource+1))
