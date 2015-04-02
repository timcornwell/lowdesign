from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=10
nnoll=60
diameter=70.0

c15opt=numpy.zeros(nsources)
c15lowbd=numpy.zeros(nsources)

ntrials=1
rmax=15

tuv=TelUV()

lowopt=TelArray()
lowopt.readLOWBD('LOW_RANDOMBOOLARDY46', l1def='LOW_RANDOMBOOLARDY46.csv')
lowopt.plot(rmax=rmax)
print lowopt.mst()
tuv.construct(lowopt);tuv.plot()

lowbd=TelArray()
lowbd.readLOWBD('LOWBD')
lowbd.plot(rmax=rmax)
print lowbd.mst()
tuv.construct(lowbd);tuv.plot()

tp=TelPiercings()

random.seed(781490893)
for nsource in range(nsources):

	sopt=numpy.zeros(nnoll)
	slowbd=numpy.zeros(nnoll)
	
	for trial in range(ntrials):
		ts.construct(nsources=nsource+1, radius=3.0/diameter)

		print "Number of sources ", nsource+1, "Trial ", trial
		tp.construct(ts,lowopt,hiono=300,rmin=0)
		if trial==0:
			tp.plot(rmax=rmax)
		sopt=sopt+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))
	
		tp.construct(ts,lowbd,hiono=300,rmin=0.0)
		if trial==0:
			tp.plot(rmax=rmax)
		slowbd=slowbd+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))

	print "sopt[0] %f slowbd[0] %f" % (sopt[0], slowbd[0])	
	plt.clf()
	plt.semilogy(slowbd, color='r')
	plt.semilogy(sopt, color='g')
	plt.title('Singular values %d sources (r:BD, g:opt)' % (nsource+1))
	plt.xlabel('Singular value index')
	plt.ylabel('Singular value')
	plt.axes().set_ylim([1e2,1e5])
	plt.savefig('Singularvalue_%d_sources.pdf' % (nsource+1))
