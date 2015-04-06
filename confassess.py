from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=6
nnoll=60
diameter=35.0

ntrials=10
rmax=40

tuv=TelUV()
tp=TelPiercings()


random.seed(781490893)
for nsource in [nsources]:

	plt.clf()

	lowbd=TelArray()
	lowbd.readCSV('LOWBD', recenter=True)
	lowbd.recenter()
	slowbd=numpy.zeros(nnoll)
	
	for trial in range(ntrials):
		ts.construct(nsources=nsource, radius=3.0/diameter)

		print "Number of sources ", nsource, "Trial ", trial
		tp.construct(ts,lowbd,hiono=300,rmin=0.0)
		slowbd=slowbd+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=False)

	plt.semilogy(slowbd, color='black')	
		
	for file in ['LOW_RANDOMBOOLARDY11', 'LOW_RANDOMBOOLARDY21', 'LOW_RANDOMBOOLARDY31', 'LOW_RANDOMBOOLARDY41', 'LOW_RANDOMBOOLARDY51']:

		print file
		lowopt=TelArray()
		lowopt.readCSV(file, l1def='%s.csv' % file)

		sopt=numpy.zeros(nnoll)
		for trial in range(ntrials):
			ts.construct(nsources=nsource, radius=3.0/diameter)

			print "Number of sources ", nsource+1, "Trial ", trial
			tp.construct(ts,lowopt,hiono=300,rmin=0)
			sopt=sopt+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=False)
	
		plt.semilogy(sopt)
	plt.title('Singular values %d sources (black:BD, colors:11,21,31,41,51)' % (nsource+1))
	plt.xlabel('Singular value index')
	plt.ylabel('Singular value')
 	plt.axes().set_ylim([1e2,1e5])
	plt.savefig('SVD_%d_sources.pdf' % (nsource+1))
