from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=6
nnoll=200
diameter=35.0
wave=3.0
rbase=30.0
rpierce=rbase+hiono*wave/35.0

ntrials=1
rmax=300.0*wave/35.0

tuv=TelUV()
tp=TelPiercings()


random.seed(781490893)
for nsource in [nsources]:

	plt.clf()

	lowbd=TelArray()
	lowbd.readCSV('LOWBD', l1def='LOW_BD.csv')
	slowbd=numpy.zeros(nnoll)

	for trial in range(ntrials):
		ts.construct(nsources=nsource, radius=wave/diameter)

		print "Number of sources ", nsource, "Trial ", trial
		tp.construct(ts,lowbd,hiono=300,rmin=1.8)
		slowbd=slowbd+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=False)

# 	plt.semilogy(numpy.sqrt(slowbd/slowbd[0]), color='black')	
	plt.plot(numpy.sqrt(slowbd/slowbd[0]), color='black')	
	for file in ['LOW_RANDOMBOOLARDY11', 'LOW_RANDOMBOOLARDY21', 'LOW_RANDOMBOOLARDY31', 'LOW_RANDOMBOOLARDY41', 'LOW_RANDOMBOOLARDY51']:

		print file
		lowopt=TelArray()
		lowopt.readCSV(file, l1def='%s.csv' % file)

		sopt=numpy.zeros(nnoll)
		for trial in range(ntrials):
			ts.construct(nsources=nsource, radius=wave/diameter)

			print "Number of sources ", nsource, "Trial ", trial
			tp.construct(ts,lowopt,hiono=300,rmin=0)
			sopt=sopt+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rmax, doplot=False)
	
# 		plt.semilogy(numpy.sqrt(sopt/slowbd[0]))
		plt.plot(numpy.sqrt(sopt/slowbd[0]))
	plt.title('Singular values %d sources (black:BD, colors:11,21,31,41,51)' % (nsource))
	plt.xlabel('Singular value index')
	plt.ylabel(r'$\sqrt{{\rm Singular value}}$')
#  	plt.axes().set_ylim([1e2,1e5])
	plt.savefig('SVD_%d_sources.pdf' % nsource)
