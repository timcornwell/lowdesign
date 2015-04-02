from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=6
nnoll=60

ntrials=100
rmax=40

lowrandbool=TelArray()
lowrandbool.randomBoolardy('LOW_RANDOMBOOLARDY46', nstations=46, nhalo=46, rhalo=40.0)

tp=TelPiercings()

random.seed(781490893)
nsource=6

srandbool=numpy.zeros(nnoll)
sbest=numpy.zeros(nnoll)

best=TelArray()
bestRmsS=0
for trial in range(ntrials):
	ts.construct(nsources=nsource+1, radius=3.0/35.0)

	lowrandbool=TelArray()
	lowrandbool.randomBoolardy('LOW_RANDOMBOOLARDY46', nstations=46, nhalo=46, rhalo=40.0)
	lowrandbool.mst(doplot=False)
	tp.construct(ts,lowrandbool,hiono=300,rmin=0)
	srandbool=tp.assess(nnoll=nnoll, rmax=rmax, doplot=(trial==0))
	rmsS=numpy.sqrt(numpy.average(srandbool*srandbool))
	if trial==0:
		plt.clf()
		plt.title('Singular values %d sources (r:Boolardy, g:Best)' % (nsource+1))
		plt.xlabel('Singular value index')
		plt.ylabel('Singular value')
	if (rmsS>bestRmsS):
		bestRmsS=rmsS
		best=lowrandbool
		sbest=srandbool
		print trial, " Found better config ", bestRmsS
		plt.plot(srandbool, color='r')

plt.plot(sbest, color='g')

lowbd=TelArray()
lowbd.readLOWBD('LOWBD')
lowbd.save('LOWBD_OPT.csv')
tp.construct(ts,lowbd,hiono=300,rmin=0)
slowbd=tp.assess(nnoll=nnoll, rmax=rmax, doplot=False)
plt.plot(slowbd, color='b')
	
plt.savefig('Singularvalue_%d_sources.pdf' % (nsource+1))

best.plot(rmax=40.0)
best.save('LOW_RANDOMBOOLARDY46_OPT.csv')
print best.mst()
