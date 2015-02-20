from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
random.seed(781490893)
ts=TelSources()
ts.construct(nsources=5, radius=6.0/35.0)

rshake=5.0

sbest=0
lowbest=TelArray()

lowrand=TelArray()
lowrand.construct('LOW_OPTIMUM', nstations=60, nhalo=60, rhalo=40)
lowfirst=lowrand

tp=TelPiercings()
tp.construct(ts,lowrand,rmin=1.0,hiono=300)
sbest = tp.assess(nnoll=50, rmax=50)
print "Initial shalf=%f" % sbest

for trial in range(100):
	random.seed(781490893)

	lowrand=lowfirst
	lowrand.shakehalo(rshake)

	tp=TelPiercings()
	tp.construct(ts,lowrand,rmin=1.0,hiono=300)
	shalf = tp.assess(nnoll=50, rmax=50)
	tp=False
	
	if shalf>sbest:
		print "Trial %d - New best array : shalf=%f, rshake=%f" % (trial, shalf, rshake)
		sbest=shalf
		lowbest=lowrand
		lowbest.plot()
	else:
		print "Trial %d - failed : shalf=%f" % (trial, shalf)
		