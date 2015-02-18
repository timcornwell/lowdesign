from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()
ts.construct(nsources=1, radius=6.0/35.0)

rshake=20.0

sbest=0
lowbest=TelArray()

lowrand=TelArray()
lowrand.construct('LOW_OPTIMUM', nstations=52, nhalo=52, rhalo=40)
lowbest=lowrand

tp=TelPiercings()
tp.construct(ts,lowrand,rmin=1.0,hiono=300)
sbest = tp.assess(nnoll=50, rmax=50)
print "Initial s=%f" % sbest

for trial in range(100):

	lowrand=lowbest
	lowrand.shakehalo(rshake)

	tp=TelPiercings()
	tp.construct(ts,lowrand,rmin=2.5,hiono=300)
	shalf = tp.assess(nnoll=50, rmax=50)
	tp=False
	
	if shalf>sbest:
		print "Trial %d - New best array : shalf=%f, rshake=%f" % (trial, shalf, rshake)
		sbest=shalf
		lowbest=lowrand
#		rshake=rshake*0.95
		lowbest.plot()
	else:
		print "Trial %d - failed : shalf=%f" % (trial, shalf)
		