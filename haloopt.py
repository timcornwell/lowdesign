from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()
ts.construct(nsources=60, radius=6.0/35.0)

lowrand=TelArray()
lowrand.construct('LOW_RANDOM', nstations=1024, nhalo=45, rhalo=40)
lowrand.plot()

low=TelArray()
low.readLOWL1('LOW_L1')
low.plot()

lofar=TelArray()
lofar.readLOFAR('LOFAR')
lofar.plot()

tp=TelPiercings()
tp.construct(ts,lowrand,hiono=400)
tp.plot()
tp.assess(nnoll=50, rmax=70)
tp=False

tplow=TelPiercings()
tplow.construct(ts,low,hiono=400)
tplow.plot()
tplow.assess(nnoll=50, rmax=70)
tplow=False

tplofar=TelPiercings()
tplofar.construct(ts,lofar,hiono=400)
tplofar.plot()
tplofar.assess(nnoll=50, rmax=70)
tplofar=False

