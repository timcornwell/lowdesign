from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()
ts.construct(nsources=20, radius=6.0/35.0)

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
tp.construct(ts,lowrand,rmin=2.5,hiono=300)
tp.plot()
tp.assess(nnoll=100, rmax=70)
tp=False

tplow=TelPiercings()
tplow.construct(ts,low,rmin=2.5,hiono=300)
tplow.plot()
tplow.assess(nnoll=100, rmax=70)
tplow=False

tplofar=TelPiercings()
tplofar.construct(ts,lofar,rmin=2.5,hiono=300)
tplofar.plot()
tplofar.assess(nnoll=100, rmax=70)
tplofar=False

