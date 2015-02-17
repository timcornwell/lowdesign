from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()
ts.construct(nsources=20, radius=6.0/35.0)

tp=TelPiercings()

lowrand=TelArray()
lowrand.circles('LOW_CIRCLES', nstations=1024, nhalo=184, rhalo=40.0)
lowrand.plot()
tp.construct(ts,lowrand,hiono=300)
tp.plot()
tp.assess(nnoll=50, rmax=50)
