from telopt import *

station=46
lowrand=TelArray()
lowrand.randomBoolardy('LOW_BOOLARDYRANDOM_%d'%station, nstations=station, nhalo=station, rhalo=40)
lowrand.plot(rmax=80)
print lowrand.mst()
