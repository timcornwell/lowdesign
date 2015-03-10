from telopt import *


size=30.0
rmax=50.0
lowrand=TelArray()
lowrand.circles('LOW_CIRCLE_diameter=%d' % size, nstations=46, nhalo=46, rhalo=rmax)
lowrand.shakehalo(rshake=5.0)
lowrand.writeKML('LOW_CIRCLE_diameter=%d.kml' % size)

lowhand=TelArray()
lowhand.readKML('boolardyhandedit.kml')
lowhand.plot()

