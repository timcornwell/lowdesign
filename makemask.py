from telopt import *
from scipy import misc

bm=TelMask()
bm.readKML()
bm.plot()

bm.readMask(maskfile='Mask_BoolardyStation.png')

for i in range(0,80,10):
	print bm.masked(0.0, float(i))
