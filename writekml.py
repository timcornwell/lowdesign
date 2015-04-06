from telopt import *

lowopt=TelArray()
lowopt.readLOWBD('LOW_RANDOMBOOLARDY46', l1def='LOW_RANDOMBOOLARDY46.csv')
lowopt.writeKML('LOW_RANDOMBOOLARDY46.kml')
