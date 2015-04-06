from telopt import *

for file in ['LOW_RANDOMBOOLARDY11', 'LOW_RANDOMBOOLARDY31', 'LOW_RANDOMBOOLARDY46', 'LOW_RANDOMBOOLARDY21', 'LOW_RANDOMBOOLARDY41', 'LOW_RANDOMBOOLARDY51']:
	lowopt=TelArray()
	print file
	lowopt.readLOWBD('file', l1def='%s.csv' % file)
	lowopt.writeKML('%s.kml' % file)
