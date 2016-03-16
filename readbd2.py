from telopt import *

	
lowboolardy=TelArray()
lowboolardy.readKML('LOW_BD2', 'v7ska1lowN1white.kml')
lowboolardy.plot(rmax=40.0)
lowboolardy.saveCSV('LOW_BD2.csv')
