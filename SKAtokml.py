
import csv


def convert(csvfile="SKA1-Mid-May2015.csv", kmlfile="SKA1-Mid-May2015.kml"):
	s=['<?xml version="1.0" encoding="UTF-8"?>', \
		'<kml xmlns="http://www.opengis.net/kml/2.2">', \
		'<Document>', \
		'<Style id="whitecirc">', \
		'<IconStyle>', \
		'<Icon>', \
		'<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>', \
		'</Icon>', \
		'</IconStyle>', \
		'</Style>', \
		'<!--name></name-->']
	l=['<Placemark>', \
		'<styleUrl>#whitecirc</styleUrl>', \
		'<name>S%d</name>', \
		'<Point>', \
		'<coordinates>%f, %f</coordinates>', \
		'</Point>', \
		'</Placemark>']
	e=['</Document>', '</kml>']
	fout=open(kmlfile, 'w')
	for ss in s:
		fout.write(ss)
	with open(csvfile, 'rU') as fin:
		reader = csv.reader(f)
		for row in reader:
			station=row[0]
			long=row[1]
			lat =row[2]
			fout.write( l[0])
			fout.write( l[1])
			fout.write( l[2] % station)
			fout.write( l[3])
			fout.write( l[4] % (long, lat))
			fout.write( l[5])
			fout.write( l[6])
	fout.write( e[0])
	fout.write( e[1])
		
convert(csvfile="SKA1-Mid-May2015.csv", kmlfile="SKA1-Mid-May2015.kml")
		
