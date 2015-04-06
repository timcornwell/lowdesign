import numpy
import math
import scipy
import matplotlib
import random
import matplotlib.pyplot as plt
import csv
import zernike
from scipy import linalg
import Image
import numpy as np
import scipy.spatial.distance as sd 
from mst import * 

#random.seed(781490893)

class TelUtils:

	def uniformcircle(self, n, rhalo=1.0):
		x=numpy.zeros(n)
		y=numpy.zeros(n)
		for i in range(n):
			phi=2*numpy.pi*random.random()
			r=rhalo*numpy.sqrt(random.random())
			x[i]=r*numpy.cos(phi)
			y[i]=r*numpy.sin(phi)
		return x, y
		
class TelMask:
	def _init_(self):
		self.name=''
		self.construct()
		
	def readMask(self, maskfile='Mask_BoolardyStation.png'):
		self.mask = scipy.misc.imread('Mask_BoolardyStation.png')
		self.center={}
		self.center['x']=self.mask.shape[0]/2
		self.center['y']=self.mask.shape[1]/2
		self.scale={}
		self.scale['x']=2*80.0/self.mask.shape[0]
		self.scale['y']=2*80.0/self.mask.shape[1]
		
	def masked(self, x, y):
		mx=+int(-y/self.scale['x']+self.center['x'])
		my=+int(+x/self.scale['y']+self.center['y'])
		if  self.mask[mx,my,0] == 255:
			return True
		else:
			return False
		

	def readKML(self, name='BoolardyStation', kmlfile="BoolardyStation2(approx).kml"):
	
		long0=116.779167
		lat0=-26.789267
		Re=6371.0
		nsegments=55
		self.segments={}
		self.segments['x1']=numpy.zeros(nsegments)
		self.segments['y1']=numpy.zeros(nsegments)
		self.segments['x2']=numpy.zeros(nsegments)
		self.segments['y2']=numpy.zeros(nsegments)
		self.name=name
		f=open(kmlfile)
		segment=0
		nextline=False
		for line in f:
			line=line.lstrip()
			if nextline:
				part=line.split(' ')[0].split(',')
				x=float(part[0])
				y=float(part[1])
				self.segments['x1'][segment]=(x-long0)*Re*numpy.pi/(180.0*numpy.cos(numpy.pi*lat0/180.0))
				self.segments['y1'][segment]=(y-lat0)*Re*numpy.pi/(180.0*numpy.cos(numpy.pi*lat0/180.0))
				part=line.split(' ')[1].split(',')
				x=float(part[0])
				y=float(part[1])
				self.segments['x2'][segment]=(x-long0)*Re*numpy.pi/(180.0*numpy.cos(numpy.pi*lat0/180.0))
				self.segments['y2'][segment]=(y-lat0)*Re*numpy.pi/(180.0*numpy.cos(numpy.pi*lat0/180.0))
				nextline=False
				segment=segment+1
			if line.find('</coordinates>') > -1:
				nextline=False
			elif line.find('<coordinates>') > -1:
				nextline=True
				
	def plot(self, rmax=20):
		plt.clf()
		plt.fill(self.segments['x1'], self.segments['y1'],fill=False,color='blue')
		plt.axes().set_xlim([-80.0,80.0])
		plt.axes().set_ylim([-80.0,80.0])
		plt.axes().set_aspect('equal')
		plt.savefig('Mask_%s_frame.png' % self.name)
		plt.fill(self.segments['x1'], self.segments['y1'],fill=True,color='blue')
		plt.axes().get_xaxis().set_visible(False)
		plt.axes().get_yaxis().set_visible(False)
# 		plt.axes().set_frame_on(False)
		plt.savefig('Mask_%s.png' % self.name)

class TelArray:

	def _init_(self):
		self.name=''
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		self.construct()
	
	def recenter(self):
		self.center={}
		self.center['x']=numpy.average(self.stations['x'])
		self.center['y']=numpy.average(self.stations['y'])
		self.stations['x']=self.stations['x']-self.center['x']
		self.stations['y']=self.stations['y']-self.center['y']

	def plot(self, rmax=40, plotfile=''):
		plt.clf()
		plt.title('Antenna locations %s' % self.name)
		plt.xlabel('X (km)')
		plt.ylabel('Y (km)')
		plt.plot(self.stations['x'], self.stations['y'], '.')
		plt.axes().set_aspect('equal')
		circ=plt.Circle((0,0), radius=rmax, color='g', fill=False)
		fig = plt.gcf()
		fig.gca().add_artist(circ)
		maxaxis=1.1*max(numpy.max(abs(self.stations['x'])), numpy.max(abs(self.stations['y'])))
		plt.axes().set_xlim([-maxaxis,maxaxis])
		plt.axes().set_ylim([-maxaxis,maxaxis])
		mask=TelMask()
		mask.readKML()
		plt.fill(mask.segments['x1'], mask.segments['y1'], fill=False)
		if plotfile=='':
			plotfile='Array_%s.pdf' % self.name
		plt.savefig(plotfile)
	
	def random(self, name='Stations', rhalo=40, rcore=1.0, nstations=512, nhalo=45, nantennas=256, fobs=1e8, diameter=35.0):
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		self.name=name
		self.rhalo=rhalo
		self.rcore=rcore
		self.nstations=nstations
		self.center={}
		self.center['x']=0.0
		self.center['y']=0.0
		ncore=self.nstations-nhalo
		self.fobs=fobs
		self.diameter=diameter
		self.stations={}
		self.stations['x'], self.stations['y']=TelUtils().uniformcircle(self.nstations, self.rhalo)
		self.stations['x'][:ncore], self.stations['y'][:ncore]=TelUtils().uniformcircle(ncore, self.rcore)
		self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)

	def randomBoolardy(self, name='RandomBoolardy', rhalo=40.0, rcore=1.0, nstations=512, nhalo=45, nantennas=256, fobs=1e8, diameter=35.0):
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		self.name=name
		self.rhalo=rhalo
		self.rcore=rcore
		self.nstations=nstations
		self.center={}
		self.center['x']=0.0
		self.center['y']=0.0
		ncore=self.nstations-nhalo
		self.fobs=fobs
		self.diameter=diameter
		self.stations={}
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
		self.stations['weight']=numpy.zeros(self.nstations)
		self.stations['x'][:ncore], self.stations['y'][:ncore]=TelUtils().uniformcircle(ncore, self.rcore)
		inhalo=ncore
		while inhalo < nstations:
			x, y=TelUtils().uniformcircle(1, self.rhalo)
			if not self.mask.masked(x,y):
				r=numpy.sqrt(x*x+y*y)
				if r<rhalo:
					self.stations['x'][inhalo]=x
					self.stations['y'][inhalo]=y
					inhalo=inhalo+1
		self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)

	def circles(self, name='Stations', rhalo=40, rcore=1.0, nstations=512, nhalo=44, fobs=1e8, diameter=35.0):
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		self.name=name
		self.rhalo=rhalo
		self.rcore=rcore
		self.nstations=nstations
		self.center={}
		self.center['x']=0.0
		self.center['y']=0.0
		ncore=self.nstations-nhalo
		self.fobs=fobs
		self.diameter=diameter
		if nhalo==60:
			self.nrings=4
			self.r=[0.0, rhalo/3.0, 2*rhalo/3.0, rhalo]
			self.nonring=[1, 9, 21, 29]
		elif nhalo==46:
			self.nrings=4
			self.r=[0.0, rhalo/3.0, 2*rhalo/3.0, rhalo]
			self.nonring=[1, 9, 13, 23]
		elif nhalo==30:
			self.nrings=4
			self.r=[0.0, rhalo/3.0, 2*rhalo/3.0, rhalo]
			self.nonring=[1, 5, 9, 15]
		elif nhalo==12:
			self.nrings=3
			self.r=[0.0, rhalo/2.0, rhalo]
			self.nonring=[1, 5, 6]
		else:
			nhalo=185
			self.nrings=7
			self.r=[0.0, rhalo/6.0, 2*rhalo/6.0, 3.0*rhalo/6.0, 4.0*rhalo/6.0, 5.0*rhalo/6.0, rhalo]
			self.nonring=[1, 9, 19, 27, 35, 43, 51]
		self.nstations=nhalo
		self.stations={}
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
		self.stations['weight']=numpy.zeros(self.nstations)
# 		self.stations['x'][:ncore], self.stations['y'][:ncore]=TelUtils().uniformcircle(ncore, self.rcore)
		station=0
		for ring in range(self.nrings):
			dphi=2*numpy.pi/self.nonring[ring]
			phi=0.0
			for spoke in range(self.nonring[ring]):
				self.stations['x'][station]=self.r[ring]*numpy.cos(phi)
				self.stations['y'][station]=self.r[ring]*numpy.sin(phi)
				self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)
				phi=phi+dphi
				station=station+1

	def shakehalo(self, rshake=5.0):
		newstations={}
		newstations['x']=self.stations['x'].copy()
		newstations['y']=self.stations['y'].copy()
		for station in range(1,self.nstations):
			cr=numpy.sqrt(self.stations['x'][station]*self.stations['x'][station]+self.stations['y'][station]*self.stations['y'][station])
			if cr>self.rcore:
				phi=2.0*numpy.pi*random.random()
				r=rshake*numpy.sqrt(random.random())
				x=newstations['x'][station]+r*numpy.cos(phi)
				y=newstations['y'][station]+r*numpy.sin(phi)
				rdr=numpy.sqrt(x*x+y*y)
 				if rdr<self.rhalo:
  					if not self.mask.masked(x,y):
  						newstations['x'][station]=x
 						newstations['y'][station]=y
 			self.stations=newstations

	def readCSV(self, name='LOWBD', rcore=0.0, l1def='SKA-low_config_baseline_design_arm_stations_2013apr30.csv', recenter=False):
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		self.name=name
		self.nstations=0
		self.stations={}
		self.rhalo=80
		self.fobs=1e8
		self.diameter=35.0
		meanx=0
		meany=0
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				meanx=meanx+float(row[1])
				meany=meany+float(row[0])
				self.nstations=self.nstations+1
		meanx=meanx/self.nstations
		meany=meany/self.nstations
		if not recenter:
			meanx=0.0
			meany=0.0
		f.close()
		self.nstations=0
		scale=0.001
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				x=scale*(float(row[0])-meanx)
				y=scale*(float(row[1])-meany)
				r=numpy.sqrt(x*x+y*y)
				self.nstations=self.nstations+1
		
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
		station=self.nstations-1
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				x=scale*(float(row[0])-meanx)
				y=scale*(float(row[1])-meany)
				r=numpy.sqrt(x*x+y*y)
				self.stations['x'][station]=x
				self.stations['y'][station]=y
				self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)
				station=station-1

	def readLOWBD(self, name='LOWBD', rcore=0.0, l1def='SKA-low_config_baseline_design_arm_stations_2013apr30.csv'):
		return self.readCSV(name, rcore, l1def)

	def readLOWL1(self, name='LOWL1', rcore=0.0, l1def='L1_configuration.csv'):
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		self.name=name
		self.nstations=0
		self.stations={}
		self.rhalo=80
		self.fobs=1e8
		self.diameter=35.0
		meanx=0
		meany=0
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				meanx=meanx+float(row[1])
				meany=meany+float(row[0])
				self.nstations=self.nstations+1
		meanx=meanx/self.nstations
		meany=meany/self.nstations
		f.close()
		self.nstations=0
		scale=6.371e6*numpy.pi/180000.0
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				x=scale*(float(row[1])-meanx)*numpy.cos(meanx*numpy.pi/180.0)
				y=scale*(float(row[0])-meany)
				r=numpy.sqrt(x*x+y*y)
				if r>rcore:
					self.nstations=self.nstations+1
		
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
		station=0
		scale=6.371e6*numpy.pi/180000.0
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				x=scale*(float(row[1])-meanx)*numpy.cos(meanx*numpy.pi/180.0)
				y=scale*(float(row[0])-meany)
				r=numpy.sqrt(x*x+y*y)
				if r>rcore:
					self.stations['x'][station]=x
					self.stations['y'][station]=y
					self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)
					station=station+1
		self.recenter()
		
	def readLOFAR(self, name='LOFAR', stationtype='S', band='HBA', lfdef='LOFAR.csv', lat=52.7):
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		cs=numpy.cos(numpy.pi*lat/180.0)
		sn=numpy.sin(numpy.pi*lat/180.0)
		self.name=name
		self.nstations=0
		self.stations={}
		self.rhalo=80
		self.fobs=1e8
		self.diameter=35.0
		meanx=0
		meany=0
		meanz=0
		with open(lfdef, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				type=row[0]
				if (type.find(stationtype) > -1) and (type.find(band) > -1):
					meanx=meanx+float(row[2])
					meany=meany+float(row[1])
					meanz=meanz+float(row[3])
					self.nstations=self.nstations+1
		meanx=meanx/self.nstations
		meany=meany/self.nstations
		meanz=meanz/self.nstations
		f.close()
		station=0
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
		self.stations['weight']=numpy.zeros(self.nstations)
		with open(lfdef, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				type=row[0]
				if (type.find(stationtype) > -1) and (type.find(band) > -1):
					x=(float(row[2])-meanx)/1000.0
					y=(float(row[1])-meany)/1000.0
					z=(float(row[3])-meanz)/1000.0
					self.stations['x'][station]=x
					self.stations['y'][station]=-cs*y+sn*z
					self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)
					station=station+1
		self.recenter()

	def save(self, filename='LOWBD.csv'):
		with open(filename, 'wb') as fp:
			rowwriter = csv.writer(fp)
			for station in range(self.nstations):
				rowwriter.writerow([1000.0*self.stations['y'][station],-1000.0*self.stations['x'][station]])
				
	def readKML(self, name='KML', kmlfile="Boolardy.kml", diameter=35.0):
	
		self.mask=TelMask()
		self.mask.readMask(maskfile='Mask_BoolardyStation.png')
		long0=116.779167
		lat0=-26.789267
		Re=6371.0
		self.stations={}
		self.stations['x']=numpy.zeros(46)
		self.stations['y']=numpy.zeros(46)
		self.name=name
		self.diameter=diameter
		f=open(kmlfile)
		self.nstations=46
		for line in f:
			line=line.lstrip()
			if line.find("name")>0:
				if line.find("Station")>0:
					station=int(line.split('Station')[1].split('<')[0])
			if line.find("coordinates")>0:
				x= float(line.split('>')[1].split('<')[0].split(',')[0])
				y= float(line.split('>')[1].split('<')[0].split(',')[1])
				self.stations['x'][station]=(x-long0)*Re*numpy.pi/(180.0*numpy.cos(numpy.pi*lat0/180.0))
				self.stations['y'][station]=(y-lat0)*Re*numpy.pi/(180.0*numpy.cos(numpy.pi*lat0/180.0))
				self.stations['weight']=self.diameter*self.diameter*self.diameter*self.diameter*float(self.nstations)
		self.center['x']=0.0
		self.center['y']=0.0
				
	def writeKML(self, kmlfile="LOW_CIRCLES.kml"):
	
		long0=116.779167
		lat0=-26.789267
		Re=6371.0
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
			'<name>Station %d</name>', \
			'<Point>', \
			'<coordinates>%f, %f</coordinates>', \
			'</Point>', \
			'</Placemark>']
		e=['</Document>', '</kml>']
		f=open(kmlfile, 'w')
		for ss in s:
			f.write(ss)
		for station in range(self.nstations):
			long= long0-180.0*(self.stations['y'][station])*numpy.cos(numpy.pi*lat0/180.0)/(Re*numpy.pi)
			lat = lat0 +180.0*(self.stations['x'][station])*numpy.cos(numpy.pi*lat0/180.0)/(Re*numpy.pi)
			f.write( l[0])
			f.write( l[1])
			f.write( l[2] % station)
			f.write( l[3])
			f.write( l[4] % (long, lat))
			f.write( l[5])
			f.write( l[6])
		f.write( e[0])
		f.write( e[1])
		
	def distance(self):
		P=numpy.zeros([self.nstations,2])
		P[...,0]=self.stations['x']
		P[...,1]=self.stations['y']
		distancemat=sd.pdist(P)
		distance=numpy.min(distancemat)
		return distance
	
	def mst(self, doplot=True):
		P=numpy.zeros([self.nstations,2])
		P[...,0]=self.stations['x']
		P[...,1]=self.stations['y']
		X=sd.squareform(sd.pdist(P))
		edge_list = minimum_spanning_tree(X)
		if doplot:
			plt.clf()
			plt.scatter(P[:, 0], P[:, 1]) 
			dist=0    
			for edge in edge_list:
				i, j = edge
				plt.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], c='r')
				dist=dist+numpy.sqrt((P[i,0]-P[j,0])*(P[i,0]-P[j,0])+(P[i,1]-P[j,1])*(P[i,1]-P[j,1]))
			plt.title('MST for %s, distance=%.1f km' % (self.name, dist))
			plt.xlabel('X (km)')
			plt.ylabel('y (km)')
			plt.axes().set_aspect('equal')
			maxaxis=numpy.max(abs(P))
			plt.axes().set_xlim([-maxaxis,maxaxis])
			plt.axes().set_ylim([-maxaxis,maxaxis])
			plt.savefig('%s_MST.pdf' %self.name)
		else:
			dist=0    
			for edge in edge_list:
				i, j = edge
				dist=dist+numpy.sqrt((P[i,0]-P[j,0])*(P[i,0]-P[j,0])+(P[i,1]-P[j,1])*(P[i,1]-P[j,1]))
			return dist

class TelUV:
	def _init_(self):
		self.name=''
		self.uv=False
	
	def construct(self, t):
		self.name=t.name
		self.nbaselines=t.nstations*t.nstations
		self.uv={}
		self.uv['x']=numpy.zeros(self.nbaselines)
		self.uv['y']=numpy.zeros(self.nbaselines)
		for station in range(t.nstations):
			self.uv['x'][station*t.nstations:(station+1)*t.nstations]=t.stations['x']-t.stations['x'][station]
			self.uv['y'][station*t.nstations:(station+1)*t.nstations]=t.stations['y']-t.stations['y'][station]
		
	
	def plot(self):
		self.plotter=True
		plt.clf()
		plt.title('UV Sampling %s' % self.name)
		plt.xlabel('X (km)')
		plt.ylabel('Y (km)')
		plt.plot(self.uv['x'], self.uv['y'], '.')
		plt.axes().set_aspect('equal')
		plt.savefig('UVcoverage_%s.pdf' % self.name)
		
	def assess(self):
		return 1.0
	
#
# Sources on the celestial sphere.
#
class TelSources:
	def _init_(self):
		self.name='Sources'
		self.nsources=100
	
	def construct(self, name='Sources', nsources=100, radius=1):
		self.name=name
		self.sources={}
		self.sources['x'], self.sources['y']=TelUtils().uniformcircle(nsources, radius)

		self.sources['x']=self.sources['x']-numpy.sum(self.sources['x'])/float(nsources)
		self.sources['y']=self.sources['y']-numpy.sum(self.sources['y'])/float(nsources)
		self.nsources=nsources
		self.radius=radius

	def plot(self):
		self.plotter=True
	
	def assess(self):
		return 1.0
#
# Piercings through the ionosphere
#
class TelPiercings:
	def _init_(self):
		self.name='Piercings'
		self.npiercings=0
		self.hiono=400
	
	def plot(self, rmax=70):
		plt.clf()
		plt.title('Piercings %s' % self.name)
		plt.xlabel('X (km)')
		plt.ylabel('Y (km)')
		plt.plot(self.piercings['x'], self.piercings['y'], '.')
		plt.axes().set_aspect('equal')
		circ=plt.Circle((0,0), radius=rmax, color='g', fill=False)
		fig = plt.gcf()
		fig.gca().add_artist(circ)
		maxaxis=max(numpy.max(abs(self.piercings['x'])), numpy.max(abs(self.piercings['y'])))
		plt.axes().set_xlim([-maxaxis,maxaxis])
		plt.axes().set_ylim([-maxaxis,maxaxis])
		plt.savefig('%s.pdf' % self.name)
	
	def construct(self, sources, array, rmin=1, hiono=300):
		self.hiono=hiono
		r2=array.stations['x']*array.stations['x']+array.stations['y']*array.stations['y']
		outside={}
		outside['x']=array.stations['x'][r2>=rmin*rmin]
		outside['y']=array.stations['y'][r2>=rmin*rmin]
		nstations=len(outside['x'])
		self.npiercings=sources.nsources*nstations
		self.name='Piercings_%d_%s_%s' % (sources.nsources, sources.name, array.name)
		self.piercings={}
		self.piercings['x']=numpy.zeros(self.npiercings)
		self.piercings['y']=numpy.zeros(self.npiercings)
		self.piercings['weight']=numpy.zeros(self.npiercings)
		for source in range(sources.nsources):
			self.piercings['x'][source*nstations:(source+1)*nstations]=self.hiono*sources.sources['x'][source]+outside['x']
			self.piercings['y'][source*nstations:(source+1)*nstations]=self.hiono*sources.sources['y'][source]+outside['y']
			self.piercings['weight'][source*nstations:(source+1)*nstations]=array.stations['weight']

	def assess(self, rmax=70.0, nnoll=20, doplot=True):
		A=numpy.zeros([self.npiercings, nnoll])
		for piercing in range(self.npiercings):
			x=self.piercings['x'][piercing]
			y=self.piercings['y'][piercing]
			weight=numpy.sqrt(self.piercings['weight'][piercing])
			r=numpy.sqrt(x*x+y*y)
			phi=numpy.arctan2(y,x)
			if(r<rmax):
				for noll in range(nnoll):
					A[piercing,noll]=weight*zernike.zernikel(noll,r/rmax,phi)
		Covar_A=numpy.zeros([nnoll, nnoll])
		for nnol1 in range(nnoll):
			for nnol2 in range(nnoll):
				Covar_A[nnol1,nnol2]=numpy.sum(A[...,nnol1]*A[...,nnol2])
		U,s,Vh = linalg.svd(Covar_A)
		s=numpy.sqrt(s)
		if doplot:
			plt.clf()
			plt.title('%s rmax=%d' % (self.name, rmax))
			plt.xlabel('Singular vector index')
			plt.ylabel('Sqrt(Singular value)')
			plt.plot(s)
			plt.savefig('%s_rmax=%d_SVD.pdf' % (self.name, rmax))
		return s
	