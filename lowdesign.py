import numpy
import scipy
import matplotlib
import random
import matplotlib.pyplot as plt
import csv

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

class TelArray:

	def _init_(self):
		self.name=''
		self.construct()
	
	def plot(self):
		plt.clf()
		plt.title('Antenna locations %s' % self.name)
		plt.xlabel('X (km)')
		plt.ylabel('Y (km)')
		plt.plot(self.stations['x'], self.stations['y'], '.')
		plt.axes().set_aspect('equal')
		plt.savefig('Array_%s.pdf' % self.name)
	
	def construct(self, name='Stations', rhalo=40, rcore=2.5, nstations=512, nhalo=45, nantennas=256, fobs=1e8, diameter=0.35):
		self.name=name
		self.rhalo=rhalo
		self.rcore=rcore
		self.nstations=nstations
		ncore=self.nstations-nhalo
		self.nantennas=nantennas
		self.fobs=fobs
		self.diameter=diameter
		self.stations={}
		self.stations['x'], self.stations['y']=TelUtils().uniformcircle(self.nstations, self.rhalo)
		self.stations['x'][:ncore], self.stations['y'][:ncore]=TelUtils().uniformcircle(ncore, self.rcore)

	def readLOWL1(self, name='LOWL1', rcore=0.0, l1def='L1_configuration.csv'):
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
		station=0
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
					station=station+1
		
	def readLOFAR(self, name='LOFAR', stationtype='S', band='HBA', lfdef='LOFAR.csv'):
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
		with open(lfdef, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				type=row[0]
				if (type.find(stationtype) > -1) and (type.find(band) > -1):
					x=(float(row[2])-meanx)
					y=(float(row[1])-meany)
					self.stations['x'][station]=x/1000.0
					self.stations['y'][station]=y/1000.0
					station=station+1

	def assess(self):
		return 1.0
	
#
# UV coverage
#
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
		self.nsources=10000
	
	def construct(self, name='Sources', nsources=10000, radius=1):
		self.name=name
		self.sources={}
		self.sources['x'], self.sources['y']=TelUtils().uniformcircle(nsources, radius)
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
	
	def plot(self):
		plt.clf()
		plt.title('Piercings %s' % self.name)
		plt.xlabel('X (km)')
		plt.ylabel('Y (km)')
		plt.plot(self.piercings['x'], self.piercings['y'], '.')
		plt.axes().set_aspect('equal')
		plt.savefig('Piercings_%s.pdf' % self.name)
	
	def construct(self, sources, array, hiono=400):
		self.name='Piercings_%s_%s' % (sources.name, array.name)
		self.hiono=hiono
		self.npiercings=sources.nsources*array.nstations
		self.piercings={}
		self.piercings['x']=numpy.zeros(self.npiercings)
		self.piercings['y']=numpy.zeros(self.npiercings)
		for source in range(sources.nsources):
			self.piercings['x'][source*array.nstations:(source+1)*array.nstations]=self.hiono*sources.sources['x'][source]+array.stations['x']
			self.piercings['y'][source*array.nstations:(source+1)*array.nstations]=self.hiono*sources.sources['y'][source]+array.stations['y']

	def assess(self):
		return 1.0
	

lowrand=TelArray()
lowrand.construct('LOW_RANDOM', nstations=512, rhalo=20)
lowrand.plot()

low=TelArray()
low.readLOWL1('LOW_L1')
low.plot()

lofar=TelArray()
lofar.readLOFAR('LOFAR')
lofar.plot()
#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()
ts.construct(nsources=33, radius=3.0/35.0)

tp=TelPiercings()
tp.construct(ts,lowrand)
tp.plot()

tplow=TelPiercings()
tplow.construct(ts,low,hiono=200)
tplow.plot()

tplofar=TelPiercings()
tplofar.construct(ts,lofar,hiono=200)
tplofar.plot()

