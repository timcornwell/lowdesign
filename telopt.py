import numpy
import scipy
import matplotlib
import random
import matplotlib.pyplot as plt
import csv
import zernike
from scipy import linalg

random.seed(781490893)


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
	
	def construct(self, name='Stations', rhalo=40, rcore=1.0, nstations=512, nhalo=45, nantennas=256, fobs=1e8, diameter=0.035):
		self.name=name
		self.rhalo=rhalo
		self.rcore=rcore
		self.nstations=nstations
		ncore=self.nstations-nhalo
		self.fobs=fobs
		self.diameter=diameter
		self.stations={}
		self.stations['x'], self.stations['y']=TelUtils().uniformcircle(self.nstations, self.rhalo)
		self.stations['x'][:ncore], self.stations['y'][:ncore]=TelUtils().uniformcircle(ncore, self.rcore)

	def circles(self, name='Stations', rhalo=40, rcore=1.0, nstations=512, nhalo=52, fobs=1e8, diameter=0.035):
		self.name=name
		self.rhalo=rhalo
		self.rcore=rcore
		self.nstations=nstations
		ncore=self.nstations-nhalo
		self.fobs=fobs
		self.diameter=diameter
		if nhalo==60:
			self.nrings=4
			self.r=[0.0, rhalo/3.0, 2*rhalo/3.0, rhalo]
			self.nonring=[1, 9, 21, 29]
		elif nhalo==44:
			self.nrings=4
			self.r=[0.0, rhalo/3.0, 2*rhalo/3.0, rhalo]
			self.nonring=[1, 7, 13, 23]
		else:
			nhalo=185
			self.nrings=7
			self.r=[0.0, rhalo/6.0, 2*rhalo/6.0, 3.0*rhalo/6.0, 4.0*rhalo/6.0, 5.0*rhalo/6.0, rhalo]
			self.nonring=[1, 9, 19, 27, 35, 43, 51]
		self.nstations=nhalo
		self.stations={}
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
# 		self.stations['x'][:ncore], self.stations['y'][:ncore]=TelUtils().uniformcircle(ncore, self.rcore)
		station=0
		for ring in range(self.nrings):
			dphi=2*numpy.pi/self.nonring[ring]
			phi=0.0
			for spoke in range(self.nonring[ring]):
				self.stations['x'][station]=self.r[ring]*numpy.cos(phi)
				self.stations['y'][station]=self.r[ring]*numpy.sin(phi)
				phi=phi+dphi
				station=station+1

	def shakehalo(self, rshake):
		for station in range(self.nstations):
			r=numpy.sqrt(self.stations['x'][station]*self.stations['x'][station]+self.stations['y'][station]*self.stations['y'][station])
			if r>self.rcore:
				self.stations['x'][station]=self.stations['x'][station]+random.uniform(-rshake,rshake)
				self.stations['y'][station]=self.stations['y'][station]+random.uniform(-rshake,rshake)
			r=numpy.sqrt(self.stations['x'][station]*self.stations['x'][station]+self.stations['y'][station]*self.stations['y'][station])
			if r>self.rhalo:
				self.stations['x'][station]=self.rhalo*self.stations['x'][station]/r
				self.stations['y'][station]=self.rhalo*self.stations['y'][station]/r
				

	def readLOWBD(self, name='LOWBD', rcore=0.0, l1def='SKA-low_config_baseline_design_arm_stations_2013apr30.csv'):
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
		scale=0.001
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				x=scale*(float(row[1])-meanx)
				y=scale*(float(row[0])-meany)
				r=numpy.sqrt(x*x+y*y)
				if r>rcore:
					self.nstations=self.nstations+1
		
		self.stations['x']=numpy.zeros(self.nstations)
		self.stations['y']=numpy.zeros(self.nstations)
		station=0
		with open(l1def, 'rU') as f:
			reader = csv.reader(f)
			for row in reader:
				x=scale*(float(row[1])-meanx)
				y=scale*(float(row[0])-meany)
				r=numpy.sqrt(x*x+y*y)
				if r>rcore:
					self.stations['x'][station]=x
					self.stations['y'][station]=y
					station=station+1
		
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
					station=station+1
		
	def readLOFAR(self, name='LOFAR', stationtype='S', band='HBA', lfdef='LOFAR.csv', lat=52.7):
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
		plt.savefig('%s.pdf' % self.name)
	
	def construct(self, sources, array, rmin=1, hiono=400):
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
		for source in range(sources.nsources):
			self.piercings['x'][source*nstations:(source+1)*nstations]=self.hiono*sources.sources['x'][source]+outside['x']
			self.piercings['y'][source*nstations:(source+1)*nstations]=self.hiono*sources.sources['y'][source]+outside['y']

	def assess(self, rmax=70.0, nnoll=20,doplot=True):
		A=numpy.zeros([self.npiercings, nnoll])
		for piercing in range(self.npiercings):
			x=self.piercings['x'][piercing]
			y=self.piercings['y'][piercing]
			r=numpy.sqrt(x*x+y*y)
			phi=numpy.arctan2(y,x)
			if(r<rmax):
				for noll in range(nnoll):
					A[piercing,noll]=zernike.zernikel(noll,r/rmax,phi)
		Covar_A=numpy.zeros([nnoll, nnoll])
		for nnol1 in range(nnoll):
			for nnol2 in range(nnoll):
				Covar_A[nnol1,nnol2]=numpy.sum(A[...,nnol1]*A[...,nnol2])
		U,s,Vh = linalg.svd(Covar_A)
		shalf=s[nnoll/2-1]
		if doplot:
			plt.clf()
			plt.title('%s rmax=%d shalf=%f' % (self.name, rmax, shalf))
			plt.xlabel('Singular vector index')
			plt.ylabel('Singular value')
			plt.semilogy(s, '.')
			plt.savefig('%s_rmax=%d_SVD.pdf' % (self.name, rmax))
		return shalf
	