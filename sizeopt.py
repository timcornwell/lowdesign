from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=6

ntrials=10
nsizesmax=14

saverandom=numpy.zeros(nsizesmax)
saveboolardy=numpy.zeros(nsizesmax)
sizes=numpy.zeros(nsizesmax)

tuv=TelUV()

for isize in range(nsizesmax):
	size=20.0+5.0*isize
	rmax=50*35.0/size
	sizes[isize]=size
	lowrand=TelArray()
	lowrand.circles('LOW_CIRCLE_diameter=%d' % size, nstations=46, nhalo=46, rhalo=rmax)
	lowrand.shakehalo(rshake=5.0)
	lowboolardy=TelArray()
	lowboolardy.readKML('boolardyhandedit.kml')

	tp=TelPiercings()

	random.seed(781490893)
#	ts.construct(nsources=int(nsources*(35.0/size)*(35.0/size)), radius=3.0/size)
	ts.construct(nsources=nsources, radius=3.0/size)
	for trial in range(ntrials):
		tp.construct(ts,lowrand,hiono=300,rmin=2.0)
		if trial==0:
			tp.plot(rmax=rmax)
		saverandom[isize]=max(saverandom[isize],tp.assess(nnoll=100, rmax=rmax, doplot=(trial==0)))
		tp.construct(ts,lowboolardy,hiono=300,rmin=2.0)
		if trial==0:
			tp.plot(rmax=rmax)
		saveboolardy[isize]=max(saveboolardy[isize],tp.assess(nnoll=100, rmax=rmax, doplot=(trial==0)))
	
plt.clf()
plt.plot(sizes, saverandom, color='r')
plt.plot(sizes, saveboolardy, color='g')
plt.title('S Average for fixed source density, r:random array, g:Boolardy array')
plt.xlabel('Size of station (m)')
plt.ylabel('S average')
plt.savefig('savesizes.pdf')
