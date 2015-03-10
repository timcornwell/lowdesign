from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=10

ntrials=10
nsizesmax=17

shalfrandom=numpy.zeros(nsizesmax)
sizes=numpy.zeros(nsizesmax)

tuv=TelUV()

for isize in range(nsizesmax):
	size=20.0+5.0*isize
	rmax=50*35.0/size
	sizes[isize]=size
	lowrand=TelArray()
	lowrand.circles('LOW_CIRCLE_diameter=%d' % size, nstations=46, nhalo=46, rhalo=rmax)
	lowrand.shakehalo(rshake=5.0)
# 	lowrand.plot()
# 	tuv.construct(lowrand);tuv.plot()

	tp=TelPiercings()

	random.seed(781490893)
	ts.construct(nsources=int(nsources*(35.0/size)*(35.0/size)), radius=3.0/size)
	for trial in range(ntrials):
		tp.construct(ts,lowrand,hiono=300,rmin=2.0)
		if trial==0:
			tp.plot(rmax=rmax)
		shalfrandom[isize]=max(shalfrandom[isize],(1.0/float(ntrials))*tp.assess(nnoll=100, rmax=rmax, doplot=(trial==0)))
	
plt.clf()
plt.plot(sizes, shalfrandom, color='r')
plt.title('Shalf for 10 sources, random array')
plt.xlabel('Size of station (m)')
plt.ylabel('Shalf')
plt.savefig('shalfsizes.pdf')
