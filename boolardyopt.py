from telopt import *

#
# Radius=1/2 FWZ ~ FWHM
#
ts=TelSources()

nsources=12
nnoll=60

ntrials=1
rmax=20

random.seed(781490893)

tuv=TelUV()

rhalo=40.0*numpy.sqrt(12.0/46.0)
lowrand12=TelArray()
lowrand12.randomBoolardy('LOW_BOOLARDYRANDOM12', nstations=12, nhalo=12, rhalo=rhalo)
lowrand12.plot(rmax=rhalo)
lowrand12.save('LOW_BOOLARDYRANDOM12_OPT.csv')
print lowrand12.mst()
tuv.construct(lowrand12);tuv.plot()
scircles12=numpy.zeros(nnoll)
for trial in range(ntrials):
	diameter=35.0*numpy.sqrt(46.0/12.0)
	ansources=int(nsources*(35.0/diameter)*(35.0/diameter))
	ts.construct(nsources=ansources, radius=3.0/diameter)
	tp=TelPiercings()
	tp.construct(ts,lowrand12,hiono=300,rmin=0)
	if trial==0:
		print rhalo, diameter, ansources
		tp.plot(rmax=20.0)
	scircles12=scircles12+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rhalo, doplot=(trial==0))	

rhalo=40.0*numpy.sqrt(46.0/46.0)
lowrand46=TelArray()
lowrand46.randomBoolardy('LOW_BOOLARDYRANDOM46', nstations=46, nhalo=46, rhalo=rhalo)
lowrand46.plot(rmax=rhalo)
lowrand46.save('LOW_BOOLARDYRANDOM46_OPT.csv')
print lowrand46.mst()
tuv.construct(lowrand46);tuv.plot()
scircles46=numpy.zeros(nnoll)
for trial in range(ntrials):
	diameter=35.0*numpy.sqrt(46.0/46.0)
	ansources=int(nsources*(35.0/diameter)*(35.0/diameter))
	ts.construct(nsources=ansources, radius=3.0/diameter)
	tp=TelPiercings()
	tp.construct(ts,lowrand46,hiono=300,rmin=0)
	if trial==0:
		print rhalo, diameter, ansources
		tp.plot(rmax=40.0)
	scircles46=scircles46+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rhalo, doplot=(trial==0))
	

rhalo=40.0*numpy.sqrt(185.0/46.0)
lowrand185=TelArray()
lowrand185.randomBoolardy('LOW_BOOLARDYRANDOM185', nstations=185, nhalo=185, rhalo=rhalo)
lowrand185.plot(rmax=rhalo)
lowrand185.save('LOW_BOOLARDYRANDOM185_OPT.csv')
print lowrand185.mst()
tuv.construct(lowrand185);tuv.plot()
scircles185=numpy.zeros(nnoll)
for trial in range(ntrials):
	diameter=35.0*numpy.sqrt(46.0/185.0)
	ansources=int(nsources*(35.0/diameter)*(35.0/diameter))
	ts.construct(nsources=ansources, radius=3.0/diameter)
	tp=TelPiercings()
	tp.construct(ts,lowrand185,hiono=300,rmin=0)
	if trial==0:
		print rhalo, diameter, ansources
		tp.plot(rmax=80.0)
	scircles185=scircles185+(1.0/float(ntrials))*tp.assess(nnoll=nnoll, rmax=rhalo, doplot=(trial==0))
	
plt.clf()
plt.semilogy(scircles12, color='r')
plt.semilogy(scircles46, color='g')
plt.semilogy(scircles185, color='b')
plt.title('Sqrt(Singular values) (68m, 35m, 17m)')
plt.xlabel('Singular value index')
plt.ylabel('Sqrt(Singular value)')
plt.axes().set_ylim([1e-5,1e1])
plt.savefig('Singularvalue.pdf')
