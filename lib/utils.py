#! /usr/bin/env  python 

#rom numarray import *  Sun Feb 10 09:11:14 EST 2008
#rom numpy import * Sun Feb 10 09:11:14 EST 2008
#mport numarray as N Sun Feb 10 09:11:14 EST 2008
import numpy as N  # Sun Feb 10 09:11:14 EST 2008
import os
import math
import re
import commands

__MYCOMPLEX = type(1.0 + 1.0j)

def uds():
	return unixdatestring()
def unixdatestring():
	(stat,unixdate) = commands.getstatusoutput('date')
	return unixdate

def createlamlist(nlam=5, lamctr=1.0, beta=0.2):

	if nlam == 1:
		return [1.0,]
	else:
		# beta frac  bandwidth
		lamspan = beta * lamctr
		lamlo = lamctr - 0.5 * lamspan
		dlam = lamspan / N.float(nlam - 1)
		lamlist = []
		for i in range(nlam):
			lamlist.append(lamlo + i * dlam)
		return lamlist


def radprof2(a, rstep=4.0):

	n = a.shape[0]
	ctr = int((N.float(n) + 1.0)/2.0) # slicewise
	#print "RADPROF2 center:", ctr
	#print n, ctr
	r_array = euclid(a.shape, c=(ctr-0.5,ctr-0.5))
	#showarray(r_array, "r_array")

	nsteps = int(n/(2*rstep))
	#print "nsteps", nsteps

	profile = N.zeros((nsteps,), type=N.float64)

	for ir in range(nsteps):
		rlo = ir * rstep
		rhi = rlo + rstep
		masklo = N.where(N.greater_equal(r_array, rlo), 1, 0)
		maskhi = N.where(         N.less(r_array, rhi), 1, 0)
		mask = masklo * maskhi
		#showarray(mask, "mask")
		indexarray = N.where(mask == 1)
		profile[ir] = a[indexarray].mean()
	radii = (N.arange(nsteps) + 0.5) * rstep 

	return radii, profile


def showvec(a, str=None):

	if str:
		print str, ": ",

	for c in range(a.shape[0]):
		print " %5.2f" % a[c],
	print


def showarray(a, str=None):

	if str:
		print str, ": ",
	print

	for c in range(a.shape[0]):
		for r in range(a.shape[1]):
			print " %5.2f" % a[r,c],
		print
	print



def test_radprof2():

	n = 6
	C = float(n)*0.5 - 0.5
	ctr = (C,C)
	im = euclid((n,n), c=ctr)
	r, p = radprof2(im, rstep=1.0)
	showarray(im)
	showvec(r, "r")
	showvec(p, "p")


	n = 12
	C = float(n)*0.5 - 0.5
	ctr = (C,C)
	im = euclid((n,n), c=ctr)
	r, p = radprof2(im, rstep=2.0)
	showarray(im)
	showvec(r, "r")
	showvec(p, "p")

	print 
	n = 50
	C = float(n)*0.5 - 0.5
	ctr = (C,C)
	im = euclid((n,n), c=ctr)
	r, p = radprof2(im, rstep=2.0)
	showvec(r, "r")
	showvec(p, "p")



def euclid2(s, c=None):

	if c is None:
		c = (0.5*float(s[0]),  0.5*float(s[1]))

	y, x = N.indices(s)
	r2 = (x - c[0])**2 + (y - c[1])**2

	return r2

#	>>> from utils import *   FFTW fftw style pixel-centered - 
#	>>> s = 4
#	>>> e = euclid((s,s), c=(float(s/2), float(s/2))) 
#	>>> e
#	array([[ 2.82842712,  2.23606798,  2.        ,  2.23606798],
#	       [ 2.23606798,  1.41421356,  1.        ,  1.41421356],
#	       [ 2.        ,  1.        ,  0.        ,  1.        ],
#	       [ 2.23606798,  1.41421356,  1.        ,  1.41421356]])
#	>>> e[2,2]
#	0.0
#	>>> 

def euclid(s, c=None):

	return sqrt(euclid2(s, c=c))


def makedisk(s=None, c=None, r=None, inside=1.0, outside=0.0, grey=None, t=None):
	
	# fft style or sft asymmetric style - center = nx/2, ny/2
	# see ellipseDriver.py for details on symm...

	disk = N.where(euclid2(s, c=c) <= r*r, inside, outside)
	return disk


def ellipse2d(s=None, c=None, ab=None):

	if c is None:
		# fft style or sft asymmetric style - center = nx/2, ny/2
		c = (0.5*float(s[0]),  0.5*float(s[1]))
	if ab is None:
		ab = (0.25*float(s[0]),  0.125*float(s[1]))

	y, x = N.indices(s)
	a2, b2 = (ab[0]*ab[0], ab[1]*ab[1])
	ell2 = ((x - c[0])**2)/a2 + ((y - c[1])**2)/b2

	return ell2


def ellmask(s=None, c=None, ab=None, inside=1.0, outside=0.0):
	
	# fft style or sft asymmetric style - center = nx/2, ny/2
	# see ellipseDriver.py for details on symm...

	ell2 =  ellipse2d(s, c=c, ab=ab)
	mask = N.where(less_equal(ell2, 1.0), inside, outside)
	return mask



def power2d(a):
	if type(a[0,0]) is __MYCOMPLEX:
		r,i  = a.real*a.real, a.imag*a.imag
		a2sum = r.sum() + i.sum()
	else:
		a2 =  a*a 
		a2sum = a2.sum()
	return a2sum

def test_power2d():

	a = arange(6)
	b = arange(6)
	c = a + 1.0j * sin(b)
	a.shape = (2,3)
	b.shape = (2,3)
	c.shape = (2,3)
	print "power2d(a), power2d(sin(b)), power2d(c)", 
	print power2d(a), power2d(sin(b)), power2d(c)




def stripe(a, width=None, value=None, nstripe=None, style=None):
	"""
	puts stripes (no partial stripes) in horiz or vert across a,
	setting stripes to be of numercal value given in 'value'
	"""
	clen,rlen = a.shape

	# n is len of row or col as appropriate
	if style == 'Vert':
		n = rlen
		singlestripe = zeros((rlen, width)) + value
	if style == 'Horiz':
		n = clen
		singlestripe = zeros((width, clen)) + value

	# periodicity requested, in pixels
	period = int(n/nstripe) # assumes all positive... avoids future

	if style == 'Vert':
		for i in range(nstripe):
			a[:, i*period : i*period + width] = singlestripe.copy()
	if style == 'Horiz':
		for i in range(nstripe):
			a[i*period : i*period + width, :] = singlestripe.copy()


#spamming funcs over arrays: 

def sine2d(x,y):
	argu  = 2.0 * math.pi * (x - sine2d.center[0]) / sine2d.spatialWavelen
	return sin(argu)

def cosine2d(x,y):
	argu  = 2.0 * math.pi * (y - cosine2d.center[0]) / cosine2d.spatialWavelen
	return cos(argu)

def coswave2d(x,y):
	argu  = 2.0 * math.pi * (y - coswave2d.center[0]) / coswave2d.spatialWavelen - coswave2d.offset
	return cos(argu)

def parabola2d(x,y):
	"""
	2-tuple 'center' is hardwired from outside me...
	normalization also set from outside
	"""
	return ((x - parabola2d.center[0]) * (x - parabola2d.center[1])  +   \
	        (y - parabola2d.center[0]) * (y - parabola2d.center[1]))

def minus3(x,y):
	"""
	2-tuple 'center' is hardwired from outside me...
	epsilonization also set from outside
	"""
	r = sqrt(((x - minus3.center[0] + 1.0) * (x - minus3.center[1])  +   \
	          (y - minus3.center[0] + 1.0) * (y - minus3.center[1]) + 1.0) + 1.0)
	return (1.0 / (r*r*r))

def kwave2d(x,y):
	argu  = 2.0 * math.pi * (kwave2d.khat[0] * (x - kwave2d.center[0])  +  \
	                         kwave2d.khat[1] * (y - kwave2d.center[1])) /  \
	        kwave2d.spatialWavelen - kwave2d.offset
	return cos(argu)



def test_fromFuncs():

	arrayshape = (64, 64)
	center = (arrayshape[0] - 0.5, arrayshape[1] - 0.5)

	parabola2d.center = center
	quadratic = fromfunction(parabola2d, arrayshape)
	
	sine2d.center = center
	sine2d.spatialWavelen = arrayshape[0]/5
	sine1height = fromfunction(sine2d, arrayshape)

	coswave2d.center = center
	coswave2d.spatialWavelen = arrayshape[0]/5
	coswave2d.offset =  math.pi / 4.0
	coswave2dheight = fromfunction(coswave2d, arrayshape)

	minus3.center = (64/2 + .5, 64/2 + .5)
	m3 = fromfunction(minus3, arrayshape)
	SimpleFitsWrite(fn="testdata/m3.fits", data=m3, clobber="y")

	arrayshape = (100, 100)
	center = (arrayshape[0]/2, arrayshape[1]/2)
	for offset in (0, 30, 45, 90):
		for angle in (0, 30, 45, 90):
			kwave2d.spatialWavelen = arrayshape[0] / 4
			kwave2d.center = center
			kwave2d.offset = float(offset) * pi/180.0
			rangle = angle*pi/180.0
			kwave2d.khat = array((sin(rangle), cos(rangle)))
			print kwave2d.khat
			kwavedata = fromfunction(kwave2d, arrayshape)
			SimpleFitsWrite(fn="testdata/kwav_16_%doff_%dang.fits"%(offset,angle), 
			                data=kwavedata, clobber="y")


	return None


def centralSection(a, npoints=None):

	if npoints%2 == 1:
		nn = npoints + 1
	else:
		nn = npoints

	sa = a.shape
	b = a[(sa[0]/2 - nn/2) : (sa[0]/2 + nn/2), \
	      (sa[1]/2 - nn/2) : (sa[1]/2 + nn/2)].copy()
	return b


def  rebin(a = None, rc=(2,2), verbose=None):

	"""  
	anand@stsci.edu
	Perform simple-minded flux-conserving binning... clip trailing
	size mismatch: eg a 10x3 array binned by 3 results in a 3x1 array
	"""

	r, c = rc

	R = a.shape[0]
	C = a.shape[1]

	nr = int(R / r)
	nc = int(C / c)

	b = a[0:nr, 0:nc].copy()
	b = b * 0

	for ri in range(0, nr):
		Rlo = ri * r
		if verbose:
			print "row loop"
		for ci in range(0, nc):
			Clo = ci * c
			b[ri, ci] = add.reduce(a[Rlo:Rlo+r, Clo:Clo+c].copy().flat)
			if verbose:
				print "    [%d:%d, %d:%d]" % (Rlo,Rlo+r, Clo,Clo+c),
				print "%4.0f"  %   add.reduce(a[Rlo:Rlo+r, Clo:Clo+c].copy().flat)
	return b


def test_rebin():

	R = 4
	C = 6
	rc = (1,2)
	a = zeros((R,C), float32)
	for i in range(R):
		for j in range(C):
			a[i,j] = 10 * i + 1 * j
	print a, "\n\n"
	b = rebin(a=a, rc=rc, verbose=0xff)
	print b

	rc = (2,1)
	b = rebin(a=a, rc=rc, verbose=0xff)
	print b

	rc = (2,2)
	b = rebin(a=a, rc=rc, verbose=0xff)
	print b


import os
import commands

def newdir(dir):

	# too much typing -- like an alias
	cmd = commands.getstatusoutput
	# if dir exists, wipe it and recreate empty one
	if os.path.isdir(dir):
		print 'Warning: %s  exists: removing it recursively now...\n' % dir
		(s,o) = cmd("rm -rf %s"  %  dir)
	os.mkdir(dir)
	return None

def newdir_ifrequired(dir):
	# if dir exists, let caller know
	# otherwise make it and inform caller

	if os.path.isdir(dir):
		print 'Info: %s  exists: using this directory...\n' % dir
	else:
		print 'Info: %s does not exist: creating this directory...\n' % dir
		os.mkdir(dir)
	return None


def test_newdir():

	dn = '/tmp/newdirtest'
	newdir(dn)
	newdir(dn)
	newdir_ifrequired(dn)



def rmifexist(fn=None):

	"""  function to clobber a file if it exists  anand@stsci.edu
	"""

	# too much typing -- like an alias
	cmd = commands.getstatusoutput
	# if dir exists, wipe it and recreate empty one
	if os.path.isfile(fn):
		##print '     Warning: %s: file or symlink %s  exists:' % (__name__, fn)
		##(s,o) = cmd("ls -l %s"  %  fn)
		##print '    ' + o
		##print '    clobbered it' 
		(s,o) = cmd("rm %s"  %  fn)
	return None


def test_rmifexist():

	cmd = commands.getstatusoutput
	fn = "/tmp/zztop_armadillo"
	(s,o) = cmd("touch %s"  %  fn)
	rmifexist(fn=fn)

	fn = "/tmp/zztop"
	rmifexist(fn=fn)



def flop(a=None):
	"""  
	Flips a 2d data array to be ft-space-compatible 
           Q4 Q3 ->  Q2 Q1
           Q1 Q2     Q3 Q4
	call with a.copy() if you want to retain a undamaged
		   anand@stsci.edu
	"""


	n = a.shape[0]
	m = a.shape[1]

	b = a - a

	# Q1 -> q,  Q3 -> Q1,  q -> Q3
	#
	quad = a[:n/2, :m/2].copy()
	b[:n/2, :m/2] = a[n/2:, m/2:].copy()
	b[n/2:, m/2:] = quad.copy()

	# Q2 -> q,  Q4 -> Q2,  q -> Q4
	#
	quad = a[n/2:, :m/2].copy()
	b[n/2:, :m/2] = a[:n/2, m/2:].copy()
	b[:n/2, m/2:] = quad.copy()

	return b



def test_flop():

	r = 6
	c = 8
	a = zeros((r,c), float32)
	for i in range(r):
		for j in range(c):
			a[i,j] = 10 * i + 1 * j
	print a, "\n\n"
	a = Flop(a)
	print a


def binarize(a=None, cut=None, t=N.int32):

	"""
	From perry@stsci.edu Tue Jan 20 14:42:07 2004
	Date: Tue, 20 Jan 2004 13:25:11 -0500
	From: Perry Greenfield <perry@stsci.edu>
	To: Anand Sivaramakrishnan <anand@stsci.edu>, Phil Hodge <hodge@stsci.
	edu>
	Subject: RE: python question
	
	why not 
	
	(x>=42.).astype(N.int32)
	
	(or whatever type you want it as if boolean isn't acceptable.
	If N.int32 is what you wanted:
	
	(x>=42.)*1
	
	would also work. It's possible that the where solution is 
	faster thought I suspect there aren't great differences.
	
	"""

	b = (a >= cut) * 1
	return b


def test_binarize():

	arr = arange(25)
	arr.shape = (5,5)
	c = 7.0
	b = binarize(a=arr, cut=c)
	print arr
	print "cut at %f"  %  c
	print b


def total2d(a=None):

	t = 0
	for i in range(a.shape[0]):
		for j in range(a.shape[1]):
			t = t + a[i,j]
	return t

def localinfo():
	print "Local info ---"
	# too much typing -- like an alias
	cmd = commands.getstatusoutput
	(s,o) = cmd('uname -snrvm')
	print "    uname -snrvm: " + o
	(s,o) = cmd('date')
	print "    today's date: " + o
	print "    current working directory: " + os.getcwd()
	print 
	return None

def test_localinfo():

	localinfo()



def binarymaskmean(a=None, mask=None):
	"""
	# where mask is good add up data
	# where mask is  0 # ignore the sum
	# return summed data over mask / mask area
	"""

	amask  = a * mask
	NN = mask.sum()
	return amask.sum()/N.float(NN)


def maskdiv(a=None, mask=None, irrelevantValue=0):
	# where mask is good data divide by it
	# where mask is zero replace the a/mask elements with 
	# the irrelevant data numerical value, typically zero itself.

	return N.where(equal(mask,0), irrelevantValue, a/mask)


def test_maskdiv():
	R, C = (4,4)
	a = zeros((R,C), N.float32)

	for i in range(R):
		for j in range(C):
			a[i,j] = 10 * i + 1 * j
	print "a"
	print a, "\n\n"
	# eg a fkat file
	b = zeros((R,C), N.float32)
	for i in range(R):
		for j in range(C):
			b[i,j] = 10 * i + 1 * j
	# diag stripe of zeroes
	for i in range(R):
		b[i,i] = 0.0
	print "b"
	print b, "\n\n"
	a = maskdiv(a=a, mask=b)
	print "a"
	print a, "\n\n"

	output = """
		a
		[[  0.   1.   2.   3.]
		 [ 10.  11.  12.  13.]
		 [ 20.  21.  22.  23.]
		 [ 30.  31.  32.  33.]] 
		
		
		b
		[[  0.   1.   2.   3.]
		 [ 10.   0.  12.  13.]
		 [ 20.  21.   0.  23.]
		 [ 30.  31.  32.   0.]] 
		
		
		Warning: Encountered invalid numeric result(s)  in divide
		Warning: Encountered divide by zero(s)  in divide
		a
		[[ 0.  1.  1.  1.]
		 [ 1.  0.  1.  1.]
		 [ 1.  1.  0.  1.]
		 [ 1.  1.  1.  0.]] 

	"""




from SimpleFits import *

def phasevariance(opd=None, pupil=None, usepupil=None, opd2rad=None, report=None):
	"""
		for 0/1 pupils only for now...
		omit opd2rad if already working in radians
		if you want variance over pupil set usepupil to anything except None
		if you want variance over whole array omit pupil and usepupil
		return zero mean phase function (radians), and the variance thereof
	"""
	if usepupil:
		if opd2rad == None:
			phi  = opd * pupil
		else:
			phi  = opd * pupil * opd2rad
		NN = add.reduce(pupil.flat)
	else:
		if opd2rad == None:
			phi = opd 
		else:
			phi = opd * opd2rad
		NN = opd.shape[0] * opd.shape[1]

	meanphi = phi.sum() / NN
	phizm = phi - meanphi
	phisq = phizm * phizm
	varphi = phisq.sum() / (NN-1)
	if report:
		if usepupil:
			print "  Using Pupil  ",
		else:
			print "  Using Full Array...  ",
		print "  mean %+.2e  var %.2e"  %   (meanphi, varphi),
		print "  total %.2e, support %.0f"  %  (phi.sum(), NN)
		if opd2rad:
			print "                          unit opd = %.4f radians" %  opd2rad
			print "                          Marechal hit = %.1f%%" %  \
			               (int(around(100.0 * varphi)))

	return phizm, varphi


def test_phasevariance():


	flist = ["aeos32l1no0000phase.fits",
	         "aeos32l1no0000SFWFS_InBand.fits",
	         "aeos32l1no0000SFWFS_Resid.fits",
		     "aeos32l1no0000Para_InBand.fits",
	         "aeos32l1no0000Para_Resid.fits"
			 ]

	for f in flist:
		print "testdata/"+f
		phi = SimpleFitsRead(fn="testdata/"+f).data
		phi, sigsq = phasevariance(opd=phi, report=0xffff)
		print "  newPhiMean = %g  variance = %g\n" % (phi.mean(), sigsq)


	#pfn = 'testdata/xaopiaperture.fits'
	#opd = (SimpleFitsRead(fn=ifn, report=42).data) * dbl
	## lambda =  H band center, opd in microns
	#opd2rad = (1.0 / 1.00) * 2.0 * math.pi # lambda =  1um, opd in microns

	#pupil = (SimpleFitsRead(fn=pfn, report=42).data) * dbl

	#print 
	#print 
	#print "raw data, repeat to see zero mean, low total second time?"
	#print 
	#phi, sig = phasevariance(opd=opd, report=0xffff)
	#phi, sig = phasevariance(opd=phi, report=0xffff)

	#print 
	#print 
	#print "opd2rad data, repeat to see zero mean, low total second time?"
	#print 
	#phi, sig = phasevariance(opd=opd, report=0xffff, opd2rad=opd2rad)
	#phi, sig = phasevariance(opd=phi, report=0xffff)

	#print 
	#print 
	#print "opd2rad data with pupil, repeat to see zero mean, low total second time?"
	#print 
	#phi, sig = phasevariance(opd=opd, report=0xffff, pupil=pupil, usepupil='y', \
	#                         opd2rad=opd2rad)
	#phi, sig = phasevariance(opd=phi, report=0xffff, pupil=pupil, usepupil='y')





def mandtparser(s):

	# you have to de-comma-ize the numerical credit or debit amounts
	# before using this program.
	# this gets rid of the trailing blank(s) after "," delimiter...
	comma = re.compile(', *')
	# now this gets rid of repeated blanks in the items already split() out
	squeezeblanks = re.compile('  *')
	dequote = re.compile('"')
	q = comma.split(s)
	l = []
	for x in q:
		y = squeezeblanks.sub(' ', x)
		z = dequote.sub('', y)
		if z != '\n' and z != '':
			l.append(z)
	return l


def test_mandparser():
	p = re.compile(r'\W+')
	q = p.split("Now is the winter of our discontent made glorious summer by this son of York")
	print q
	

	txt = 'a, b b,   c c c, "D  D         ..."'
	tokens = mandtparser(txt)
	for t in tokens:
		print t
	return None



def rotsub(P, M):
	""" flip and rotate two images then sub m-fliprot from p 
	I assume the central pixel is at NN/2, NN/2 in an 
	even-number-sided square NNxNN, so I trim the first row and col
	off to 'centralize' the zero-tilt pixel before rotsubbing it
	"""
	p = P[1:, 1:].copy()
	m = M[1:, 1:].copy()

	rs =  p - m[::-1, ::-1]
	return rs



def rotsub_test():

	"""
	flip-rotate two images and subtract them for MGS-style signal viewing
	anand@stsci.edu April 2005 
	"""

	if sys.argv[3:]:

		of = sys.argv[3]
		print "output file will be " + of
		localinfo()
		lo = 0
		hi = 1024
		P = SimpleFitsRead(fn=sys.argv[1], report=0x55, verbose=0xff).data[lo:hi, lo:hi].copy()
		M = SimpleFitsRead(fn=sys.argv[2], report=0x55, verbose=0xff).data[lo:hi, lo:hi].copy()
	
		D = rotsub(P, M)
		SimpleFitsWrite(fn=of, data=D, verbose=0xff, clobber='y')

	else:
		print """
			"__moi__ DefocP.fits DefocM.fits  DefocSignal.fits
			
			which creates DefocSignal.fits
		"""



def fact(n):
	""" factorial function a la guido
	"""
	result = 1
	while n > 1:
		result = result * n
		n = n - 1
	return result

def test_fact():
	for n in range(6):
		print " %d! = %d"  %  (n, fact(n))

def linmac(c600=None):

	whocares, machine =  commands.getstatusoutput('uname -m')
	if machine == "Power Macintosh":
		if c600:
			ddir = '/Users/anand/C600Linux/galois/data/'
		else:
			ddir = '/Users/anand/data/'
	else:
		ddir = '/home/anand/data/'
	return ddir


def makeLambdas(lamCentral=1.65e-6, fracBW=0.2, Nintervals=2):
	lamlo, lamhi  = (lamCentral * (1 - fracBW/2.0), lamCentral * (1 + fracBW/2.0))
	dlam = (lamhi - lamlo)
	lamstep = dlam / Nintervals
	lamEps = lamstep / 10000.0
	lambdas = arange(lamlo, lamhi+lamEps, lamstep)
	return lambdas


def scientific(x, decplaces=1):

	sign = 1
	if x < 0:
		sign = -1
		
	X = abs(x)
	tenpow = floor(log10(X))
	factor = pow(10.0, tenpow)
	sigs = X / factor

	rounder = pow(10, decplaces)
	rsigs = around(sigs*rounder)
	dsigs = sign * rsigs / rounder

	#print "X %g, tenpow %g, factor %g, sigs %g: rounder %g, rsigs %g dsigs %g" % \
	#(X, tenpow, factor, sigs, rounder, rsigs, dsigs )

	return r"%g \times 10^{%d}" % (dsigs, tenpow)

def test_scientific():

	print
	for x in (1.0e-0, -2.6e-7, -2.3e-7, 3.142e2, 3.542e2):
		print None, scientific(x), x
	print

	for x in (1.0e-0, -2.6e-7, -2.3e-7, 3.142e2, 3.542e2):
		print 0, scientific(x, decplaces=0), x
	print

	for x in (1.0e-0, -2.6e-7, -2.3e-7, 3.142e2, 3.542e2):
		print 6, scientific(x, decplaces=6), x
	print

	for x in (1.0e-0, -2.6e-7, -2.3e-7, 3.142e2, 3.542e2):
		print 4, scientific(x, decplaces=4), x
	print






def quantumdisplay(a, str=None, nc=5):

	if str:
		print
		print "**",  str

	af = a.flat
	ctr = 0
	for j in range(len(af)):
		print "%6.2f" % af[ctr],
		ctr += 1
		if (ctr%9 == 0):
			print
	print


def quantize(a, nquanta=10, verbose=None):
	"""
	Create nquanta levels above the ground state, maintaining
	the canonical measure on Dedekind's famous cuts, as measure
	was defined first by Lebesgue.
	Thus asking for 10 levels on an array possessing values
	only between 0 and 100 produces zero, and ten greys, 
	viz. 0, 10, 20,... 100, which is 11 levels.
	"""

	az = a - a.min()
	span = a.max() - a.min()
	if verbose is not None:
		quantumdisplay(az, str="AZ")

	azmax = span
	azUnitspan = az / azmax
	if verbose is not None:
		quantumdisplay(azUnitspan, str="AZ unit span")

	azscaled = azUnitspan * nquanta
	if verbose is not None:
		quantumdisplay(azscaled , str="AZ unit span, scaled")

	azint = floor(azscaled + 0.5)
	if verbose is not None:
		quantumdisplay(azint , str="AZ unit span, scaled and floored")

	azunscaled = azint / N.float(nquanta)
	if verbose is not None:
		quantumdisplay(azunscaled , str="AZ unit span, unscaled")

	azmaxed = azunscaled * azmax
	if verbose is not None:
		quantumdisplay(azmaxed , str="AZ unit span, unscaled, maxpanded")

	adezeroed = azmaxed + a.min()
	if verbose is not None:
		quantumdisplay(adezeroed , str="AZ unit span, unscaled, maxpanded, dezeroed")

	return  adezeroed


def test_quantize():

	a = arange(26)
	a = a - 5
	h = 4.0
	a = a * h
	quantumdisplay(a, str = "h(A - 5.0)")

	q = quantize(a, nquanta=10)
	quantumdisplay(q, str="q")



def showarray(a, str=None):

	if str:
		print str, ": ",
	print
	for c in range(a.shape[0]):
		for r in range(a.shape[1]):
			print " %5.2f" % a[r,c],
		print
	print


# n-D array AA of r values
# make sure your r values in the interp array cover your input
# range of r values, though... esp at high end
def FofR(f, r, AA):

	if r.max() < AA.max():
		print "Invalid range of interpolation: r.max < AA.max"
		print "utils.FofR may crash"

	s = AA.shape
	A = AA.flat.copy()

	indexarray = N.searchsorted(r, A) - 1
	frac = (A - r[indexarray]) / (r[indexarray+1] - r[indexarray])
	C = f[indexarray] + frac * (f[indexarray+1] - f[indexarray])
	C.shape = s
	return C



def test_FofR():

	# f(r) represented numerically at various values of r...
	# a small arbitrary example used here.  
	r = N.array((-1.0,  0.0, 1.0, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 10.0, 15.0))
	f = r + 10.0
	print "    r = ", r
	print " f(r) = ", f
	print
	

	n = 4
	s = (n,n)

	# cen = N.float(n)/2.0 - 0.5
	# c = (cen,cen)
	# SS = euclid(s, c=c)


	SS = N.arange(16)
	AA = SS + 0.1 * N.sin(N.pi * SS/7.5)
	AA.shape = (4,4)
	SS.shape = (4,4)
	showarray(SS, str=" SS")
	showarray(AA, str=" AA")
	A = AA.flat.copy()

	indexarray = N.searchsorted(r, A) - 1
	print " indexarray = ", indexarray
	frac = (A - r[indexarray]) / (r[indexarray+1] - r[indexarray])
	C = f[indexarray] + frac * (f[indexarray+1] - f[indexarray])

	C.shape = AA.shape
	showarray(C, str="\n C")


	print" ...and all in one call without explanation...:"
	C = FofR(f, r, AA)
	showarray(C, str="\n C")



def l_z(F, f, z=0.0):

	return F * N.sqrt(1.0 + 0.25/(f*f) - 2.0*z/F + z*z/(F*F))

def p2vDefoc_vs_z(D=5.0, f=15.8, z=1.0e-3, wavelen=None):

	m = 1.0
	mm = 1.0e-3
	um = 1.0e-6
	nm = 1.0e-9
	rad = 1.0


	M2RAD = 2.0 * N.pi / wavelen
	F = f * D

	z0 = 0.0
	OPDctr = F * (1 - z0/F)
	OPDedg = l_z(F, f, z=z0)
	P2V = OPDctr - OPDedg

	P2V0 = P2V

	OPDctr = F * (1 - z/F)
	OPDedg = l_z(F, f, z=z)
	P2V = OPDctr - OPDedg - P2V0
	return (P2V,  P2V/wavelen, P2V*M2RAD, P2V/nm)




def test_p2vDefocus_vs_z():

	print """

	Lens dia D   f-ratio f   focal length F
	Move detector (plane perp to Z axis) by z towards the lens
	What is the P-V OPD of the wavefront forming the PSF for
	a given z translation of the detector?

	"""

	m = 1.0
	mm = 1.0e-3
	um = 1.0e-6
	nm = 1.0e-9
	rad = 1.0

	wavelen = 1.65 * um

	D = 1.0 * m
	f = 50.0

	D = 8.1 * m
	f = 64.0

	F = D * f
	print "\tD = ", D, "m"
	print "\tF = ", F, "m"
	print "\tf-ratio = ", f
	print "\twavelen = ", wavelen/um, "um"
	print "\t10nm = ", 10.0*nm/wavelen, "waves"


	print"""
   z       P2V         P2V          P2V         P2V
   mm       m          wav          rad         nm 
	"""

	for z_mm in (0.0, 1.0, 6.758, 2*1.65e-6*64*64/mm,  30.0, 100.0):

		z = z_mm * mm
		(P2V,  P2Vwav, P2Vrad, P2Vnm) = p2vDefoc_vs_z(D=D, f=f, z=z, wavelen=wavelen)

		print " %3.0f   %+.3e   %+.3e   %+.3e   %+.1e"  %  \
		      (z/mm,  P2V,  P2Vwav, P2Vrad, P2Vnm)
	print

def gauss(sigma=5.0, fwhm=None, normalize='yes',  \
	               c=None, s=(128,128), t=N.float32):
	"""
	  makes isotropic Gaussian pdf ie normalized to unity (anands@stsci.edu)
	  normalize=yes returns the value of a pdf even if total power in the
	  array is less than unity.  Otherwise exp-r^2/sigma^2 is returned.
	  anand@stsci.edu
	"""

	# Set center of Gaussian in pixel coordinates...
	if c is None:
		c = (0.5*N.float(s[0]),  0.5*N.float(s[1]))

	if fwhm:
		sigma = fwhm/2.35

	y, x = N.indices(s)
	r2 = (x - c[0])**2 + (y - c[1])**2
	g  = N.exp(-1.0 * r2 / (sigma*sigma))

	if normalize == 'yes':
		g = g / (N.pi * sigma * sigma)

	return g



def testgauss():

	sigma=100.0
	FWHM = 64.0
	s=(400,400)
	c = (N.float(s[0])/2.0 - 0.5, N.float(s[1])/2.0 - 0.5)
	g = gauss(sigma=sigma, c=c, s=s, normalize="no")
	SimpleFitsWrite(fn="testdata/gauss.fits", data=g, verbose='v', clobber='y')
	g = gauss(sigma=sigma, c=c, s=s, normalize="yes")
	print g.sum()
	SimpleFitsWrite(fn="testdata/gaussN.fits", data=g, verbose='v', clobber='y')


def center (x):
    """Shift the contents of x by half its width.

    x is assumed to be rank 2, and the length of each axis is assumed
    to be even.  x will be modified in-place.
    """

    shape = x.shape
    nx = shape[1]
    ny = shape[0]
    hnx = nx // 2
    hny = ny // 2

    temp = x[0:hny,0:hnx].copy()
    x[0:hny,0:hnx] = x[hny:ny,hnx:nx].copy()
    x[hny:ny,hnx:nx] = temp

    temp = x[0:hny,hnx:nx].copy()
    x[0:hny,hnx:nx] = x[hny:ny,0:hnx].copy()
    x[hny:ny,0:hnx] = temp
    
def frange(start, end=None, inc=None):
    """A range function, that does accept float increments..."""
   
    if end == None:
        end = start + 0.0
        start = 0.0
    else: start += 0.0 # force it to be a float

    if inc == None:
        inc = 1.0
    count = int(math.ceil((end - start) / inc))

    L = [None,] * count
    L[0] = start
    for i in xrange(1,count):
        L[i] = L[i-1] + inc
    return L


