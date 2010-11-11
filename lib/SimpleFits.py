#! /usr/bin/env python

import numpy as N
from utils import rmifexist
import types
import pyfits
import string

def addkwarray(ctrstr, kwdfmt, a1d, cmntfmt, kwdict, kwlist, verbose=None):

	kvc = (ctrstr, a1d.shape[0], "number of entries")
	for i in range(a1d.shape[0]):

		kvc = (kwdfmt%i, a1d[i], cmntfmt%i)
		addto_kwds(kvc, kwdict, kwlist, verbose=verbose)



def rm_kwd(kk, kwdict, kwlist, verbose=None):

	k = kk.upper()

	if kwdict.has_key(k):
		del kwdict[k]
	if k in kwlist:
		kwlist.remove(k)

	if verbose:
		print k, "removed"


def addto_kwds(kvc, kwdict, kwlist, verbose=None):
	kk, v, c = kvc
	k = kk.upper()

	if kwdict.has_key(k):
		del kwdict[k]
	if k in kwlist:
		kwlist.remove(k)

	if c:
		kwdict[k] = (v,c)
	else:
		kwdict[k] = v

	kwlist.append(k)

	if verbose:
		print kvc

class ExtFitsRead:

	def __init__(self, fn=None, extnum=0, report=None, verbose=None):
		"""reads a simple fits file from disk, 2d """
		try:
			if verbose:
				print "Trying to open " + "'"+fn+"'"
			fimg = pyfits.open (fn)
			self.hdr = fimg[extnum].header
			self.data = fimg[extnum].data
			self.fn = fn
		except IOError:
			print ': FATAL: ReadSimpleFits: Cannot open ' + "'"+fn+"'"

		self.dimension = self.hdr['NAXIS']
		self.numCols = self.hdr['NAXIS1']
		self.numRows = self.hdr['NAXIS2']

		self.immax = self.data.flat.max()
		self.immin = self.data.min()

		if report or verbose:
			print "Opened " + '"'+fn+'"'
		if verbose:
			for ascard in self.hdr.ascardlist():
				print ascard

		return None



class SimpleFitsWrite:

	def __init__(self, fn=None, data=None, kwdict=None, kwlist=None, verbose=None, clobber=None):
		"""Writes a simple fits file to disk --- any acceptable type of data"
				--- Phil Hodge gave the original to me"""


		if clobber:
			rmifexist(fn=fn)

		fitsobj = pyfits.HDUList ()

		hdu = pyfits.PrimaryHDU ()
		hdu.data = data
		hdr = hdu.header

		fitsobj.append(hdu)



		if kwdict:

			if kwlist is not None:
				for kwd in kwlist:
					if kwdict.has_key(kwd):
						v = kwdict[kwd]
						if type(v) is types.TupleType:
							hdr.update(kwd, v[0], comment=v[1])
						else:
							v = kwdict[kwd]
							hdr.update(kwd, v, comment=None)

			else: # no keyword list supplied...
				for k,v in kwdict.items():
					if type(v) is types.TupleType:
						hdr.update(k, v[0], comment=v[1])
					else:
						hdr.update(k, v, comment=None)

		if verbose:
			print "  %s ---"  %  fn
			if verbose == 'vv':
				for kvpair in  hdr.items():
					print "    %8s = "  %  kvpair[0],
					print kvpair[1]

		fitsobj.writeto(fn)
		fitsobj.close()
		if verbose:
			print "WriteSimpleFits: '%s' written"  %  fn




class SimpleFitsRead:

	def __init__(self, fn=None, report=None, verbose=None):
		"""reads a simple fits file from disk, 2d """
		try:
			if verbose:
				print "Trying to open " + "'"+fn+"'"
			fimg = pyfits.open (fn)
			self.hdr = fimg[0].header
			self.data = fimg[0].data
			self.fn = fn
		except IOError:
			print ': FATAL: ReadSimpleFits: Cannot open ' + "'"+fn+"'"

		self.dimension = self.hdr['NAXIS']
		self.numCols = self.hdr['NAXIS1']
		self.numRows = self.hdr['NAXIS2']

		self.immax = self.data.max()
		self.immin = self.data.min()

		if report or verbose:
			print "Opened " + "'"+fn+"'"
		if verbose:
			for ascard in self.hdr.ascardlist():
				print ascard

		return None


def testWrite():

	# N.Float64 is an alternate data type...
	d = N.zeros((128,128), N.float32) + 42.0
	kd = {}
	kd['answer'] =  (42, "the answer")
	kd['nswer'] =  ('0 + 42', 'a string strung')
	kd['swer'] =  (42.00001, "floating by")
	kd['truth'] =  (True, "boolinary delight")
	kd['lies'] =  (False, "schweinhund")
	kd['lyes'] =  False

	SimpleFitsWrite(fn='testdata/SimpleFits.fits',    data=d, clobber='y')
	SimpleFitsWrite(fn='testdata/SimpleFitsKWD_CMT.fits', data=d, kwdict = kd, \
	                clobber='y')
	SimpleFitsWrite(fn='testdata/SimpleFitsV.fits',   data=d, verbose='v', clobber='y') 
	SimpleFitsWrite(fn='testdata/SimpleFitsVV.fits',  data=d, kwdict = kd, verbose='vv', clobber='y')


def testRead():

	print 'silent...'
	indata = SimpleFitsRead(fn="testdata/SimpleFitsKWD.fits")
	print 'report...'
	indata = SimpleFitsRead(fn="testdata/SimpleFitsKWD.fits",  \
                         	report=0xffff)
	print 'verbose...'
	indata = SimpleFitsRead(fn="testdata/SimpleFitsKWD.fits",  \
                         	verbose=0xffff)

	print 'report verbose...'
	indata = SimpleFitsRead(fn="testdata/SimpleFitsKWD.fits",  \
                         	report=0xffff,  verbose=0xffff)
	print indata.data[0,0]
	print 'Done'

def testExtRead():

	print 'silent...'
	indata = ExtFitsRead(fn="testdata/ExtFitsKWD.fits")
	print 'report...'
	indata = ExtFitsRead(fn="testdata/ExtFitsKWD.fits", extnum=1,  \
                         	report=0xffff)
	print 'verbose...'
	indata = ExtFitsRead(fn="testdata/ExtFitsKWD.fits", extnum=1,  \
                         	verbose=0xffff)

	print 'report verbose...'
	indata = ExtFitsRead(fn="testdata/ExtFitsKWD.fits", extnum=1,  \
                         	report=0xffff,  verbose=0xffff)
	print indata.data[0,0]
	print 'Done'


if __name__ == '__main__':

	def test():

		testWrite()
		#estRead()
		#estExtRead()

	test()
