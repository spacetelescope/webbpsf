import os
import pysynphot
import webbpsf
import atpy


WebbPSF_basepath = os.getenv('WEBBPSF_PATH', default= os.path.dirname(os.path.dirname(os.path.abspath(webbpsf.__file__))) +os.sep+"data" )


def norm_one_filter(instrument, filter_, clobber=False):
    try:
        bp = pysynphot.ObsBandpass('%s,im,%s' %(instrument.lower(), filter_.lower()))
        #normalized_throughput = bp.throughput / bp.throughput.max() # set max to 1.0
        bp = bp/bp.throughput.max() # set max to 1

        t = atpy.Table()
        t.add_column('WAVELENGTH', bp.wave, unit=bp.waveunits.name)
        t.add_column('THROUGHPUT', bp.throughput)
        t.add_keyword('TELESCOP', 'JWST')
        t.add_keyword('INSTRUME', instrument)
        t.add_keyword('FILTER', filter_)
        t.add_keyword('SOURCE', 'pysynphot, normalized to peak=1')


        t.write("%s/%s/filters/%s_throughput.fits" % (WebbPSF_basepath, instrument, filter_))
        print "Wrote throughput.fits for %s %s" % (instrument, filter_)
    except:
        print "Error for %s %s" % (instrument, filter_)



def norm_all():
    for inst_name in ['MIRI','NIRCam','NIRSpec']:
        inst = webbpsf.Instrument(inst_name)
        filts = [f for f in inst.filter_list if f[0] == 'F']
        for f in filts:
            #print inst_name, f
            norm_one_filter(inst_name,f)



