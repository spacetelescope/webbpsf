import sys
import multiprocessing
import logging
import logging.handlers
import datetime
import time
from os.path import abspath, sep, join, exists, isdir, split
import os
from itertools import product, chain
import matplotlib
matplotlib.use('Agg')
if not os.environ.get('PYSYN_CDBS'):
    os.environ['PYSYN_CDBS'] = '/grp/hst/cdbs'
assert exists(os.environ['PYSYN_CDBS']), "Can't load synthetic photometry files!"

if not os.environ.get('WEBBPSF_PATH'):
    os.environ['WEBBPSF_PATH'] = '/grp/jwst/ote/webbpsf-data'
import webbpsf


N_PROCESSES = 16

def _worker_logging_setup(queue_instance):
    queue_handler = logging.handlers.QueueHandler(queue_instance)
    root = logging.getLogger()
    root.addHandler(queue_handler)
    root.setLevel(logging.DEBUG)

INSTRUMENTS = (webbpsf.NIRCam, webbpsf.NIRSpec, webbpsf.MIRI, webbpsf.NIRISS, webbpsf.FGS)

STELLAR_SPECTRAL_TYPE = 'G0V'

# NIRCam Coronagraph Ops v.3
# John Stansberry, Nov 4 2015
# Table 3. Allowed Filters for Coronagraphic Science
NIRCAM_ALLOWED_FILTERS_FOR_MASKS = {
    'MASKSWB':  ('F182M', 'F187N', 'F210M', 'F212N', 'F200W'),
    'MASK210R':  ('F182M', 'F187N', 'F210M', 'F212N', 'F200W'),
    'MASKLWB':  ('F250M', 'F300M', 'F335M', 'F360M', 'F410M', 'F430M', 'F460M', 'F480M',
                 'F277W', 'F356W', 'F444W'),
    'MASK335R': ('F300M', 'F335M', 'F360M', 'F410M', 'F356W'),
    'MASK430R': ('F410M', 'F360M', 'F430M', 'F460M', 'F444W'),
}

NIRCAM_IMAGE_MASKS_FOR_PUPILS = {
    'CIRCLYOT': ('MASK210R', 'MASK335R', 'MASK430R'),
    'WEDGELYOT': ('MASKSWB', 'MASKLWB'),
    'WEAK LENS +4': (None,),
    'WEAK LENS +8': (None,),
    'WEAK LENS -8': (None,),
    'WEAK LENS +12 (=4+8)': (None,),
    'WEAK LENS -4 (=4-8)': (None,),
}

MIRI_ALLOWED_FILTER_FOR_MASKS = {
    'FQPM1065': 'F1065C',
    'FQPM1140': 'F1140C',
    'FQPM1550': 'F1550C',
    'LYOT2300': 'F2300C',
}

MIRI_IMAGE_MASKS_FOR_PUPILS = {
    'MASKFQPM': ('FQPM1065', 'FQPM1140', 'FQPM1550'),
    'MASKLYOT': ('LYOT2300',),
    # 'P750L LRS grating': ('LRS slit',),
}

NIRSPEC_ABBREVIATED_MASK_NAMES = {
    'MSA all open': 'msa_all',
    'Single MSA open shutter': 'msa_single',
    'Three adjacent MSA open shutters': 'msa_three',
    'NIRSpec grating': 'with_grating',
}

NIRISS_LONG_FILTERS = ('F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M')

# Picked roughly the 'middle' WFE realization from the Rev. V files
# (but these can/will be updated for measured WFE maps)
OPD_TO_USE = {
    'NIRCam': 'OPD_RevW_ote_for_NIRCam_requirements.fits.gz',
    'MIRI': 'OPD_RevW_ote_for_MIRI_requirements.fits.gz',
    'NIRSpec': 'OPD_RevW_ote_for_NIRSpec_requirements.fits.gz',
    'FGS': 'OPD_RevW_ote_for_FGS_requirements.fits.gz',
    'NIRISS': 'OPD_RevW_ote_for_NIRISS_requirements.fits.gz',
}

def ensure_dir(dirpath):
    try:
        os.makedirs(dirpath)
    except OSError:
        if not isdir(dirpath):
            raise
    return dirpath

def apply_configuration(instrument_instance, configuration):
    for key, value in configuration.items():
        setattr(instrument_instance, key, value)

def make_file_path(instrument_instance, output_directory):
    output_file_path = abspath(join(output_directory, instrument_instance.name))
    parts = ['PSF', instrument_instance.name]
    if 'requirements' in instrument_instance.pupilopd:
        parts.append('requirements_opd')
    elif 'predicted' in instrument_instance.pupilopd:
        parts.append('predicted_opd')
    else:
        parts.append('perfect_opd')

    for attribute in ('filter', 'image_mask', 'pupil_mask'):
        value = getattr(instrument_instance, attribute)

        # Special case for space-filled NIRSpec mask names
        if (instrument_instance.name == 'NIRSpec' and
            value in NIRSPEC_ABBREVIATED_MASK_NAMES):
                value = NIRSPEC_ABBREVIATED_MASK_NAMES[value]

        if value is not None:
            parts.append('{}_{}'.format(attribute, value))
    return join(output_file_path, '_'.join(parts) + '.fits')

def _validate(opd, filter_name, image_mask, pupil_mask, instrument_name):
    def both_or_neither(a, b):
        if (a is None and b is None) or (a is not None and b is not None):
            return True
        else:
            return False

    if instrument_name == 'NIRCam':
        if image_mask is not None and pupil_mask is not None:
            if filter_name not in NIRCAM_ALLOWED_FILTERS_FOR_MASKS[image_mask]:
                return False
            if image_mask not in NIRCAM_IMAGE_MASKS_FOR_PUPILS[pupil_mask]:
                return False
        if not both_or_neither(pupil_mask, image_mask):
            return False
    elif instrument_name == 'MIRI':
        # Cannot simulate LRS slit in broadband mode because there's
        # no bandpass for it, and that might not make sense anyway
        if image_mask == 'LRS slit' or pupil_mask == 'P750L LRS grating':
            return False
        if image_mask is not None and pupil_mask is not None:
            if filter_name != MIRI_ALLOWED_FILTER_FOR_MASKS[image_mask]:
                return False
            if image_mask not in MIRI_IMAGE_MASKS_FOR_PUPILS[pupil_mask]:
                return False
        if not both_or_neither(pupil_mask, image_mask):
            return False
    elif instrument_name == 'NIRSpec':
        if image_mask is not None and pupil_mask != 'NIRSpec grating':
            # fixed slits get dispersed through the same gratings, so must
            # include the rectangular pupil stop from the grating
            return False
        if image_mask is None and pupil_mask is not None:
            # Even a target acq PSF should use the MSA pupil mask
            return False
        if filter_name == 'IFU':
            return False  # not yet implemented in WebbPSF
    elif instrument_name == 'FGS':
        return True
    elif instrument_name == 'NIRISS':
        if image_mask is not None:
            # Coronagraph spots are basically vestigial
            return False
        if pupil_mask in ('MASK_NRM', 'GR700XD'):
            # These should be available as cubes, but not broadband
            # PSFs. We'll wait and see if anyone requests them.
            return False
        if filter_name == 'CLEAR':
            # CLEAR would be a wide-open bandpass, which doesn't make sense
            # except for PSF cubes
            return False
        if filter_name in NIRISS_LONG_FILTERS and pupil_mask != 'CLEARP':
            # long wavelength filters cannot be configured without the modified
            # CLEARP pupil in the pupil wheel
            return False
        if filter_name not in NIRISS_LONG_FILTERS and pupil_mask == 'CLEARP':
            # likewise, short wavelength filters can only be configured with
            # the true 'CLEAR' position in the filter wheel
            return False
    else:
        return False
    return True


def generate_configuration(instrument_instance):
    selections = (
        # For all different OPDs:
        # chain([None,], instrument_instance.opd_list),
        # For only one OPD map:
        [OPD_TO_USE[instrument_instance.name],],
        instrument_instance.filter_list,
        chain([None,], instrument_instance.image_mask_list),
        chain([None,], instrument_instance.pupil_mask_list),
    )
    for opd, filter_name, image_mask, pupil_mask in product(*selections):
        if not _validate(opd, filter_name, image_mask, pupil_mask, instrument_instance.name):
            continue
        else:
            yield {
                'pupilopd': opd,
                'filter': filter_name,
                'image_mask': image_mask,
                'pupil_mask': pupil_mask,
            }

def _do_one_psf(InstrumentClass, configuration, output_directory):
    log = logging.getLogger("worker." + multiprocessing.current_process().name)
    inst = InstrumentClass()
    apply_configuration(inst, configuration)
    log.debug("Computing PSF for {} ({})".format(inst.name, configuration))
    output_file_path = make_file_path(inst, output_directory)
    log.debug("Checking for existing PSF: {}".format(output_file_path))
    if not exists(output_file_path):
        dirs, name = split(output_file_path)
        ensure_dir(dirs)
        # TODO contemplate better way to handle special case
        fov_arcsec = 20.0 if inst.name == 'MIRI' else 10.0
        spectrum = webbpsf.specFromSpectralType(STELLAR_SPECTRAL_TYPE, catalog='ck04')
        psf = inst.calc_psf(fov_arcsec=fov_arcsec, source=spectrum)
        psf.writeto(output_file_path)
        log.debug("Computed PSF\n\t{}\nand wrote to: {}".format(configuration, output_file_path))
    else:
        log.debug("Got existing PSF at {}".format(output_file_path))
    return output_file_path

def compute_library(output_directory, pool, instrument_classes=INSTRUMENTS):
    def submit_job(InstrumentClass, configuration):
        return pool.apply_async(
            _do_one_psf,
            (InstrumentClass, configuration, output_directory)
        )

    results = []

    for InstrumentClass in instrument_classes:
        inst = InstrumentClass()  # the attributes used to generate configs are set at instantiation
        for configuration in generate_configuration(inst):
            result = submit_job(InstrumentClass, configuration)
            results.append(result)

    for result in results:
        logging.info(result.get())

    return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute a PSF library for JWST")
    parser.add_argument("--nircam", action='store_true')
    parser.add_argument("--niriss", action='store_true')
    parser.add_argument("--nirspec", action='store_true')
    parser.add_argument("--fgs", action='store_true')
    parser.add_argument("--miri", action='store_true')
    parser.add_argument("--all", action='store_true')
    parser.add_argument("output_dir")
    args = parser.parse_args()
    instrument_classes = []
    if args.all:
        instrument_classes = INSTRUMENTS
    else:
        if args.nircam:
            instrument_classes.append(webbpsf.NIRCam)
        if args.niriss:
            instrument_classes.append(webbpsf.NIRISS)
        if args.nirspec:
            instrument_classes.append(webbpsf.NIRSpec)
        if args.fgs:
            instrument_classes.append(webbpsf.FGS)
        if args.miri:
            instrument_classes.append(webbpsf.MIRI)
        if len(instrument_classes) == 0:
            print("You must specify at least one instrument to compute PSFs for!")
            sys.exit(1)

    multiprocessing.set_start_method('forkserver')

    q = multiprocessing.Queue()
    stdout_handler = logging.StreamHandler(stream=sys.stdout)

    y, m, d, hh, mm, ss = datetime.datetime.now().utctimetuple()[:6]
    logging_filename = 'compute_psf_library_{y}-{m}-{d}-{hh}-{mm}-{ss}.log'.format(y=y, m=m, d=d, hh=hh, mm=mm, ss=ss)
    log = logging.getLogger()
    file_handler = logging.FileHandler(logging_filename)
    log.addHandler(file_handler)
    listener = logging.handlers.QueueListener(q, file_handler, stdout_handler, respect_handler_level=True)
    listener.start()

    pool = multiprocessing.Pool(
        processes=N_PROCESSES,
        initializer=_worker_logging_setup,
        initargs=(q,),
    )
    with pool:
        compute_library(args.output_dir, pool, instrument_classes)
    sys.exit()
