# This file contains constants and data that do not fit in the data
# package for one reason or another. It could be that they are in a
# data structure that doesn't serialize well or they're too small to
# make a separate file.

# It's hard to make anything truly immutable in Python, but here
# tuples are preferred and numpy arrays should have
# flags.writeable = False

import numpy as np

__all__ = (
    'JWST_PRIMARY_SEGMENTS',
    'JWST_PRIMARY_STRUTS',
    'JWST_PRIMARY_SEGMENT_CENTERS',
    'JWST_SEGMENT_RADIUS',
    'JWST_CIRCUMSCRIBED_DIAMETER',
    'SEGNAMES',
    'SEGNAMES_WSS'
)

SEGNAMES = tuple([letter + str(number) for letter in ['A', 'B', 'C'] for number in range(1, 7)])

#
# JWST Primary segment and obscuration shapes and centers
#
# Provenance:
#
# Segment coordinates from "2010.03.16 Transmission X Area Budget.xls"
# spreadsheet by Paul Lightsey, which was based in turn on
# Ball Aerospace drawing # 2220169 Rev B and the OTE Cryogenic Optics
# Interface Control Document, Ball Aerospace doc # C327693

JWST_PRIMARY_SEGMENTS = (
    ('A1-1', np.array([
        [-0.38101, 0.667604],
        [-0.758826, 1.321999],
        [-0.38101, 1.976407],
        [0.38101, 1.976407],
        [0.758826, 1.321999],
        [0.38101, 0.667604]])),
    ('A2-2', np.array([
        [0.38765702, 0.66376634],
        [0.76547172, 1.31816209],
        [1.52111367, 1.31816784],
        [1.90212367, 0.65823916],
        [1.52429772, 0.00383691],
        [0.76866702, 0.00383766]])),
    ('A3-3', np.array([
        [0.76866702, -0.00383766],
        [1.52429772, -0.00383691],
        [1.90212367, -0.65823916],
        [1.52111367, -1.31816784],
        [0.76547172, -1.31816209],
        [0.38765702, -0.66376634]])),
    ('A4-4', np.array([
        [0.38101, -0.667604],
        [0.758826, -1.321999],
        [0.38101, -1.976407],
        [-0.38101, -1.976407],
        [-0.758826, -1.321999],
        [-0.38101, -0.667604]])),
    ('A5-5', np.array([
        [-0.38765702, -0.66376634],
        [-0.76547172, -1.31816209],
        [-1.52111367, -1.31816784],
        [-1.90212367, -0.65823916],
        [-1.52429772, -0.00383691],
        [-0.76866702, -0.00383766]])),
    ('A6-6', np.array([
        [-0.76866702, 0.00383766],
        [-1.52429772, 0.00383691],
        [-1.90212367, 0.65823916],
        [-1.52111367, 1.31816784],
        [-0.76547172, 1.31816209],
        [-0.38765702, 0.66376634]])),
    ('B1-7', np.array([
        [0.38101, 3.279674],
        [0.758826, 2.631791],
        [0.38101, 1.98402],
        [-0.38101, 1.98402],
        [-0.758826, 2.631791],
        [-0.38101, 3.279674]])),
    ('B2-9', np.array([
        [3.030786, 1.30987266],
        [2.65861086, 0.65873291],
        [1.90871672, 0.66204566],
        [1.52770672, 1.32197434],
        [1.89978486, 1.97305809],
        [2.649776, 1.96980134]])),
    ('B3-11', np.array([
        [2.649776, -1.96980134],
        [1.89978486, -1.97305809],
        [1.52770672, -1.32197434],
        [1.90871672, -0.66204566],
        [2.65861086, -0.65873291],
        [3.030786, -1.30987266]])),
    ('B4-13', np.array([
        [-0.38101, -3.279674],
        [-0.758826, -2.631791],
        [-0.38101, -1.98402],
        [0.38101, -1.98402],
        [0.758826, -2.631791],
        [0.38101, -3.279674]])),
    ('B5-15', np.array([
        [-3.030786, -1.30987266],
        [-2.65861086, -0.65873291],
        [-1.90871672, -0.66204566],
        [-1.52770672, -1.32197434],
        [-1.89978486, -1.97305809],
        [-2.649776, -1.96980134]])),
    ('B6-17', np.array([
        [-2.649776, 1.96980134],
        [-1.89978486, 1.97305809],
        [-1.52770672, 1.32197434],
        [-1.90871672, 0.66204566],
        [-2.65861086, 0.65873291],
        [-3.030786, 1.30987266]])),
    ('C1-8', np.array([
        [0.765201, 2.627516],
        [1.517956, 2.629178],
        [1.892896, 1.976441],
        [1.521076, 1.325812],
        [0.765454, 1.325807],
        [0.387649, 1.980196]])),
    ('C2-10', np.array([
        [2.6580961, 0.651074495],
        [3.03591294, 5.42172989e-07],
        [2.65809612, -0.651075523],
        [1.90872487, -0.654384457],
        [1.53090954, 8.90571587e-07],
        [1.90872454, .654384118]])),
    ('C3-12', np.array([
        [1.8928951, -1.97644151],
        [1.51795694, -2.62917746],
        [0.76520012, -2.62751652],
        [0.38764887, -1.98019646],
        [0.76545554, -1.32580611],
        [1.52107554, -1.32581188]])),
    ('C4-14', np.array([
        [-0.765201, -2.627516],
        [-1.517956, -2.629178],
        [-1.892896, -1.976441],
        [-1.521076, -1.325812],
        [-0.765454, -1.325807],
        [-0.387649, -1.980196]])),
    ('C5-16', np.array([
        [-2.6580961, -.651074495],
        [-3.03591294, -5.42172990e-07],
        [-2.65809612, .651075523],
        [-1.90872487, .654384457],
        [-1.53090954, -8.90571587e-07],
        [-1.90872454, -.654384118]])),
    ('C6-18', np.array([
        [-1.8928951, 1.97644151],
        [-1.51795694, 2.62917746],
        [-0.76520012, 2.62751652],
        [-0.38764887, 1.98019646],
        [-0.76545554, 1.32580611],
        [-1.52107554, 1.32581188]])),
)

for name, arr in JWST_PRIMARY_SEGMENTS:
    arr.flags.writeable = False

# A1-6,B1-6,C1-6
SEGNAMES_WSS = tuple(name for name, arr in JWST_PRIMARY_SEGMENTS)

# Sort same names by another order: A1-6,B1,C1,B2,C2,etc
SEGNAMES_WSS_ORDER = tuple(np.asarray(SEGNAMES_WSS)[
                               np.argsort([int(a.split('-')[1]) for a in SEGNAMES_WSS])])

JWST_PRIMARY_STRUTS = (
    ("strut1", np.array([
        [-0.05301375, -0.0306075],
        [1.59698625, -2.88849133],
        [1.70301375, -2.82727633],
        [0.05301375, 0.0306075],
        [-0.05301375, -0.0306075]])),
    ("strut2", np.array([
        [-0.05301375, 0.0306075],
        [-1.70301375, -2.82727633],
        [-1.59698625, -2.88849133],
        [0.05301375, -0.0306075],
        [-0.05301375, 0.0306075]])),
    ("strut3", np.array([
        [5.94350000e-02, -1.45573765e-17],
        [5.94350000e-02, 3.30000000e+00],
        [-5.94350000e-02, 3.30000000e+00],
        [-5.94350000e-02, 1.45573765e-17],
        [5.94350000e-02, -1.45573765e-17]])),
    ("strut3_bumps", np.array([
        [0.059435, 0.666],
        [0.059435, 2.14627],
        [0.082595, 2.14627],
        [0.082595, 2.3645],
        [0.059435, 2.3645],
        [0.059435, 2.48335],
        [0.069795, 2.48335],
        [0.069795, 2.54445],
        [0.059435, 2.54445],
        [0.059435, 3.279674],
        [-0.059435, 3.279674],
        [-0.059435, 2.54445],
        [-0.069795, 2.54445],
        [-0.069795, 2.48335],
        [-0.059435, 2.48335],
        [-0.059435, 2.3645],
        [-0.082595, 2.3645],
        [-0.082595, 2.14627],
        [-0.059435, 2.14627],
        [-0.059435, 0.666]]))
)

for name, arr in JWST_PRIMARY_STRUTS:
    arr.flags.writeable = False

JWST_PRIMARY_SEGMENT_CENTERS = (
    ('A1-1', (0.000000, 1.323500)),
    ('A2-2', (1.146185, 0.661750)),
    ('A3-3', (1.146185, -0.661750)),
    ('A4-4', (0.000000, -1.323500)),
    ('A5-5', (-1.146185, -0.661750)),
    ('A6-6', (-1.146185, 0.661750)),
    ('B1-7', (0.000000, 2.634719)),
    ('B2-9', (2.281734, 1.317360)),
    ('B3-11', (2.281734, -1.317359)),
    ('B4-13', (0.000000, -2.634719)),
    ('B5-15', (-2.281734, -1.317360)),
    ('B6-17', (-2.281734, 1.317360)),
    ('C1-8', (1.142963, 1.979670)),
    ('C2-10', (2.285926, 0.000000)),
    ('C3-12', (1.142963, -1.979670)),
    ('C4-14', (-1.142963, -1.979670)),
    ('C5-16', (-2.285926, -0.000000)),
    ('C6-18', (-1.142963, 1.979670))
)

# TODO - add in V1 positions to the above? 0.055154 for As, 0.218578 Bs, 0.164535 Cs

JWST_SEGMENT_RADIUS = 1.517 / 2
JWST_CIRCUMSCRIBED_DIAMETER = 6.603464  # meters. Outer corners of B segments
JWST_INSCRIBED_DIAMETER = 5.47334  # meters. Middle corners of C segments

JWST_TYPICAL_LOS_JITTER_PER_AXIS = 0.0008 # milliarcseconds jitter, 1 sigma per axis. = approx 1 mas rms radial, typically

# Alignment information about instrument internal pupil masks (
INSTRUMENT_PUPIL_MASK_DEFAULT_POSITIONS = {
    'NIRCam_MASKSWB': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
    'NIRCam_MASKLWB': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
    'NIRCam_MASKRND_SW': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
    'NIRCam_MASKRND_LW': {'pupil_shift_x': -0.012, 'pupil_shift_y': -0.023, 'pupil_rotation': -0.60},  # from K. Lawson, fits to ERS progid 1386 data
    'MIRI_MASKFQPM_F1065C': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
    'MIRI_MASKFQPM_F11140': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
    'MIRI_MASKFQPM_F1550C': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
    'MIRI_MASKLYOT': {'pupil_shift_x': None, 'pupil_shift_y': None, 'pupil_rotation': None},
}

# ad hoc, highly simplified models for charge diffusion within detectors
# These values are PLACEHOLDERS and should be updated based on comparisons with data and ePSFs (ongoing)
# Note, these are parameterized as arcseconds for convenience (and consistency with the jitter paramater)
# but the underlying physics cares more about detector pixel pitch.
INSTRUMENT_DETECTOR_CHARGE_DIFFUSION_DEFAULT_PARAMETERS = {
    'NIRCAM_SW': 0.0062,    # Fit by Marcio to WFS TA ePSFs, and by Marshall to prelim NIRCam SW ePSFs by J. Anderson
    'NIRCAM_LW': 0.018,     # Fit by Marshall to prelim LW ePSFs by J. Anderson
    'NIRISS': 0.0202,       # Fit by Marcio to MIMF-3 F158M (ePSF), and by Marshall to NIRISS ePSFs by Anderson & Libralato
    'FGS': 0.07,            # Fit by Marcio to FGS_ID images
    'NIRSPEC': 0.036,
    'MIRI': 0.001,          # Fit by Marshall + Marcio to ePSFs, after adding IPC
                            #  0.070 Based on user reports, see issue #674. However, this was before adding IPC effects
}
# add Interpixel capacitance (IPC) effects. These are the parameters for each detector kernel
# For NIRCam we  use CV3/Flight convolution kernels from Jarron Leisenring, see detectors.apply_detector_ipc for details
# NIRISS has different kernels provided by Kevin Volk (STScI), see detectors.apply_detector_ipc for details
INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS = {
    'MIRI': (0.033, 0.024, 0.013),          # Based on JWST-STScI-002925 by Mike Engesser
}

# How many detector pixels to mask out for the inner "hole" in the cruciform?
# See Gaspar et al. 2021 for illustrative figures.
# This is a rough approximation of a detector-position-dependent phenomenon
MIRI_CRUCIFORM_INNER_RADIUS_PIX = 12
