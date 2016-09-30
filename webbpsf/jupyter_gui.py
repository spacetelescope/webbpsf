from __future__ import division, print_function, absolute_import, unicode_literals
import logging
import matplotlib
import astropy.io.fits as fits
import six

import poppy
from . import webbpsf_core
from . import wfirst

_log = logging.getLogger('webbpsf')



def show_notebook_interface(instrumentname):
    """
    Show Jupyter notebook widget interface

    Parameters
    -----------
    instrumentname : string
        one of 'NIRCam','NIRSpec','NIRISS','FGS','MIRI'
        or 'WFI'
    """

    if instrumentname.upper()=='WFI':
        instrument = wfirst.WFI()
        show_notebook_interface_wfi(instrument)
    else:
        instrument = webbpsf_core.Instrument(instrumentname)
        show_notebook_interface_jwst(instrument)



def show_notebook_interface_wfi(instrument):
    # Widget related imports.
    # (Currently not a hard dependency for the full webbpsf package, so we import
    # within the function.)
    import ipywidgets as widgets
    from IPython.display import display, clear_output
    from matplotlib import pyplot as plt

    try:
        import pysynphot
    except ImportError:
        raise ImportError("For now, PySynphot must be installed to use the notebook interface")

    # Clean up some warnings we know about so as not to scare the users
    import warnings
    from matplotlib.cbook import MatplotlibDeprecationWarning
    warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
    warnings.simplefilter('ignore', fits.verify.VerifyWarning)

    def make_binding_for_attribute(attribute):
        def callback(trait_name, new_value):
            setattr(instrument, attribute, new_value)
        return callback

    filter_selection = widgets.ToggleButtons(
        options=instrument.filter_list,
        value=instrument.filter,
        description='Filter:'
    )
    filter_selection.on_trait_change(
        make_binding_for_attribute('filter'),
        name='selected_label'
    )
    display(filter_selection)

    monochromatic_wavelength = widgets.BoundedFloatText(
        value=0.76,
        min=0.6,
        max=2.0,
    )
    monochromatic_wavelength.disabled = True
    monochromatic_toggle = widgets.Checkbox(description='Monochromatic calculation?')

    def update_monochromatic(trait_name, new_value):
        filter_selection.disabled = new_value
        monochromatic_wavelength.disabled = not new_value

    monochromatic_toggle.on_trait_change(update_monochromatic, name='value')

    display(widgets.HTML(value='''<p style="padding: 1em 0;">
    <span style="font-style:italic; font-size:1.0em">
    Monochromatic calculations can be performed for any wavelength in the 0.6 to 2.0 &micro;m range.
    </span></p>'''))  # kludge
    monochromatic_controls = widgets.HBox(children=(
            monochromatic_toggle,
            widgets.HTML(value='<span style="display: inline-block; width: 0.6em;"></span>'),
            monochromatic_wavelength,
            widgets.HTML(value='<span style="display: inline-block; width: 0.25em;"></span> &micro;m '),
    ))
    display(monochromatic_controls)

    display(widgets.HTML(value="<hr>"))

    source_selection = widgets.Select(
        options=poppy.specFromSpectralType('', return_list=True),
        value='G0V',
        description="Source spectrum"
    )
    display(source_selection)
    display(widgets.HTML(value="<hr>"))

    sca_selection = widgets.Dropdown(
        options=instrument.detector_list,
        value=instrument.detector,
        description='Detector:'
    )
    sca_selection.on_trait_change(
        make_binding_for_attribute('detector'),
        name='selected_label'
    )
    display(sca_selection)

    detector_field_points = [
        ('Top left', (4.0, 4092.0)),
        ('Bottom left', (4.0, 4.0)),
        ('Center', (2048.0, 2048.0)),
        ('Top right', (4092.0, 4092.0)),
        ('Bottom right', (4092.0, 4.0)),
    ]
    # enforce ordering of buttons
    detector_field_point_labels = [a[0] for a in detector_field_points]
    detector_field_points = dict(detector_field_points)

    def set_field_position(trait_name, new_value):
        instrument.detector_position = detector_field_points[new_value]

    field_position = widgets.ToggleButtons(options=detector_field_point_labels, value='Center', description='Detector field point:')
    field_position.on_trait_change(set_field_position, name='selected_label')
    display(field_position)

    calculate_button = widgets.Button(
        description="Calculate PSF",
        width='10em',
        color='white',
        background_color='#00c403',
        border_color='#318732'
    )
    display_osys_button = widgets.Button(
        description="Display Optical System",
        width='13em',
        color='white',
        background_color='#005fc4',
        border_color='#224A75'
    )
    clear_button = widgets.Button(
        description="Clear Output",
        width='10em',
        color='white',
        background_color='#ed4747',
        border_color='#911C1C'
    )
    progress = widgets.HTML(value='<progress>')

    OUTPUT_FILENAME = 'psf.fits'
    DOWNLOAD_BUTTON_HTML = """
    <a class="btn btn-info" href="files/{}" target="_blank">
        Download FITS image from last calculation
    </a>
    """
    download_link = widgets.HTML(value=DOWNLOAD_BUTTON_HTML.format(OUTPUT_FILENAME))

    def disp(*args):
        progress.visible = True
        plt.figure(figsize=(12, 8))
        instrument.display()
        progress.visible = None

    def calc(*args):
        progress.visible = True
        if monochromatic_toggle.value is True:
            psf = instrument.calcPSF(
                monochromatic=monochromatic_wavelength.value * 1e-6,
                display=True,
                outfile=OUTPUT_FILENAME,
                clobber=True
            )
        else:
            source = poppy.specFromSpectralType(source_selection.value)
            _log.debug("Got source type {}: {}".format(source_selection.value, source))
            psf = instrument.calcPSF(
                source=source,
                display=True,
                outfile=OUTPUT_FILENAME,
                clobber=True
            )
        fig, (ax_oversamp, ax_detsamp) = plt.subplots(1, 2)
        poppy.display_PSF(psf, ax=ax_oversamp)
        poppy.display_PSF(psf, ax=ax_detsamp, ext='DET_SAMP')
        progress.visible = None
        download_link.visible = True

    def clear(*args):
        clear_output()
        progress.visible = None
        download_link.visible = None

    calculate_button.on_click(calc)
    display_osys_button.on_click(disp)
    clear_button.on_click(clear)
    display(widgets.HTML(value="<br/>"))  # kludge
    buttons = widgets.HBox(children=[calculate_button, display_osys_button, clear_button])
    display(buttons)

    # Insert the progress bar, hidden by default
    display(progress)
    progress.visible = None
    # and the download link
    display(download_link)
    download_link.visible = None

def show_notebook_interface_jwst(instrument):
    """ Show Jupyter Notebook interface, for a JWST instrument

    Parameters
    -------------
    instrument : string or object
        either a webbpsf instrument object, e.g.  `NIRCam()`
        or the string name of an instrument.


    Example
    --------
    nc = webbpsf.NIRCam()
    webbpsf.show_notebook_interface(nc)
    """
    # Widget related imports.
    # (Currently not a hard dependency for the full webbpsf package, so we import
    # within the function.)
    import ipywidgets as widgets
    from IPython.display import display, clear_output
    from matplotlib import pyplot as plt


    if isinstance(instrument, six.string_types):
        instrument = Instrument(instrument)

    try:
        import pysynphot
    except ImportError:
        raise ImportError("For now, PySynphot must be installed to use the notebook interface")

    # Clean up some warnings we know about so as not to scare the users
    import warnings
    from matplotlib.cbook import MatplotlibDeprecationWarning
    warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
    warnings.simplefilter('ignore', fits.verify.VerifyWarning)

    def make_binding_for_attribute(attribute):
        def callback(trait_name, new_value):
            if new_value == 'None':
                setattr(instrument, attribute, None)
            else:
                setattr(instrument, attribute, new_value)
        return callback

    display(widgets.HTML(value='''<p style="padding: 1em 0;">
    <span style="font-weight:bold; font-size:1.0em">
    Notebook Interface for {} PSF sims
    </span></p>'''.format(instrument.name)))


    filter_selection = widgets.Dropdown(
        options=instrument.filter_list,
        value=instrument.filter,
        description='Filter:')

    filter_selection.on_trait_change(
        make_binding_for_attribute('filter'),
        name='selected_label'
    )
    display(filter_selection)


    wl_bounds = (5., 30., 10.0) if instrument.name=='MIRI' else (0.6, 5.3, 2.0)
    monochromatic_wavelength = widgets.BoundedFloatText(
        value=wl_bounds[2],
        min=wl_bounds[0],
        max=wl_bounds[1],
    )
    monochromatic_wavelength.disabled = True
    monochromatic_toggle = widgets.Checkbox(description='Monochromatic calculation?')

    def update_monochromatic(trait_name, new_value):
        filter_selection.disabled = new_value
        monochromatic_wavelength.disabled = not new_value

    monochromatic_toggle.on_trait_change(update_monochromatic, name='value')

    display(widgets.HTML(value='''<p style="padding: 1em 0;">
    <span style="font-style:italic; font-size:1.0em">
    Monochromatic calculations can be performed for any wavelength in the {} to {} &micro;m range.
    </span></p>'''.format(*wl_bounds)))  # kludge
    monochromatic_controls = widgets.HBox(children=(
            monochromatic_toggle,
            widgets.HTML(value='<span style="display: inline-block; width: 0.6em;"></span>'),
            monochromatic_wavelength,
            widgets.HTML(value='<span style="display: inline-block; width: 0.25em;"></span> &micro;m '),
    ))
    display(monochromatic_controls)
    display(widgets.HTML(value="<hr>"))


    if instrument.name != 'FGS':
        image_selection = widgets.Dropdown(
            options=['None'] + instrument.image_mask_list,
            value=str(instrument.image_mask),
            description='Image Mask:')

        image_selection.on_trait_change(
            make_binding_for_attribute('image_mask'),
            name='selected_label'
        )
        display(image_selection)


        pupil_selection = widgets.Dropdown(
            options=['None'] + instrument.pupil_mask_list,
            value=str(instrument.pupil_mask),
            description='Pupil Mask: ')

        pupil_selection.on_trait_change(
            make_binding_for_attribute('pupil_mask'),
            name='selected_label'
        )
        display(pupil_selection)


    display(widgets.HTML(value="<hr>"))

    source_selection = widgets.Dropdown(
        options=poppy.specFromSpectralType('', return_list=True),
        value='G0V',
        description="Source spectrum"
    )
    display(source_selection)
    display(widgets.HTML(value="<hr>"))

    calculate_button = widgets.Button(
        description="Calculate PSF",
        width='10em',
        color='white',
        background_color='#00c403',
        border_color='#318732'
    )
    display_osys_button = widgets.Button(
        description="Display Optical System",
        width='13em',
        color='white',
        background_color='#005fc4',
        border_color='#224A75'
    )
    clear_button = widgets.Button(
        description="Clear Output",
        width='10em',
        color='white',
        background_color='#ed4747',
        border_color='#911C1C'
    )
    progress = widgets.HTML(value='<progress>')

    OUTPUT_FILENAME = 'psf.fits'
    DOWNLOAD_BUTTON_HTML = """
    <a class="btn btn-info" href="files/{}" target="_blank">
        Download FITS image from last calculation
    </a>
    """
    download_link = widgets.HTML(value=DOWNLOAD_BUTTON_HTML.format(OUTPUT_FILENAME))

    def disp(*args):
        progress.visible = True
        plt.figure(figsize=(12, 8))
        instrument.display()
        progress.visible = None

    def calc(*args):
        progress.visible = True
        if monochromatic_toggle.value is True:
            psf = instrument.calcPSF(
                monochromatic=monochromatic_wavelength.value * 1e-6,
                display=True,
                outfile=OUTPUT_FILENAME,
                clobber=True
            )
        else:
            source = poppy.specFromSpectralType(source_selection.value)
            _log.debug("Got source type {}: {}".format(source_selection.value, source))
            psf = instrument.calcPSF(
                source=source,
                display=True,
                outfile=OUTPUT_FILENAME,
                clobber=True
            )
        fig, (ax_oversamp, ax_detsamp) = plt.subplots(1, 2,figsize=(12, 4))
        title1 = "PSF sim for {}, {}\n".format(instrument.name, instrument.filter)
        poppy.display_PSF(psf, ax=ax_oversamp,
                          title=title1+"Oversampled PSF")
        poppy.display_PSF(psf, ax=ax_detsamp, ext='DET_SAMP',
                          title=title1+'Detector pixel sampled PSF')
        progress.visible = None
        download_link.visible = True

    def clear(*args):
        clear_output()
        progress.visible = None
        download_link.visible = None

    calculate_button.on_click(calc)
    display_osys_button.on_click(disp)
    clear_button.on_click(clear)
    display(widgets.HTML(value="<br/>"))  # kludge
    buttons = widgets.HBox(children=[calculate_button, display_osys_button, clear_button])
    display(buttons)

    # Insert the progress bar, hidden by default
    display(progress)
    progress.visible = None
    # and the download link
    display(download_link)
    download_link.visible = None
