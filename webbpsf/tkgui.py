#!/usr/bin/env python
from __future__ import division, print_function, absolute_import, unicode_literals
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.io.fits as fits
from six.moves import tkinter as tk
from six.moves import tkinter_tkfiledialog as tkFileDialog
from six.moves import tkinter_messagebox as tkMessageBox
import logging
import logging
_log = logging.getLogger('webbpsf')


try:
    from six.moves import tkinter_ttk as ttk
except ImportError:
    raise RuntimeError("Python 2.7.0 (or newer) with ttk widget support is required")



try:
    import pysynphot
    path = os.getenv('PYSYN_CDBS')
    if path is not None and os.path.isdir(path):
        _HAS_PYSYNPHOT=True
    else:
        _HAS_PYSYNPHOT=False
except ImportError:
    _HAS_PYSYNPHOT=False


import poppy
from . import webbpsf_core

class WebbPSF_GUI(object):
    """ A GUI for the PSF Simulator

    Documentation TBD!

    """
    def __init__(self, opdserver=None):
        # init the object and subobjects
        self.instrument = {}
        self.widgets = {}
        self.vars = {}
        self.advanced_options = {'parity': 'any', 'force_coron': False, 'no_sam': False, 'psf_vmin': 1e-8, 'psf_vmax': 1.0, 'psf_scale': 'log', 'psf_cmap_str': 'Jet (blue to red)' , 'psf_normalize': 'Total', 'psf_cmap': matplotlib.cm.jet}
        insts = ['NIRCam', 'NIRSpec','NIRISS', 'MIRI', 'FGS']
        for i in insts:
            self.instrument[i] = webbpsf_core.Instrument(i)
        #invoke link to ITM server if provided?
        if opdserver is not None:
            self._enable_opdserver = True
            self._opdserver = opdserver
        else:
            self._enable_opdserver = False

        
        # create widgets & run
        self._create_widgets()
        self.root.update()


    def _add_labeled_dropdown(self, name, root,label="Entry:", values=None, default=None, width=5, position=(0,0), **kwargs):
        "convenient wrapper for adding a Combobox"

        ttk.Label(root, text=label).grid(row=position[0],  column=position[1], sticky='W')

        self.vars[name] = tk.StringVar()
        self.widgets[name] = ttk.Combobox(root, textvariable=self.vars[name], width=width, state='readonly')
        self.widgets[name].grid(row=position[0], column=position[1]+1, **kwargs)
        self.widgets[name]['values'] = values

        if default is None: default=values[0]
        self.widgets[name].set(default)
 

    def _add_labeled_entry(self, name, root,label="Entry:", value="", width=5, position=(0,0), postlabel=None, **kwargs):
        "convenient wrapper for adding an Entry"
        ttk.Label(root, text=label).grid(row=position[0],  column=position[1], sticky='W')

        self.vars[name] = tk.StringVar()
        self.widgets[name] = ttk.Entry(root, textvariable=self.vars[name], width=width)
        self.widgets[name].insert(0,value)
        self.widgets[name].grid(row=position[0], column=position[1]+1, **kwargs)

        if postlabel is not None:
            ttk.Label(root, text=postlabel).grid(row=position[0],  column=position[1]+2, sticky='W')


    def _create_widgets(self):
        """Create a nice GUI using the enhanced widget set provided by 
        the ttk extension to Tkinter, available in Python 2.7 or newer
        """
        #---- create the GUIs
        insts = ['NIRCam', 'NIRSpec','NIRISS', 'MIRI',  'FGS']
        self.root = tk.Tk()
        self.root.geometry('+50+50')
        self.root.title("James Webb Space Telescope PSF Calculator")

        frame = ttk.Frame(self.root)
        #frame = ttk.Frame(self.root, padx=10,pady=10)

        #ttk.Label(frame, text='James Webb PSF Calculator' ).grid(row=0)

        #-- star
        lf = ttk.LabelFrame(frame, text='Source Properties')

        if _HAS_PYSYNPHOT:
            self._add_labeled_dropdown("SpType", lf, label='    Spectral Type:', values=poppy.specFromSpectralType("",return_list=True), default='G0V', width=25, position=(0,0), sticky='W')
            ttk.Button(lf, text='Plot spectrum', command=self.ev_plotspectrum).grid(row=0,column=2,sticky='E',columnspan=4)

        r = 1
        fr2 = ttk.Frame(lf)

        self._add_labeled_entry("source_off_r", fr2, label='    Source Position: r=', value='0.0', width=5, position=(r,0), sticky='W')
        self._add_labeled_entry("source_off_theta", fr2, label='arcsec,  PA=', value='0', width=3, position=(r,2), sticky='W')

        self.vars["source_off_centerpos"] = tk.StringVar()
        self.vars["source_off_centerpos"].set('corner')

        ttk.Label(fr2, text='deg, centered on ' ).grid(row=r, column=4)
        pixel = ttk.Radiobutton(fr2, text='pixel', variable=self.vars["source_off_centerpos"], value='pixel')
        pixel.grid(row=r, column=5)
        corner = ttk.Radiobutton(fr2, text='corner', variable=self.vars["source_off_centerpos"], value='corner')
        corner.grid(row=r, column=6)
        fr2.grid(row=r, column=0, columnspan=5, sticky='W')



        lf.columnconfigure(2, weight=1)
        lf.grid(row=1, sticky='E,W', padx=10,pady=5)

        #-- instruments
        lf = ttk.LabelFrame(frame, text='Instrument Config')
        notebook = ttk.Notebook(lf)
        self.widgets['tabset'] = notebook
        notebook.pack(fill='both')
        for iname,i in zip(insts, range(len(insts))):
            page = ttk.Frame(notebook)
            notebook.add(page,text=iname) 
            notebook.select(i)  # make it active
            self.widgets[notebook.select()] = iname # save reverse lookup from meaningless widget "name" to string name
            if iname =='NIRCam':
                lframe = ttk.Frame(page)

                ttk.Label(lframe, text='Configuration Options for '+iname+',     module: ').grid(row=0, column=0, sticky='W')
                mname='NIRCam module'
                self.vars[mname] = tk.StringVar()
                self.widgets[mname] = ttk.Combobox(lframe, textvariable=self.vars[mname], width=2, state='readonly')
                self.widgets[mname].grid(row=0,column=1, sticky='W')
                self.widgets[mname]['values'] = ['A','B']
                self.widgets[mname].set('A')

                lframe.grid(row=0, columnspan=2, sticky='W')
            else:
                ttk.Label(page, text='Configuration Options for '+iname+"                      ").grid(row=0, columnspan=2, sticky='W')

            ttk.Button(page, text='Display Optics', command=self.ev_displayOptics ).grid(column=2, row=0, sticky='E', columnspan=3)


            #if  iname != 'TFI':
            self._add_labeled_dropdown(iname+"_filter", page, label='    Filter:', values=self.instrument[iname].filter_list, default=self.instrument[iname].filter, width=12, position=(1,0), sticky='W')
            #else:
                #ttk.Label(page, text='Etalon wavelength: ' , state='disabled').grid(row=1, column=0, sticky='W')
                #self.widgets[iname+"_wavelen"] = ttk.Entry(page, width=7) #, disabledforeground="#A0A0A0")
                #self.widgets[iname+"_wavelen"].insert(0, str(self.instrument[iname].etalon_wavelength))
                #self.widgets[iname+"_wavelen"].grid(row=1, column=1, sticky='W')
                #ttk.Label(page, text=' um' ).grid(row=1, column=2, sticky='W')
 
            #self.vars[iname+"_filter"] = tk.StringVar()
            #self.widgets[iname+"_filter"] = ttk.Combobox(page,textvariable =self.vars[iname+"_filter"], width=10, state='readonly')
            #self.widgets[iname+"_filter"]['values'] = self.instrument[iname].filter_list
            #self.widgets[iname+"_filter"].set(self.instrument[iname].filter)
            #self.widgets[iname+"_filter"]['readonly'] = True
            #ttk.Label(page, text='    Filter: ' ).grid(row=1, column=0)
            #self.widgets[iname+"_filter"].grid(row=1, column=1)


            #if hasattr(self.instrument[iname], 'ifu_wavelength'):
            if iname == 'NIRSpec' or iname =='MIRI':
                fr2 = ttk.Frame(page)
                #label = 'IFU' if iname !='TFI' else 'TF'
                ttk.Label(fr2, text='   IFU wavelen: ', state='disabled').grid(row=0, column=0)
                self.widgets[iname+"_ifu_wavelen"] = ttk.Entry(fr2, width=5) #, disabledforeground="#A0A0A0")
                self.widgets[iname+"_ifu_wavelen"].insert(0, str(self.instrument[iname].monochromatic))
                self.widgets[iname+"_ifu_wavelen"].grid(row=0, column=1)
                self.widgets[iname+"_ifu_wavelen"].state(['disabled'])
                ttk.Label(fr2, text=' um' , state='disabled').grid(row=0, column=2)
                fr2.grid(row=1,column=2, columnspan=6, sticky='E')

                iname2 = iname+"" # need to make a copy so the following lambda function works right:
                self.widgets[iname+"_filter"].bind('<<ComboboxSelected>>', lambda e: self.ev_update_ifu_label(iname2))


            if len(self.instrument[iname].image_mask_list) >0 :
                masks = self.instrument[iname].image_mask_list
                masks.insert(0, "")
 
                self._add_labeled_dropdown(iname+"_coron", page, label='    Coron:', values=masks,  width=12, position=(2,0), sticky='W')
                #self.vars[iname+"_coron"] = tk.StringVar()
                #self.widgets[iname+"_coron"] = ttk.Combobox(page,textvariable =self.vars[iname+"_coron"], width=10, state='readonly')
                #self.widgets[iname+"_coron"]['values'] = masks
                #ttk.Label(page, text='    Coron: ' ).grid(row=2, column=0)
                #self.widgets[iname+"_coron"].set(self.widgets[iname+"_coron"]['values'][0])
                #self.widgets[iname+"_coron"].grid(row=2, column=1)

                #fr2 = ttk.Frame(page)
                #self.vars[iname+"_cor_off_r"] = tk.StringVar()
                #self.vars[iname+"_cor_off_theta"] = tk.StringVar()
                #ttk.Label(fr2, text='target offset:  r=' ).grid(row=2, column=4)
                #self.widgets[iname+"_cor_off_r"] = ttk.Entry(fr2,textvariable =self.vars[iname+"_cor_off_r"], width=5)
                #self.widgets[iname+"_cor_off_r"].insert(0,"0.0")
                #self.widgets[iname+"_cor_off_r"].grid(row=2, column=5)
                #ttk.Label(fr2, text='arcsec,  PA=' ).grid(row=2, column=6)
                #self.widgets[iname+"_cor_off_theta"] = ttk.Entry(fr2,textvariable =self.vars[iname+"_cor_off_theta"], width=3)
                #self.widgets[iname+"_cor_off_theta"].insert(0,"0")
                #self.widgets[iname+"_cor_off_theta"].grid(row=2, column=7)
                #ttk.Label(fr2, text='deg' ).grid(row=2, column=8)
                #fr2.grid(row=2,column=3, sticky='W')


            if len(self.instrument[iname].image_mask_list) >0 :
                masks = self.instrument[iname].pupil_mask_list
                masks.insert(0, "")
                self._add_labeled_dropdown(iname+"_pupil", page, label='    Pupil:', values=masks,  width=12, position=(3,0), sticky='W')

                fr2 = ttk.Frame(page)
                self._add_labeled_entry(iname+"_pupilshift_x", fr2, label='  pupil shift in X:', value='0', width=3, position=(3,4), sticky='W')
                self._add_labeled_entry(iname+"_pupilshift_y", fr2, label=' Y:', value='0', width=3, position=(3,6), sticky='W')

                ttk.Label(fr2, text='% of pupil' ).grid(row=3, column=8)
                fr2.grid(row=3,column=3, sticky='W')


            ttk.Label(page, text='Configuration Options for the OTE').grid(row=4, columnspan=2, sticky='W')
            fr2 = ttk.Frame(page)

            opd_list =  self.instrument[iname].opd_list
            opd_list.insert(0,"Zero OPD (perfect)")
            #if os.getenv("WEBBPSF_ITM") or 1:  
            if self._enable_opdserver:
                opd_list.append("OPD from ITM Server")
            default_opd = self.instrument[iname].pupilopd if self.instrument[iname].pupilopd is not None else "Zero OPD (perfect)"
            self._add_labeled_dropdown(iname+"_opd", fr2, label='    OPD File:', values=opd_list, default=default_opd, width=21, position=(0,0), sticky='W')

            self._add_labeled_dropdown(iname+"_opd_i", fr2, label=' # ', values= [str(i) for i in range(10)], width=3, position=(0,2), sticky='W')

            self.widgets[iname+"_opd_label"] = ttk.Label(fr2, text=' 0 nm RMS            ', width=35)
            self.widgets[iname+"_opd_label"].grid( column=4,sticky='W', row=0)

            self.widgets[iname+"_opd"].bind('<<ComboboxSelected>>', 
                    lambda e: self.ev_update_OPD_labels() )
                    # The below code does not work, and I can't tell why. This only ever has iname = 'FGS' no matter which instrument.
                    # So instead brute-force it with the above to just update all 5. 
                    #lambda e: self.ev_update_OPD_label(self.widgets[iname+"_opd"], self.widgets[iname+"_opd_label"], iname) )
            ttk.Button(fr2, text='Display', command=self.ev_displayOPD).grid(column=5,sticky='E',row=0)

            fr2.grid(row=5, column=0, columnspan=4,sticky='S')



            # ITM interface here - build the widgets now but they will be hidden by default until the ITM option is selected
            fr2 = ttk.Frame(page)
            self._add_labeled_entry(iname+"_coords", fr2, label='    Source location:', value='0, 0', width=12, position=(1,0), sticky='W')
            units_list = ['V1,V2 coords', 'detector pixels']
            self._add_labeled_dropdown(iname+"_coord_units", fr2, label='in:', values=units_list, default=units_list[0], width=11, position=(1,2), sticky='W')
            choose_list=['', 'SI center', 'SI upper left corner', 'SI upper right corner', 'SI lower left corner', 'SI lower right corner']
            self._add_labeled_dropdown(iname+"_coord_choose", fr2, label='or select:', values=choose_list, default=choose_list[0], width=21, position=(1,4), sticky='W')


            ttk.Label(fr2, text='    ITM output:').grid(row=2, column=0, sticky='W')
            self.widgets[iname+"_itm_output"] = ttk.Label(fr2, text='    - no file available yet -')
            self.widgets[iname+"_itm_output"].grid(row=2, column=1, columnspan=4, sticky='W')
            ttk.Button(fr2, text='Access ITM...', command=self.ev_launch_ITM_dialog).grid(column=5,sticky='E',row=2)


            fr2.grid(row=6, column=0, columnspan=4,sticky='SW')
            self.widgets[iname+"_itm_coords"] = fr2


        self.ev_update_OPD_labels()
        lf.grid(row=2, sticky='E,W', padx=10, pady=5)
        notebook.select(0)

        lf = ttk.LabelFrame(frame, text='Calculation Options')
        r =0
        self._add_labeled_entry('FOV', lf, label='Field of View:',  width=3, value='5', postlabel='arcsec/side', position=(r,0))
        r+=1
        self._add_labeled_entry('detector_oversampling', lf, label='Output Oversampling:',  width=3, value='2', postlabel='x finer than instrument pixels       ', position=(r,0))

        #self.vars['downsamp'] = tk.BooleanVar()
        #self.vars['downsamp'].set(True)
        #self.widgets['downsamp'] = ttk.Checkbutton(lf, text='Save in instr. pixel scale, too?', onvalue=True, offvalue=False,variable=self.vars['downsamp'])
        #self.widgets['downsamp'].grid(row=r, column=4, sticky='E')

        output_options=['Oversampled PSF only', 'Oversampled + Detector Res. PSFs', 'Mock full image from JWST DMS']
        self._add_labeled_dropdown("output_type", fr2, label='Output format:', values=output_options, default=output_options[1], width=31, position=(r,4), sticky='W')


        r+=1
        self._add_labeled_entry('fft_oversampling', lf, label='Coronagraph FFT Oversampling:',  width=3, value='2', postlabel='x finer than Nyquist', position=(r,0))
        r+=1
        self._add_labeled_entry('nlambda', lf, label='# of wavelengths:',  width=3, value='', position=(r,0), postlabel='Leave blank for autoselect')
        r+=1

        self._add_labeled_dropdown("jitter", lf, label='Jitter model:', values=  ['Just use OPDs', 'Gaussian - 7 mas rms', 'Gaussian - 30 mas rms'], width=20, position=(r,0), sticky='W', columnspan=2)
        r+=1
        self._add_labeled_dropdown("output_format", lf, label='Output Format:', values=  ['Oversampled image','Detector sampled image','Both as FITS extensions', 'Mock JWST DMS Output' ], width=30, position=(r,0), sticky='W', columnspan=2)
        #self._add_labeled_dropdown("jitter", lf, label='Jitter model:', values=  ['Just use OPDs', 'Gaussian blur', 'Accurate yet SLOW grid'], width=20, position=(r,0), sticky='W', columnspan=2)

        lf.grid(row=4, sticky='E,W', padx=10, pady=5)

        lf = ttk.Frame(frame)

        def addbutton(self,lf, text, command, pos, disabled=False):
            self.widgets[text] = ttk.Button(lf, text=text, command=command )
            self.widgets[text].grid(column=pos, row=0, sticky='E')
            if disabled:
                self.widgets[text].state(['disabled'])

 
        addbutton(self,lf,'Compute PSF', self.ev_calcPSF, 0)
        addbutton(self,lf,'Display PSF', self.ev_displayPSF, 1, disabled=True)
        addbutton(self,lf,'Display profiles', self.ev_displayProfiles, 2, disabled=True)
        addbutton(self,lf,'Save PSF As...', self.ev_SaveAs, 3, disabled=True)
        addbutton(self,lf,'More options...', self.ev_options, 4, disabled=False)

        ttk.Button(lf, text='Quit', command=self.quit).grid(column=5, row=0)
        lf.columnconfigure(2, weight=1)
        lf.columnconfigure(4, weight=1)
        lf.grid(row=5, sticky='E,W', padx=10, pady=15)

        frame.grid(row=0, sticky='N,E,S,W')
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def quit(self):
        " Quit the GUI"
        if tkMessageBox.askyesno( message='Are you sure you want to quit WebbPSF?', icon='question', title='Confirm quit') :
            self.root.destroy()

    def ev_SaveAs(self):
        "Event handler for Save As of output PSFs"
        filename = tkFileDialog.asksaveasfilename(
                initialfile='PSF_%s_%s.fits' %(self.iname, self.filter), 
                filetypes=[('FITS', '.fits')],
                defaultextension='.fits',
                parent=self.root)
        if len(filename) > 0:
            self.PSF_HDUlist.writeto(filename) 
            print("Saved to {}".format(filename))

    def ev_options(self):
        d = WebbPSFOptionsDialog(self.root, input_options = self.advanced_options)
        if d.results is not None: # none means the user hit 'cancel'
            self.advanced_options = d.results

    def ev_plotspectrum(self):
        "Event handler for Plot Spectrum "
        self._updateFromGUI()

        print("Spectral type is "+self.sptype)
        print("Selected instrument tab is "+self.iname)
        print("Selected instrument filter is "+self.filter)


        plt.clf()

        ax1 = plt.subplot(311)
        spectrum = poppy.specFromSpectralType(self.sptype)
        synplot(spectrum)
        ax1.set_ybound(1e-6, 1e8) # hard coded for now
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=1000))
        legend_font = matplotlib.font_manager.FontProperties(size=10)
        ax1.legend(loc='lower right', prop=legend_font)


        ax2 = plt.subplot(312, sharex=ax1)
        ax2.set_ybound(0,1.1)
        band = self.inst._getSynphotBandpass(self.inst.filter) #pysynphot.ObsBandpass(obsname)
        band.name = "%s %s" % (self.iname, self.inst.filter)
        synplot(band) #, **kwargs)
        legend_font = matplotlib.font_manager.FontProperties(size=10)
        plt.legend(loc='lower right', prop=legend_font)

        ax2.set_ybound(0,1.1)


        ax3 = plt.subplot(313, sharex=ax1)
        if self.nlambda is None:
            # Automatically determine number of appropriate wavelengths.
            # Make selection based on filter configuration file
            try:
                nlambda = self.inst._filters[self.filter].default_nlambda
            except KeyError:
                nlambda=10
        else:
            nlambda = self.nlambda
        ax1.set_xbound(0.1, 100)
        plt.draw()
        waves, weights = self.inst._getWeights(spectrum, nlambda=nlambda)

        wave_step = waves[1]-waves[0]
        plot_waves = np.concatenate( ([waves[0]-wave_step], waves, [waves[-1]+wave_step])) * 1e6
        plot_weights = np.concatenate(([0], weights,[0]))


        plt.ylabel("Weight")
        plt.xlabel("Wavelength [$\mu$m]")

        ax3.plot(plot_waves, plot_weights,  drawstyle='steps-mid')

        ax1.set_xbound(0.1, 100)

        self._refresh_window()

    def _refresh_window(self):
        """ Force the window to refresh, and optionally to show itself if hidden (for recent matplotlibs)"""
        plt.draw()
        from distutils.version import StrictVersion
        if StrictVersion(matplotlib.__version__) >= StrictVersion('1.1'):
            plt.show(block=False)

    def ev_calcPSF(self):
        "Event handler for PSF Calculations"
        self._updateFromGUI()

        if _HAS_PYSYNPHOT:
            source = poppy.specFromSpectralType(self.sptype)
        else:
            source=None # generic flat spectrum

        self.PSF_HDUlist = self.inst.calcPSF(source=source, 
                detector_oversample= self.detector_oversampling,
                fft_oversample=self.fft_oversampling,
                fov_arcsec = self.FOV,  nlambda = self.nlambda, display=True)
        #self.PSF_HDUlist.display()
        for w in ['Display PSF', 'Display profiles', 'Save PSF As...']:
           self.widgets[w].state(['!disabled'])
        self._refresh_window()
        _log.info("PSF calculation complete")

    def ev_displayPSF(self):
        "Event handler for Displaying the PSF"
        #self._updateFromGUI()
        #if self.PSF_HDUlist is not None:
        plt.clf()
        poppy.display_PSF(self.PSF_HDUlist, vmin = self.advanced_options['psf_vmin'], vmax = self.advanced_options['psf_vmax'], 
                scale = self.advanced_options['psf_scale'], cmap= self.advanced_options['psf_cmap'], normalize=self.advanced_options['psf_normalize'])
        self._refresh_window()

    def ev_displayProfiles(self):
        "Event handler for Displaying the PSF"
        #self._updateFromGUI()
        poppy.display_profiles(self.PSF_HDUlist)        
        self._refresh_window()

    def ev_displayOptics(self):
        "Event handler for Displaying the optical system"
        self._updateFromGUI()
        _log.info("Selected OPD is "+str(self.opd_name))

        plt.clf()
        self.inst.display()
        self._refresh_window()

    def ev_displayOPD(self):
        self._updateFromGUI()
        if self.inst.pupilopd is None:
            tkMessageBox.showwarning( message="You currently have selected no OPD file (i.e. perfect telescope) so there's nothing to display.", title="Can't Display") 
        else:
            if self._enable_opdserver and 'ITM' in self.opd_name:
                opd = self.inst.pupilopd   # will contain the actual OPD loaded in _updateFromGUI just above
            else:
                opd = fits.getdata(self.inst.pupilopd[0])     # in this case self.inst.pupilopd is a tuple with a string so we have to load it here. 

            if len(opd.shape) >2:
                opd = opd[self.opd_i,:,:] # grab correct slice

            masked_opd = np.ma.masked_equal(opd,  0) # mask out all pixels which are exactly 0, outside the aperture
            cmap = matplotlib.cm.jet
            cmap.set_bad('k', 0.8)
            plt.clf()
            plt.imshow(masked_opd, cmap=cmap, interpolation='nearest', vmin=-0.5, vmax=0.5)
            plt.title("OPD from %s, #%d" %( os.path.basename(self.opd_name), self.opd_i))
            cb = plt.colorbar(orientation='vertical')
            cb.set_label('microns')

            f = plt.gcf()
            plt.text(0.4, 0.02, "OPD WFE = %6.2f nm RMS" % (masked_opd.std()*1000.), transform=f.transFigure)

        self._refresh_window()

    def ev_launch_ITM_dialog(self):
        tkMessageBox.showwarning( message="ITM dialog box not yet implemented", title="Can't Display") 

    def ev_update_OPD_labels(self):
        "Update the descriptive text for all OPD files"
        for iname in self.instrument:
            self.ev_update_OPD_label(self.widgets[iname+"_opd"], self.widgets[iname+"_opd_label"], iname)

    def ev_update_OPD_label(self, widget_combobox, widget_label, iname):
        "Update the descriptive text displayed about one OPD file"
        showitm=False # default is do not show
        filename = self.instrument[iname]._datapath +os.sep+ 'OPD'+ os.sep+widget_combobox.get()
        if filename.endswith(".fits"):
            header_summary = fits.getheader(filename)['SUMMARY']
            self.widgets[iname+"_opd_i"]['state'] = 'readonly'
        else:  # Special options for non-FITS file inputs
            self.widgets[iname+"_opd_i"]['state'] = 'disabled'
            if 'Zero' in widget_combobox.get():
                header_summary = " 0 nm RMS"
            elif 'ITM' in widget_combobox.get() and self._enable_opdserver:
                header_summary= "Get OPD from ITM Server"
                showitm=True
            elif 'ITM' in widget_combobox.get() and not self._enable_opdserver:
                header_summary = "ITM Server is not running or otherwise unavailable."
            else: # other??
                header_summary = "   "

        widget_label.configure(text=header_summary, width=30)


        if showitm:
            self.widgets[iname+"_itm_coords"].grid() # re-show ITM options
        else:
            self.widgets[iname+"_itm_coords"].grid_remove()  # hide ITM options

    def _updateFromGUI(self):
        # get GUI values
        if _HAS_PYSYNPHOT:
            self.sptype = self.widgets['SpType'].get()
        self.iname = self.widgets[self.widgets['tabset'].select()]


 

        try:
            self.nlambda= int(self.widgets['nlambda'].get())
        except ValueError:
            self.nlambda = None # invoke autoselect for nlambda
        self.FOV= float(self.widgets['FOV'].get())
        self.fft_oversampling= int(self.widgets['fft_oversampling'].get())
        self.detector_oversampling= int(self.widgets['detector_oversampling'].get())

        self.output_type = self.widgets['output_type'].get()

        options = {}
        #options['downsample'] = bool(self.vars['downsamp'])
        options['rebin'] = not (self.output_type == 'Oversampled PSF only')  #was downsample, which seems wrong?
        options['mock_dms'] = (self.output_type == 'Mock full image from JWST DMS')
        jitter_choice = self.widgets['jitter'].get()
        if 'Gaussian' in jitter_choice:
            options['jitter'] = 'gaussian'
            if '7 mas' in jitter_choice:
                options['jitter_sigma'] = 0.007
            elif '30 mas' in jitter_choice:
                options['jitter_sigma'] = 0.030
        else:
            options['jitter'] = None

        # and get the values that may have previously been set by the 'advanced options' dialog
        if self.advanced_options is not None:
            for a in self.advanced_options:
                options[a] = self.advanced_options[a]


        # configure the relevant instrument object
        self.inst = self.instrument[self.iname]
        self.filter = self.widgets[self.iname+"_filter"].get() # save for use in default filenames, etc.
        self.inst.filter = self.filter

        self.opd_name = self.widgets[self.iname+"_opd"].get()
        if self._enable_opdserver and 'ITM' in self.opd_name:
            # get from ITM server
            self.opd_i= 0
            self.inst.pupilopd = self._opdserver.get_OPD(return_as="FITS")
            self.opd_name = "OPD from ITM OPD GUI"

        elif self.opd_name == "Zero OPD (perfect)": 
            # perfect OPD
            self.opd_name = "Perfect"
            self.inst.pupilopd = None
        else:
            # Regular FITS file version
            self.opd_name= self.widgets[self.iname+"_opd"].get()
            self.opd_i= int(self.widgets[self.iname+"_opd_i"].get())
            self.inst.pupilopd = (self.inst._datapath+os.sep+"OPD"+os.sep+self.opd_name,self.opd_i)  #filename, slice

        _log.info("Selected OPD is "+str(self.opd_name))


        if self.iname+"_coron" in self.widgets:
            self.inst.image_mask = self.widgets[self.iname+"_coron"].get()
            self.inst.pupil_mask = self.widgets[self.iname+"_pupil"].get()
            # TODO read in mis-registration options here.


            options['source_offset_r'] = float(self.widgets["source_off_r"].get())
            options['source_offset_theta'] = float(self.widgets["source_off_theta"].get())
            options['pupil_shift_x'] = float(self.widgets[self.iname+"_pupilshift_x"].get())/100. # convert from percent to fraction
            options['pupil_shift_y'] = float(self.widgets[self.iname+"_pupilshift_y"].get())/100. # convert from percent to fraction

        self.inst.options = options


    def mainloop(self):
        self.root.mainloop()

#-------------------------------------------------------------------------

class Dialog(tk.Toplevel):
    """ Base class for a modal dialog box. 
    From example code at  http://effbot.org/tkinterbook/tkinter-dialog-windows.htm

    """

    def __init__(self, parent, title = None, input_options=None):

        tk.Toplevel.__init__(self, parent)
        self.transient(parent)
        self.input_options = input_options

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = ttk.Frame(self)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    #
    # construction hooks
    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        pass


    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons

        box = ttk.Frame(self)

        w = ttk.Button(box, text="OK", width=10, command=self.ok, default=tk.ACTIVE)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=tk.LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    #
    # standard button semantics

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):

        return True # override

    def apply(self):

        pass # override


class WebbPSFOptionsDialog(Dialog):
    def _add_labeled_dropdown(self, name, root,label="Entry:", values=None, default=None, width=5, position=(0,0), **kwargs):
        "convenient wrapper for adding a Combobox"

        ttk.Label(root, text=label).grid(row=position[0],  column=position[1], sticky='W')

        self.vars[name] = tk.StringVar()
        self.widgets[name] = ttk.Combobox(root, textvariable=self.vars[name], width=width, state='readonly')
        self.widgets[name].grid(row=position[0], column=position[1]+1, **kwargs)
        self.widgets[name]['values'] = values

        if default is None: default=values[0]
        self.widgets[name].set(default)
 

    def _add_labeled_entry(self, name, root,label="Entry:", value="", width=5, position=(0,0), postlabel=None, **kwargs):
        "convenient wrapper for adding an Entry"
        ttk.Label(root, text=label).grid(row=position[0],  column=position[1], sticky='W')

        self.vars[name] = tk.StringVar()
        self.widgets[name] = ttk.Entry(root, textvariable=self.vars[name], width=width)
        self.widgets[name].insert(0,value)
        self.widgets[name].grid(row=position[0], column=position[1]+1, **kwargs)

        if postlabel is not None:
            ttk.Label(root, text=postlabel).grid(row=position[0],  column=position[1]+2, sticky='W')



    def body(self, master):
        self.results = None # in case we cancel this gets returned
        self.results = None # in case we cancel this gets returned
        self.vars = {}
        self.widgets = {}
        self.values = {}

        colortables = [
         ('Jet (blue to red)',matplotlib.cm.jet),
         ('Gray', matplotlib.cm.gray),
         ('Heat (black-red-yellow)', matplotlib.cm.gist_heat),
         ('Copper (black to tan)',matplotlib.cm.copper),
         ('Stern',matplotlib.cm.gist_stern),
         ('Prism (repeating rainbow)', matplotlib.cm.prism)]

        self.colortables = dict(colortables)


        lf = ttk.LabelFrame(master, text='WebbPSF Advanced Options')

        r=1
        fr2 = ttk.Frame(lf)

        self.values['force_coron'] = ['regular propagation (MFT)', 'full coronagraphic propagation (FFT/SAM)']

        self._add_labeled_dropdown("force_coron", lf, label='Direct imaging calculations use', values=self.values['force_coron'], 
                default = self.values['force_coron'][ 1 if self.input_options['force_coron'] else 0]  ,  width=30, position=(r,0), sticky='W')
        r+=1
        self.values['no_sam'] = ['semi-analytic method if possible', 'basic FFT method always']
        self._add_labeled_dropdown("no_sam", lf, label='Coronagraphic calculations use', values=self.values['no_sam'],
                default=self.values['no_sam'][ 1 if self.input_options['no_sam'] else 0] , width=30, position=(r,0), sticky='W')
        r+=1
        self._add_labeled_dropdown("parity", lf, label='Output pixel grid parity is', values=['odd', 'even', 'either'], default=self.input_options['parity'], width=10, position=(r,0), sticky='W')



        lf.grid(row=1, sticky='E,W', padx=10,pady=5)

        lf = ttk.LabelFrame(master, text='PSF Display Options')
        r=0
        self._add_labeled_dropdown("psf_scale", lf, label='    Display scale:', values=['log','linear'],default=self.input_options['psf_scale'], width=5, position=(r,0), sticky='W')
        r+=1
        self._add_labeled_entry("psf_vmin", lf, label='    Min scale value:', value="%.2g" % self.input_options['psf_vmin'], width=7, position=(r,0), sticky='W')
        r+=1
        self._add_labeled_entry("psf_vmax", lf, label='    Max scale value:', value="%.2g" % self.input_options['psf_vmax'], width=7, position=(r,0), sticky='W')
        r+=1
        self._add_labeled_dropdown("psf_normalize", lf, label='    Normalize PSF to:', values=['Total', 'Peak'], default=self.input_options['psf_normalize'], width=5, position=(r,0), sticky='W')
        r+=1
        self._add_labeled_dropdown("psf_cmap", lf, label='    Color table:', values=[a[0] for a in colortables],  default=self.input_options['psf_cmap_str'], width=20, position=(r,0), sticky='W')
        lf.grid(row=2, sticky='E,W', padx=10,pady=5)


        return self.widgets['force_coron']# return widget to have initial focus

    def apply(self, test=False):
        try:
            results = {}
            results['force_coron'] = self.vars['force_coron'].get() == 'full coronagraphic propagation (FFT/SAM)'
            results['no_sam'] = self.vars['no_sam'].get() == 'basic FFT method always'
            results['parity'] = self.vars['parity'].get() 
            results['psf_scale'] = self.vars['psf_scale'].get() 
            results['psf_vmax'] = float(self.vars['psf_vmax'].get())
            results['psf_vmin'] = float(self.vars['psf_vmin'].get())
            results['psf_cmap_str'] = self.vars['psf_cmap'].get()
            results['psf_normalize'] = self.vars['psf_normalize'].get()
            results['psf_cmap'] = self.colortables[ str(self.vars['psf_cmap'].get() ) ]
        except ValueError:
            return False

        if not test:
            self.results = results
        return True

    def validate(self):
        can_apply = self.apply(test=True)
        if not can_apply:
            _log.error("Invalid entries in one or more fields. Please re-enter!")
        return can_apply



#-------------------------------------------------------------------------

def synplot(thing, waveunit='micron', label=None, **kwargs):
    """ Plot a single PySynPhot object (either SpectralElement or SourceSpectrum)
    versus wavelength, with nice axes labels.

    Really just a simple convenience function.
    """

    # convert to requested display unit.
    wave = thing.waveunits.Convert(thing.wave,waveunit)


    if label is None:
        label = thing.name


    if isinstance(thing, pysynphot.spectrum.SourceSpectrum):
        artist = plt.loglog(wave, thing.flux, label=label, **kwargs)
        plt.xlabel("Wavelength [%s]" % waveunit)
        if str(thing.fluxunits) == 'flam':
            plt.ylabel("Flux [%s]" % ' erg cm$^{-2}$ s$^{-1}$ Ang$^{-1}$' )
        else:
            plt.ylabel("Flux [%s]" % thing.fluxunits)
    elif isinstance(thing, pysynphot.spectrum.SpectralElement):
        artist = plt.plot(wave, thing.throughput,label=label, **kwargs)
        plt.xlabel("Wavelength [%s]" % waveunit)
        plt.ylabel("Throughput")
        plt.gca().set_ylim(0,1)
    else:
        _log.error( "Don't know how to plot that object...")
        artist = None
    return artist




def tkgui(fignum=1):
    # enable log message printout
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    # start the GUI
    gui = WebbPSF_GUI()
    plt.figure(fignum)
    #plt.show(block=False)
    gui.mainloop()



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    gui()



