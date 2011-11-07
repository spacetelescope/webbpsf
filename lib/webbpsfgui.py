#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pyfits
import Tkinter as tk
import tkMessageBox
import tkFileDialog
import logging
#from Tkinter import N,E,S,W


try:
    import ttk
    _use_ttk = True
except:
    _use_ttk = False


try:
    __IPYTHON__
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    pass


try:
    import pysynphot
    _HAS_PYSYNPHOT=True
except:
    _HAS_PYSYNPHOT=False


import webbpsf
__version__ = webbpsf.__version__

class WebbPSF_GUI(object):
    """ A GUI for the PSF Simulator 

    Documentation TBD!

    """
    def __init__(self):
        # init the object and subobjects
        self.instrument = {}
        self.widgets = {}
        self.vars = {}
        self.advanced_options = {'parity': 'any', 'force_coron': False, 'no_sam': False, 'psf_vmin': 1e-8, 'psf_vmax': 1.0, 'psf_scale': 'log', 'psf_cmap_str': 'Jet (blue to red)' , 'psf_normalize': 'Total', 'psf_cmap': matplotlib.cm.jet}
        insts = ['NIRCam', 'NIRSpec','NIRISS', 'MIRI', 'FGS']
        for i in insts:
            self.instrument[i] = webbpsf.Instrument(i)

        if _use_ttk:
            self._create_widgets_py27()
        else:
            self._create_widgets_py26()
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


    def _create_widgets_py27(self):
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
            self._add_labeled_dropdown("SpType", lf, label='    Spectral Type:', values=specFromSpectralType("",return_list=True), default='G0V', width=20, position=(0,0), sticky='W')
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
            self._add_labeled_dropdown(iname+"_filter", page, label='    Filter:', values=self.instrument[iname].filter_list, default=self.instrument[iname].filter, width=10, position=(1,0), sticky='W')
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
 
                self._add_labeled_dropdown(iname+"_coron", page, label='    Coron:', values=masks,  width=10, position=(2,0), sticky='W')
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
                self._add_labeled_dropdown(iname+"_pupil", page, label='    Pupil:', values=masks,  width=10, position=(3,0), sticky='W')

                fr2 = ttk.Frame(page)
                self._add_labeled_entry(iname+"_pupilshift_x", fr2, label='  pupil shift in X:', value='0', width=3, position=(3,4), sticky='W')
                self._add_labeled_entry(iname+"_pupilshift_y", fr2, label=' Y:', value='0', width=3, position=(3,6), sticky='W')

                ttk.Label(fr2, text='% of pupil' ).grid(row=3, column=8)
                fr2.grid(row=3,column=3, sticky='W')


            ttk.Label(page, text='Configuration Options for the OTE').grid(row=4, columnspan=2, sticky='W')
            fr2 = ttk.Frame(page)

            opd_list =  self.instrument[iname].opd_list
            opd_list.insert(0,"Zero OPD (perfect)")
            if os.getenv("WEBBPSF_ITM") or 1:  
                opd_list.append("ITM Server")
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

            # ITM interface here
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
        self._add_labeled_entry('calc_oversampling', lf, label='Coronagraph Oversampling:',  width=3, value='2', postlabel='x finer than Nyquist', position=(r,0))
        r+=1
        self._add_labeled_entry('nlambda', lf, label='# of wavelengths:',  width=3, value='', position=(r,0), postlabel='Leave blank for autoselect')
        r+=1

        self._add_labeled_dropdown("jitter", lf, label='Jitter model:', values=  ['Just use OPDs' ], width=20, position=(r,0), sticky='W', columnspan=2)
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

        #ttk.Button(lf, text='Display Optics', command=self.ev_displayOptics ).grid(column=4, row=0)
        ttk.Button(lf, text='Quit', command=self.quit).grid(column=5, row=0)
        lf.columnconfigure(2, weight=1)
        lf.columnconfigure(4, weight=1)
        lf.grid(row=5, sticky='E,W', padx=10, pady=15)

        frame.grid(row=0, sticky='N,E,S,W')
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def _create_widgets_py26(self):
        """Create a mediocre GUI using the unimpressive widget set 
        available in Python 2.6
        """

        # Query the user what instrument is selected, using
        # a modal dialog widget.

        #---- create the GUIs
        self.root = tk.Tk()
        self.root.geometry('+50+50')
        self.root.title("JWPSF GUI 4alpha")

        frame = ttk.Frame(self.root)

        ttk.Label(frame, text='JWPSF new GUI' ).grid(row=0)

        #-- star
        lf = ttk.LabelFrame(frame, text='Source Properties')
        ttk.Label(lf, text='    Spectral Type' ).grid(row=0, column=0)
        self.vars['SpType'] = tk.StringVar()
        self.widgets['SpType'] = ttk.Combobox(lf, textvariable = self.vars['SpType'],width=6, state='readonly')
        self.widgets['SpType'].grid(row=0, column=1)
        self.widgets['SpType']['values'] = self.stars.sptype_list
        self.widgets['SpType'].set('G0V')
        ttk.Label(lf, text='    ' ).grid(row=0, column=2)
        ttk.Button(lf, text='Plot spectrum', command=self.ev_plotspectrum).grid(row=0,column=2,sticky='E',columnspan=4)

        r = 1
        ttk.Label(lf, text='    Source Position: r=' ).grid(row=r, column=0)
        fr2 = ttk.Frame(lf)

        self.vars["source_off_r"] = tk.StringVar()
        self.vars["source_off_theta"] = tk.StringVar()
        self.vars["source_off_centerpos"] = tk.StringVar()
        self.vars["source_off_centerpos"].set('corner')
        self.widgets["source_off_r"] = ttk.Entry(fr2,textvariable =self.vars["source_off_r"], width=5)
        self.widgets["source_off_r"].insert(0,"0.0")
        self.widgets["source_off_r"].grid(row=r, column=1, sticky='W')
        ttk.Label(fr2, text='arcsec,  PA=' ).grid(row=r, column=2)
        self.widgets["source_off_theta"] = ttk.Entry(fr2,textvariable =self.vars["source_off_theta"], width=3)
        self.widgets["source_off_theta"].insert(0,"0")
        self.widgets["source_off_theta"].grid(row=r, column=3)
        ttk.Label(fr2, text='deg, centered on ' ).grid(row=r, column=4)
        pixel = ttk.Radiobutton(fr2, text='pixel', variable=self.vars["source_off_centerpos"], value='pixel')
        pixel.grid(row=r, column=5)
        corner = ttk.Radiobutton(fr2, text='corner', variable=self.vars["source_off_centerpos"], value='corner')
        corner.grid(row=r, column=6)

        fr2.grid(row=r, column=1, columnspan=5, sticky='W')



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
            ttk.Label(page, text='Configuration Options for '+iname+"                      ").grid(row=0, columnspan=2, sticky='W')

            self.vars[iname+"_filter"] = tk.StringVar()
            self.widgets[iname+"_filter"] = ttk.Combobox(page,textvariable =self.vars[iname+"_filter"], width=10, state='readonly')
            self.widgets[iname+"_filter"]['values'] = self.instrument[iname].filter_list
            self.widgets[iname+"_filter"].set(self.instrument[iname].filter)
            #self.widgets[iname+"_filter"]['readonly'] = True
            ttk.Label(page, text='    Filter: ' ).grid(row=1, column=0)
            self.widgets[iname+"_filter"].grid(row=1, column=1)


            if hasattr(self.instrument[iname], 'ifu_wavelength'):
                fr2 = ttk.Frame(page)
                #label = 'IFU' if iname !='TFI' else 'TFI'
                ttk.Label(fr2, text='   IFU wavelen: ' , state='disabled').grid(row=0, column=0)
                self.widgets[iname+"_ifu_wavelen"] = ttk.Entry(fr2, width=5) #, disabledforeground="#A0A0A0")
                self.widgets[iname+"_ifu_wavelen"].insert(0, str(self.instrument[iname].ifu_wavelength))
                self.widgets[iname+"_ifu_wavelen"].grid(row=0, column=1)
                self.widgets[iname+"_ifu_wavelen"].state(['disabled'])
                ttk.Label(fr2, text=' um' , state='disabled').grid(row=0, column=2)
                fr2.grid(row=1,column=2, columnspan=6, sticky='E')

            if len(self.instrument[iname].image_mask_list) >0 :
                self.vars[iname+"_coron"] = tk.StringVar()
                self.widgets[iname+"_coron"] = ttk.Combobox(page,textvariable =self.vars[iname+"_coron"], width=10, state='readonly')
                masks = self.instrument[iname].image_mask_list
                masks.insert(0, "")
                self.widgets[iname+"_coron"]['values'] = masks
                ttk.Label(page, text='    Coron: ' ).grid(row=2, column=0)
                self.widgets[iname+"_coron"].set(self.widgets[iname+"_coron"]['values'][0])
                self.widgets[iname+"_coron"].grid(row=2, column=1)

            if len(self.instrument[iname].image_mask_list) >0 :
                self.vars[iname+"_pupil"] = tk.StringVar()
                self.widgets[iname+"_pupil"] = ttk.Combobox(page,textvariable =self.vars[iname+"_pupil"], width=10, state='readonly')
                masks = self.instrument[iname].pupil_mask_list
                masks.insert(0, "")
                self.widgets[iname+"_pupil"]['values'] = masks
                ttk.Label(page, text='    Lyot: ' ).grid(row=3, column=0)
                self.widgets[iname+"_pupil"].set(self.widgets[iname+"_pupil"]['values'][0])
                self.widgets[iname+"_pupil"].grid(row=3, column=1)

                fr2 = ttk.Frame(page)
                self.vars[iname+"_pupilshift_x"] = tk.StringVar()
                self.vars[iname+"_pupilshift_y"] = tk.StringVar()
                ttk.Label(fr2, text='  pupil shift in X:' ).grid(row=3, column=4)
                self.widgets[iname+"_pupilshift_x"] = ttk.Entry(fr2,textvariable =self.vars[iname+"_pupilshift_x"], width=3)
                self.widgets[iname+"_pupilshift_x"].insert(0,"0")
                self.widgets[iname+"_pupilshift_x"].grid(row=3, column=5)
                ttk.Label(fr2, text='Y:' ).grid(row=3, column=6)
                self.widgets[iname+"_pupilshift_y"] = ttk.Entry(fr2,textvariable =self.vars[iname+"_pupilshift_y"], width=3)
                self.widgets[iname+"_pupilshift_y"].insert(0,"0")
                self.widgets[iname+"_pupilshift_y"].grid(row=3, column=7)
                ttk.Label(fr2, text='% of pupil' ).grid(row=3, column=8)
                fr2.grid(row=3,column=3, sticky='W')


            ttk.Label(page, text='Configuration Options for the OTE').grid(row=4, columnspan=2, sticky='W')

            fr2 = ttk.Frame(page)
            ttk.Label(fr2, text='   OPD file: ' ).grid(row=0, column=0)
            self.vars[iname+"_opd"] = tk.StringVar()
            self.widgets[iname+"_opd"] = ttk.Combobox(fr2,textvariable =self.vars[iname+"_opd"], width=21, state='readonly')
            opd_list =  self.instrument[iname].opd_list
            opd_list.insert(0,"Zero OPD (perfect)")
            self.widgets[iname+"_opd"]['values'] = opd_list
            #print iname
            #print self.widgets[iname+"_opd"]['values'][0]
            #print self.instrument[iname].pupilopd
            #print "-----"
            self.widgets[iname+"_opd"].set( self.instrument[iname].pupilopd )
            self.widgets[iname+"_opd"].grid(column=1,row=0)

            ttk.Label(fr2, text=' # ' ).grid(column=2,row=0)
            self.vars[iname+"_opd_i"] = tk.StringVar()
            self.widgets[iname+"_opd_i"] = ttk.Combobox(fr2,textvariable =self.vars[iname+"_opd_i"], width=3, state='disabled')
            self.widgets[iname+"_opd_i"]['values'] = [str(i+1) for i in range(10)]
            self.widgets[iname+"_opd_i"].set(self.widgets[iname+"_opd_i"]['values'][0])
            self.widgets[iname+"_opd_i"].grid(column=3,row=0)

            self.widgets[iname+"_opd_label"] = ttk.Label(fr2, text=' 0 nm RMS            ', width=30)
            self.widgets[iname+"_opd_label"].grid( column=4,sticky='W', row=0)


            self.widgets[iname+"_opd"].bind('<<ComboboxSelected>>', 
                    lambda e: self.ev_update_OPD_labels() )
                    # The below code does not work, and I can't tell why. This only ever has iname = 'FGS' no matter which instrument.
                    # So instead brute-force it with the above to just update all 5. 
                    #lambda e: self.ev_update_OPD_label(self.widgets[iname+"_opd"], self.widgets[iname+"_opd_label"], iname) )
            ttk.Button(fr2, text='Display', command=self.ev_displayOPD).grid(column=5,sticky='E',row=0)

            fr2.grid(row=5, column=0, columnspan=4,sticky='S')




        self.ev_update_OPD_labels()
        lf.grid(row=2, sticky='E,W', padx=10, pady=5)
        notebook.select(0)
        #print notebook.tabs()





        #notebook.tab('NIRCam').focus_set()

        #lf = ttk.LabelFrame(frame, text='OTE Config')
        #ttk.Label(lf, text='OTE OPD' ).grid(row=0)
        #lf.grid(row=3, sticky='E'+W, padx=10, pady=5)


        lf = ttk.LabelFrame(frame, text='Calculation Options')
        r =0
        ttk.Label(lf, text='Field of View:' ).grid(row=r, sticky='W')
        self.widgets['FOV'] = ttk.Entry(lf, width=3)
        self.widgets['FOV'].grid(row=r,column=1, sticky='E')
        self.widgets['FOV'].insert(0,'5')
        ttk.Label(lf, text='arcsec/side' ).grid(row=r, column=2, sticky='W')
        r+=1
        ttk.Label(lf, text='Output Oversampling:').grid(row=r, sticky='W')
        self.widgets['detector_oversampling'] = ttk.Entry(lf, width=3)
        self.widgets['detector_oversampling'].grid(row=r,column=1, sticky='E')
        self.widgets['detector_oversampling'].insert(0,'2')
        ttk.Label(lf, text='x finer than instrument pixels       ' ).grid(row=r, column=2, sticky='W', columnspan=2)

        self.vars['downsamp'] = tk.BooleanVar()
        self.vars['downsamp'].set(True)
        self.widgets['downsamp'] = ttk.Checkbutton(lf, text='Save in instr. pixel scale, too?', onvalue=True, offvalue=False,variable=self.vars['downsamp'])
        self.widgets['downsamp'].grid(row=r, column=4, sticky='E')


        r+=1
        ttk.Label(lf, text='Coronagraph Oversampling:').grid(row=r, sticky='W')
        self.widgets['calc_oversampling'] = ttk.Entry(lf, width=3)
        self.widgets['calc_oversampling'].grid(row=r,column=1, sticky='E')
        self.widgets['calc_oversampling'].insert(0,'2')
        ttk.Label(lf, text='x finer than Nyquist' ).grid(row=r, column=2, sticky='W', columnspan=2)


        r+=1
        ttk.Label(lf, text='# of wavelengths:').grid(row=r, sticky='W')
        self.widgets['nlambda'] = ttk.Entry(lf, width=3)
        self.widgets['nlambda'].grid(row=r,column=1, sticky='E')
        #self.widgets['nlambda'].insert(0,'1')
        r+=1
        ttk.Label(lf, text='Jitter model:').grid(row=r, sticky='W')
        self.widgets['jitter'] = ttk.Combobox(lf, width=20, state='readonly')
        self.widgets['jitter']['values'] = ['Just use OPDs', 'Gaussian blur', 'Accurate yet SLOW grid']
        self.widgets['jitter'].set('Just use OPDs')
        self.widgets['jitter'].grid(row=r,column=1, columnspan=2, sticky='E')
 
        #ttk.Label(lf, text='x finer than instrument pixels' ).grid(row=r, column=2, sticky='W', columnspan=2)
 
        #r=r+1
        #ttk.Label(lf, text='Output filename:',  ).grid(row=r, sticky='W')
        #self.widgets['outfile'] = ttk.Entry(lf, width=40)
        #self.widgets['outfile'].grid(row=r,column=1, sticky='E', columnspan=2)
        #ttk.Label(lf, text='.fits' ).grid(row=r, column=3, sticky='W')


        lf.grid(row=4, sticky='E,W', padx=10, pady=5)

        lf = ttk.Frame(frame)
        ttk.Button(lf, text='Compute PSF', command=self.ev_calcPSF ).grid(column=0, row=0)
        self.widgets['SaveAs'] = ttk.Button(lf, text='Save PSF...', command=self.ev_SaveAs )
        self.widgets['SaveAs'].grid(column=1, row=0, sticky='E')
        self.widgets['SaveAs'].state(['disabled'])
        #ttk.Button(lf, text='Display PSF', command=self.ev_displayPSF).grid(column=1, row=0)
        ttk.Button(lf, text='Display Optics', command=self.ev_displayOptics ).grid(column=3, row=0)
        ttk.Button(lf, text='Quit', command=self.quit).grid(column=5, row=0)
        lf.columnconfigure(2, weight=1)
        lf.columnconfigure(4, weight=1)
        lf.grid(row=5, sticky='E,W', padx=10, pady=15)



        #frame.pack(padx=10, pady=10)
        #frame.grid(row=0, sticky=N+E+S+W, padx=10, pady=10)
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
            print "Saved to %s" % filename

    def ev_options(self):
        d = WebbPSFOptionsDialog(self.root, input_options = self.advanced_options)
        if d.results is not None: # none means the user hit 'cancel'
            self.advanced_options = d.results

    def ev_plotspectrum(self):
        "Event handler for Plot Spectrum "
        self._updateFromGUI()

        #sptype = self.widgets['SpType'].get()
        #iname = self.widgets[self.widgets['tabset'].select()]
        print "Spectral type is "+self.sptype
        print "Selected instrument tab is "+self.iname
        #if iname != 'TFI':
            #filter = self.widgets[self.iname+"_filter"].get()
        print "Selected instrument filter is "+self.filter


        plt.clf()

        ax1 = plt.subplot(311)
        spectrum = specFromSpectralType(self.sptype)
        synplot(spectrum)
        ax1.set_ybound(1e-6, 1e8) # hard coded for now
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=1000))

        ax2 = plt.subplot(312, sharex=ax1)
        ax2.set_ybound(0,1.1)
        #try:
        band = self.inst._getSynphotBandpass() #pysynphot.ObsBandpass(obsname)
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
                #if self.inst.name=='TFI':    # filter names are irrelevant for TFI.
                    #nlambda=5
                #else:
                    #filt_width = self.filter[-1]
                    #lookup_table = {'NIRCam': {'2': 10, 'W':20,'M':3,'N':1}, 
                                    #'NIRSpec':{'W':5,'M':3,'N':1}, 
                                    #'MIRI':{'W':5,'M':3,'N':1}, 
                                    #'FGS':{'W':5,'M':3,'N':1}}

                    #nlambda = lookup_table[self.name][filt_width]
                nlambda = self.inst._filter_nlambda_default[self.filter]
            except:
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


        #speclib.specFromSpectralType(sptype)
        #plt.loglog(spectrum['wavelength_um'],spectrum['flux']/spectrum['flux'].max(),label=sptype+" star")
        #plt.xlabel("Wavelength [$\mu$m]")
        #plt.ylabel("Flux [$erg s^{-1} cm^{-2} \AA^{-1} \\times$ arbitrary scale factor ]")

        #filt = self.instrument[iname].getFilter(filter)
        #plt.plot(filt.WAVELENGTH, filt.THROUGHPUT+1e-9, "--" ,label=filter+" filter")
        #plt.gca().set_ybound(1e-6,10)
        #plt.legend(loc="upper left")

    def ev_calcPSF(self):
        "Event handler for PSF Calculations"
        self._updateFromGUI()

        if _HAS_PYSYNPHOT:
            source = specFromSpectralType(self.sptype)
        else:
            source=None # generic flat spectrum

        self.PSF_HDUlist = self.inst.calcPSF(source=source, 
                detector_oversample= self.detector_oversampling,
                calc_oversample=self.calc_oversampling,
                fov_arcsec = self.FOV,  nlambda = self.nlambda, display=True)
        #self.PSF_HDUlist.display()
        for w in ['Display PSF', 'Display profiles', 'Save PSF As...']:
           self.widgets[w].state(['!disabled'])

    def ev_displayPSF(self):
        "Event handler for Displaying the PSF"
        #self._updateFromGUI()
        #if self.PSF_HDUlist is not None:
        plt.clf()
        webbpsf.display_PSF(self.PSF_HDUlist, vmin = self.advanced_options['psf_vmin'], vmax = self.advanced_options['psf_vmax'], 
                scale = self.advanced_options['psf_scale'], cmap= self.advanced_options['psf_cmap'], normalize=self.advanced_options['psf_normalize'])

    def ev_displayProfiles(self):
        "Event handler for Displaying the PSF"
        #self._updateFromGUI()
        webbpsf.display_profiles(self.PSF_HDUlist)        

    def ev_displayOptics(self):
        "Event handler for Displaying the optical system"
        self._updateFromGUI()
        plt.clf()
        self.inst.display()

    def ev_displayOPD(self):
        self._updateFromGUI()
        if self.inst.pupilopd is None:
            tkMessageBox.showwarning( message="You currently have selected no OPD file (i.e. perfect telescope) so there's nothing to display.", title="Can't Display") 
        else:
            opd = pyfits.getdata(self.inst.pupilopd[0])

            opd = opd[self.opd_i,:,:] # grab correct slice

            masked_opd = np.ma.masked_equal(opd,  0) # mask out all pixels which are exactly 0, outside the aperture
            cmap = matplotlib.cm.jet
            cmap.set_bad('k', 0.8)
            plt.clf()
            plt.imshow(masked_opd, cmap=cmap, interpolation='nearest', vmin=-0.5, vmax=0.5)
            plt.title("OPD from %s, #%d" %( os.path.basename(self.opd), self.opd_i))
            cb = plt.colorbar(orientation='vertical')
            cb.set_label('microns')

            f = plt.gcf()
            plt.text(0.4, 0.02, "OPD WFE = %6.2f nm RMS" % (masked_opd.std()*1000.), transform=f.transFigure)

    def ev_update_ifu_label(self, iname):
        pass
        # FIXME
        # I no longer remember what the point of this function was... it doesn't appear to do anything??
        # Disabling for now - remove in future version?

        #newfilt = self.widgets[iname+"_filter"].get()
        #print "Updating IFU label for "+iname+", filt="+newfilt

    def ev_launch_ITM_dialog(self):
        tkMessageBox.showwarning( message="ITM dialog box not yet implemented", title="Can't Display") 


    def ev_update_OPD_labels(self):
        "Update the descriptive text for all OPD files"
        for iname in self.instrument.keys():
            self.ev_update_OPD_label(self.widgets[iname+"_opd"], self.widgets[iname+"_opd_label"], iname)

    def ev_update_OPD_label(self, widget_combobox, widget_label, iname):
        "Update the descriptive text displayed about one OPD file"
        #print "Here! for "+iname
        filename = self.instrument[iname]._datapath +os.sep+ 'OPD'+ os.sep+widget_combobox.get()
        if filename.endswith(".fits"):
            header_summary = pyfits.getheader(filename)['SUMMARY']
            self.widgets[iname+"_opd_i"]['state'] = 'readonly'
            self.widgets[iname+"_itm_coords"].grid_remove()  # hide ITM options
        else:  # Special options for non-FITS file inputs
            self.widgets[iname+"_opd_i"]['state'] = 'disabled'
            if 'Zero' in widget_combobox.get():
                header_summary = " 0 nm RMS"
                self.widgets[iname+"_itm_coords"].grid_remove()  # hide ITM options
            elif 'ITM' in widget_combobox.get():
                header_summary= "Get OPD from ITM Server"

                self.widgets[iname+"_itm_coords"].grid() # re-show
            else:
                header_summary = "   "
                self.widgets[iname+"_itm_coords"].grid_remove()  # hide

        widget_label.configure(text=header_summary, width=30)

    def _updateFromGUI(self):
        # get GUI values
        if _HAS_PYSYNPHOT:
            self.sptype = self.widgets['SpType'].get()
        self.iname = self.widgets[self.widgets['tabset'].select()]
        self.opd= self.widgets[self.iname+"_opd"].get()
        self.opd_i= int(self.widgets[self.iname+"_opd_i"].get())
        try:
            self.nlambda= int(self.widgets['nlambda'].get())
        except:
            self.nlambda = None # invoke autoselect for nlambda
        self.FOV= float(self.widgets['FOV'].get())
        self.calc_oversampling= int(self.widgets['calc_oversampling'].get())
        self.detector_oversampling= int(self.widgets['detector_oversampling'].get())

        self.output_type = self.widgets['output_type'].get()

        options = {}
        #options['downsample'] = bool(self.vars['downsamp'])
        options['rebin'] = not (self.output_type == 'Oversampled PSF only')  #was downsample, which seems wrong?
        options['mock_dms'] = (self.output_type == 'Mock full image from JWST DMS')
        options['jitter'] = self.widgets['jitter'].get()
        #print "Downsamp value: ",  options['downsample']

        # and get the values that may have previously been set by the 'advanced options' dialog
        if self.advanced_options is not None:
            for a in self.advanced_options.keys():
                options[a] = self.advanced_options[a]


        # configure the relevant instrument object
        self.inst = self.instrument[self.iname]
        if self.iname != 'TFI':
            self.filter = self.widgets[self.iname+"_filter"].get()
            self.inst.filter = self.filter
        else:
            self.wavelen = float(self.widgets[self.iname+"_wavelen"].get())
            self.inst.etalon_wavelength = self.wavelen
            self.filter = '%.3fum' % self.wavelen # save for use in save as filename

            
        if self.opd == "Zero OPD (perfect)": 
            self.opd = None
            self.inst.pupilopd = None
        else:
            self.inst.pupilopd = (self.inst._datapath+os.sep+"OPD"+os.sep+self.opd,self.opd_i)  #filename, slice
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
        except:
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

def synplot(thing, waveunit='micron', **kwargs):
    """ Plot a single PySynPhot object (either SpectralElement or SourceSpectrum)
    versus wavelength, with nice axes labels.

    Really just a simple convenience function.
    """

    # convert to requested display unit.
    wave = thing.waveunits.Convert(thing.wave,waveunit)


    if isinstance(thing, pysynphot.spectrum.SourceSpectrum):
        plt.loglog(wave, thing.flux, label=thing.name, **kwargs)
        plt.xlabel("Wavelength [%s]" % waveunit)
        if str(thing.fluxunits) == 'flam':
            plt.ylabel("Flux [%s]" % ' erg cm$^{-2}$ s$^{-1}$ Ang$^{-1}$' )
        else:
            plt.ylabel("Flux [%s]" % thing.fluxunits)
    elif isinstance(thing, pysynphot.spectrum.SpectralElement):
        plt.plot(wave, thing.throughput,label=thing.name, **kwargs)
        plt.xlabel("Wavelength [%s]" % waveunit)
        plt.ylabel("Throughput")
        plt.gca().set_ylim(0,1)
    else:
        print "Don't know how to plot that object..."



def specFromSpectralType(sptype, return_list=False):
    """Get Pysynphot Spectrum object from a spectral type string.

    """
    lookuptable = {
        "O3V":   (50000, 0.0, 5.0),
        "O5V":   (45000, 0.0, 5.0),
        "O6V":   (40000, 0.0, 4.5),
        "O8V":   (35000, 0.0, 4.0),
        "O5I":   (40000, 0.0, 4.5),
        "O6I":   (40000, 0.0, 4.5),
        "O8I":   (34000, 0.0, 4.0),
        "B0V":   (30000, 0.0, 4.0),
        "B3V":   (19000, 0.0, 4.0),
        "B5V":   (15000, 0.0, 4.0),
        "B8V":   (12000, 0.0, 4.0),
        "B0III": (29000, 0.0, 3.5),
        "B5III": (15000, 0.0, 3.5),
        "B0I":   (26000, 0.0, 3.0),
        "B5I":   (14000, 0.0, 2.5),
        "A0V":   (9500, 0.0, 4.0),
        "A5V":   (8250, 0.0, 4.5),
        "A0I":   (9750, 0.0, 2.0),
        "A5I":   (8500, 0.0, 2.0),
        "F0V":   (7250, 0.0, 4.5),
        "F5V":   (6500, 0.0, 4.5),
        "F0I":   (7750, 0.0, 2.0),
        "F5I":   (7000, 0.0, 1.5),
        "G0V":   (6000, 0.0, 4.5),
        "G5V":   (5750, 0.0, 4.5),
        "G0III": (5750, 0.0, 3.0),
        "G5III": (5250, 0.0, 2.5),
        "G0I":   (5500, 0.0, 1.5),
        "G5I":   (4750, 0.0, 1.0),
        "K0V":   (5250, 0.0, 4.5),
        "K5V":   (4250, 0.0, 4.5),
        "K0III": (4750, 0.0, 2.0),
        "K5III": (4000, 0.0, 1.5),
        "K0I":   (4500, 0.0, 1.0),
        "K5I":   (3750, 0.0, 0.5),
        "M0V":   (3750, 0.0, 4.5),
        "M2V":   (3500, 0.0, 4.5),
        "M5V":   (3500, 0.0, 5.0),
        "M0III": (3750, 0.0, 1.5),
        "M0I":   (3750, 0.0, 0.0),
        "M2I":   (3500, 0.0, 0.0)}


    if return_list:
        sptype_list = lookuptable.keys()
        def sort_sptype(typestr):
            letter = typestr[0]
            lettervals = {'O':0, 'B': 10, 'A': 20,'F': 30, 'G':40, 'K': 50, 'M':60}
            value = lettervals[letter]*1.0
            value += int(typestr[1])
            if "III" in typestr: value += .3
            elif "I" in typestr: value += .1
            elif "V" in typestr: value += .5
            return value
        sptype_list.sort(key=sort_sptype)
        sptype_list.insert(0,"Flat spectrum in F_nu")
        sptype_list.insert(0,"Flat spectrum in F_lambda")
        return sptype_list

    try:
        keys = lookuptable[sptype]
        return pysynphot.Icat('ck04models',keys[0], keys[1], keys[2])
    except:
        if "Flat" in sptype:
            if sptype == "Flat spectrum in F_nu":    spec = pysynphot.FlatSpectrum( 1, fluxunits = 'fnu')
            elif sptype == "Flat spectrum in F_lambda":  spec= pysynphot.FlatSpectrum( 1, fluxunits = 'flam')
            spec.convert('flam')
            return spec/spec.flux.mean()
        else: raise LookupError("Lookup table does not include spectral type %s" % sptype)




def gui():
    gui = WebbPSF_GUI()
    gui.mainloop()



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    gui()



