#!/usr/bin/env python
from __future__ import division, print_function, absolute_import, unicode_literals
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.io.fits as fits

from threading import Thread

from . import config
from . import conf

__doc__ = """
Graphical Interface for WebbPSF 

Developed in wxpython by Marshall Perrin

Some code based on/borrowed from Scamp by Klaus Pontoppidan

"""

try:
    import wx, wx.html
except ImportError:
    raise ImportError("The wxPython module is required to run this program.")


import logging
import logging
_log = logging.getLogger('webbpsf')



def _default_options():
    import poppy
    return {'force_coron': False, 'no_sam': False, 'parity':'Either',
                'psf_scale':'log', 'psf_normalize':'Peak', 
                'psf_cmap_str': 'Jet (blue to red)', 'psf_cmap': matplotlib.cm.jet,
                'psf_vmin':1e-4, 'psf_vmax':1.0, 'monochromatic': False, 'fov_in_arcsec': True,
                'parallelization': poppy.conf.use_multiprocessing }



try:
    import pysynphot
    _HAS_PYSYNPHOT=True
except:
    _HAS_PYSYNPHOT=False


import poppy
import webbpsf_core

class WebbPSF_GUI(wx.Frame):
    """ A GUI for the PSF Simulator 

    Documentation TBD!

    """
    def __init__(self, parent=None, id=-1, title="WebbPSF: JWST PSF Calculator"): 
        wx.Frame.__init__(self,parent,id=id,title=title)
        self.parent=parent
        opdserver=None
        self.log_window=None
        # init the object and subobjects
        self.instrument = {}
        self.widgets = {}
        self.vars = {}
        self.advanced_options = _default_options()
        insts = ['NIRCam', 'NIRSpec','NIRISS', 'MIRI', 'FGS']
        for i in insts:
            self.instrument[i] = webbpsf_core.Instrument(i)

        self.inst = self.instrument['NIRCam'] # default

        #invoke link to ITM server if provided?
        if opdserver is not None:
            self._enable_opdserver = True
            self._opdserver = opdserver
        else:
            self._enable_opdserver = False

        
        # create widgets & run
        self._create_widgets_wx()

        # Set up event handler for any worker thread results
        EVT_RESULT(self,self.OnCalcDone)

        self.Show(True)

    def _add_labeled_dropdown(self, name, parent, parentsizer, label="Entry:", choices=None, default=0, width=None, position=(0,0), columnspan=1, expand=True, **kwargs):
        "convenient wrapper for adding a Combobox"

        mylabel = wx.StaticText(parent, -1,label=label)
        parentsizer.Add( mylabel, position,  (1,1), wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)

        if choices is None:
            try:
                choices = self.values[name]
            except:
                choices = ['Option A','Option B']

        if isinstance(default,int):
            value = choices[default]
        else:
            value=default

        #size = (width,25) if width is not None else None
        size=None

        mycombo = wx.ComboBox(parent, -1,value=value, choices=choices, style=wx.CB_DROPDOWN|wx.CB_READONLY, size=size)
        style = wx.EXPAND if expand else 0
        parentsizer.Add( mycombo, (position[0],position[1]+1),  (1,columnspan), style)

        self.widgets[name] = mycombo
 

    def _add_labeled_entry(self, name, parent, parentsizer, label="Entry:", value=None, format="%.2g", 
            width=None, position=(0,0), postlabel=None, **kwargs):
        "convenient wrapper for adding an Entry"
        mylabel = wx.StaticText(parent, -1,label=label)
        parentsizer.Add( mylabel, position,  (1,1), wx.EXPAND)
        self.widgets[name+"_pre_label"] = mylabel



        if value is None:
            try:
                value=format % self.input_options[name]
            except:
                value=""

        size = (width,23) if width is not None else None

        mytext = wx.TextCtrl(parent, -1,value=value, size=size)
        parentsizer.Add( mytext, (position[0],position[1]+1),  (1,1), wx.EXPAND)

        self.widgets[name] = mytext

        if postlabel is not None:
            mylabel2 = wx.StaticText(parent, -1,label=postlabel)
            parentsizer.Add( mylabel2, (position[0],position[1]+2),  (1,1), wx.EXPAND)
            self.widgets[name+"_post_label"] = mylabel2


    def _create_widgets_wx(self):
        """Create a nice GUI using the enhanced widget set provided by 
        the ttk extension to Tkinter, available in Python 2.7 or newer
        """
        #---- create the GUIs
        insts = ['NIRCam', 'NIRSpec','NIRISS', 'MIRI',  'FGS']
        self.sb =self.CreateStatusBar()

        menuBar = WebbPSFMenuBar(self)
        self.SetMenuBar(menuBar)

        self.Bind(wx.EVT_MENU, self.ev_AboutBox, id=wx.ID_ABOUT )
        self.Bind(wx.EVT_MENU, self.ev_Preferences, id=menuBar.ids['preferences'] )
        self.Bind(wx.EVT_MENU, self.ev_SaveAs, id=menuBar.ids['save_psf'] )
        self.Bind(wx.EVT_MENU, self.ev_SaveProfiles, id=menuBar.ids['save_profile'] )
        self.Bind(wx.EVT_MENU, self.ev_calcPSF, id=menuBar.ids['calcPSF'] )
        self.Bind(wx.EVT_MENU, self.ev_options, id=menuBar.ids['calc_options'] )
        self.Bind(wx.EVT_MENU, self.ev_showDocs, id=menuBar.ids['documentation'] )
        self.Bind(wx.EVT_MENU, self.ev_plotspectrum, id=menuBar.ids['display_spectrum'] )
        self.Bind(wx.EVT_MENU, self.ev_displayOptics, id=menuBar.ids['display_optics'] )
        self.Bind(wx.EVT_MENU, self.ev_displayOPD, id=menuBar.ids['display_opd'] )
        self.Bind(wx.EVT_MENU, self.ev_displayPSF, id=menuBar.ids['display_psf'] )
        self.Bind(wx.EVT_MENU, self.ev_displayProfiles, id=menuBar.ids['display_profiles'] )


        top_panel = wx.Panel(self)


        topSizer = wx.BoxSizer(wx.VERTICAL)


        #frame = ttk.Frame(self.root)
        #frame = ttk.Frame(self.root, padx=10,pady=10)

        #ttk.Label(frame, text='James Webb PSF Calculator' ).grid(row=0)

        #===== source properties panel ======
        sb1 = wx.StaticBox(top_panel, label='Source Properties')
        sb1Sizer = wx.StaticBoxSizer(sb1,wx.VERTICAL)


        if _HAS_PYSYNPHOT:
            spectrumPanel = wx.Panel(top_panel)
            spectrumSizer = wx.GridBagSizer()
            
            try:
                choices = poppy.specFromSpectralType("",return_list=True)
                default='G0V'
            except:
                choices = ['Error: $PYSYN_CDBS does not have any spectral models']
                default = choices[0]

            self._add_labeled_dropdown("SpType", spectrumPanel,spectrumSizer, label='    Spectral Type:     ', 
                    choices=choices, default=default, 
                    position=(0,0))
            self.ButtonPlotSpec = wx.Button(spectrumPanel, label='Plot Spectrum')
            self.Bind(wx.EVT_BUTTON, self.ev_plotspectrum, self.ButtonPlotSpec)
            spectrumSizer.Add(self.ButtonPlotSpec, (0,3),(1,1))
            spectrumSizer.AddGrowableCol(2)

            spectrumPanel.SetSizerAndFit(spectrumSizer)
            sb1Sizer.Add(spectrumPanel, 0, wx.ALL|wx.EXPAND, 2)

        posPanel = wx.Panel(top_panel)
        posSizer = wx.GridBagSizer()
        r=0
        self._add_labeled_entry("source_off_r", posPanel,posSizer, label='    Source Position: r=', value='0.0', width=60,  position=(r,0))
        self._add_labeled_entry("source_off_theta", posPanel,posSizer, label='arcsec,  PA=', value='0', width=60, position=(r,2))
        posPanel.SetSizerAndFit(posSizer)
        sb1Sizer.Add(posPanel,0, wx.ALL|wx.EXPAND, 2)


        #===== instruments panels ========
        nb = wx.Notebook(top_panel)
        for iname in insts:
            inst_panel = wx.Panel(nb)
            inst_panel.WebbPSFInstrumentName = iname
            panelSizer = wx.GridBagSizer()


            panelSizer.Add(wx.StaticText(inst_panel, label='Configuration Options for '+iname+"     "), (0,0),(1,3), flag=wx.ALIGN_LEFT)

            instButton = wx.Button(inst_panel, label='Display Optics')
            panelSizer.Add(instButton, (0,4),(1,2), flag=wx.ALIGN_RIGHT)
            self.Bind(wx.EVT_BUTTON, self.ev_displayOptics, instButton)


            self._add_labeled_dropdown(iname+"_filter", inst_panel,panelSizer, label='    Filter:', choices=self.instrument[iname].filter_list, 
                default=self.instrument[iname].filter,  position=(1,0))

            if len(self.instrument[iname].image_mask_list) >0 :
                masks = self.instrument[iname].image_mask_list
                masks.insert(0, "")

                label = '    Slit:    ' if iname=='NIRSpec' else '    Image Stop:'
                self._add_labeled_dropdown(iname+"_coron", inst_panel,panelSizer, label=label, choices=masks,  position=(2,0))


            if len(self.instrument[iname].image_mask_list) >0 :
                masks = self.instrument[iname].pupil_mask_list
                masks.insert(0, "")
                self._add_labeled_dropdown(iname+"_pupil", inst_panel,panelSizer, label='    Pupil Stop:', choices=masks,  position=(3,0))

                fr2 = wx.Panel(inst_panel) 
                fr2Sizer = wx.GridBagSizer()
                self._add_labeled_entry(iname+"_pupilshift_x", fr2,fr2Sizer, label='  pupil shift in X:', value='0', width=40, position=(0,4))
                self._add_labeled_entry(iname+"_pupilshift_y", fr2,fr2Sizer, label=' Y:', value='0', width=40, position=(0,6))

                fr2Sizer.Add(wx.StaticText(fr2, label='% of pupil,' ), (0,8), (1,1))
                self._add_labeled_entry(iname+"_pupil_rot", fr2,fr2Sizer, label=' rotated', value='0', width=40, position=(0,9))
                fr2Sizer.Add(wx.StaticText(fr2, label='deg' ), (0,11), (1,1))
                fr2.SetSizerAndFit(fr2Sizer)

                panelSizer.Add(fr2, (3,3),(1,3))

            

            panelSizer.Add(wx.StaticText(inst_panel, label='Configuration Options for the Telescope (OTE) '), (5,0),(1,5), flag=wx.ALIGN_LEFT)
            opdPanel = wx.Panel(inst_panel)
            opdSizer = wx.GridBagSizer()
            opd_list = self.instrument[iname].opd_list
            opd_list.insert(0,"Zero OPD (perfect)")
            if self._enable_opdserver:
                opd_list.append("OPD from ITM Server")
            default_opd = self.instrument[iname].pupilopd if self.instrument[iname].pupilopd is not None else "Zero OPD (perfect)"

            self._add_labeled_dropdown(iname+"_opd", opdPanel,opdSizer, label='    OPD File:', choices=opd_list, default=default_opd, width=220, position=(0,0), expand=False)
            self.Bind(wx.EVT_COMBOBOX, self.ev_update_OPD_labels, self.widgets[iname+"_opd"])

            self._add_labeled_dropdown(iname+"_opd_i", opdPanel,opdSizer, label=' # ', choices= [str(i) for i in range(10)], width=3, position=(0,2), expand=False)

            self.widgets[iname+"_opd_label"] = wx.StaticText(opdPanel, label=' 0 nm RMS                                            ', style=wx.ALIGN_LEFT|wx.ST_NO_AUTORESIZE)
            opdSizer.Add(self.widgets[iname+"_opd_label"], (0,5),(1,1))
            opdSizer.AddGrowableCol(5)

            instDispButton = wx.Button(opdPanel,label='Display OPD')
            opdSizer.Add(instDispButton, (0,6),(1,1), flag=wx.ALIGN_RIGHT)
            self.Bind(wx.EVT_BUTTON, self.ev_displayOPD, instDispButton)

            opdPanel.SetSizerAndFit(opdSizer)

            panelSizer.Add(opdPanel, (6,0), (1,6), flag=wx.EXPAND|wx.ALL)


            #spacerPanel = wx.Panel(inst_panel)
            #panelSizer.Add(spacerPanel, (7,0), (1,6), flag=wx.EXPAND|wx.ALL)

            panelSizer.AddGrowableCol(2)
            panelSizer.AddGrowableRow(4)
            panelSizer.AddGrowableRow(6)

            inst_panel.SetSizerAndFit(panelSizer)
        
            nb.AddPage(inst_panel, iname)

        self.widgets['tabset'] = nb
        self.ev_update_OPD_labels(None)

        #===== Calculation Options ======

        sb3 = wx.StaticBox(top_panel, label='Calculation Options')
        sb3Sizer = wx.StaticBoxSizer(sb3,wx.VERTICAL)

        calcPanel = wx.Panel(top_panel)
        calcSizer = wx.GridBagSizer()


        r=0 
        self._add_labeled_entry('FOV', calcPanel,calcSizer, label='Field of View:',  value=str(conf.default_fov_arcsec), postlabel='arcsec/side', position=(0,0))
        r+=1
        self._add_labeled_entry('detector_oversampling', calcPanel,calcSizer, label='Output Oversampling:',  width=3, value=str(conf.default_oversampling), postlabel='x finer than instrument pixels       ', position=(r,0))




        r+=1
        self._add_labeled_entry('fft_oversampling', calcPanel,calcSizer, label='Coronagraph FFT Oversampling:',  width=3, value=str(conf.default_oversampling), postlabel='x finer than Nyquist', position=(r,0))
        r+=1
        self._add_labeled_entry('nlambda', calcPanel,calcSizer, label='# of wavelengths:',  width=3, value='', position=(r,0), postlabel='Leave blank for autoselect')
        r+=1

        self._add_labeled_dropdown("jitter", calcPanel,calcSizer, label='Jitter model:', choices=  ['Just use OPDs', 'Gaussian jitter with 7 mas rms', 'Gaussian jitter with 30 mas rms' ], width=20, position=(r,0), columnspan=2)
        r+=1
        output_options=['Oversampled PSF only', 'Oversampled + Detector Res. PSFs', 'Mock full image from JWST DMS']
        self._add_labeled_dropdown("output_format", calcPanel,calcSizer, label='Output Format:', choices=  ['Oversampled image','Detector sampled image','Both as FITS extensions', 'Mock JWST DMS Output' ], width=30, position=(r,0), columnspan=2)

        calcPanel.SetSizerAndFit(calcSizer)
        sb3Sizer.Add(calcPanel,0, wx.ALL|wx.EXPAND, 2)


        #====== button bar ===========

        bbar = wx.Panel(top_panel)
        bbarSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ButtonCompute = wx.Button(bbar, label='Compute PSF')
        bbarSizer.Add(self.ButtonCompute, 1, flag=wx.ALL|wx.EXPAND, border=3)
        self.ButtonSavePSF = wx.Button(bbar, label='Save PSF As...')
        self.ButtonDisplayPSF = wx.Button(bbar, label='Display PSF')
        self.ButtonDisplayProf =wx.Button( bbar, label='Display Profiles')
        self.ButtonSavePSF.Enable(False)
        self.ButtonDisplayPSF.Enable(False)
        self.ButtonDisplayProf.Enable(False)
        self.ButtonOptions = wx.Button(bbar, label='More Options...')
        self.ButtonQuit = wx.Button(bbar, label='Quit')
        self.widgets['Display PSF'] = self.ButtonDisplayPSF
        self.widgets['Display profiles'] = self.ButtonDisplayProf
        self.widgets['Save PSF As...'] = self.ButtonSavePSF
        bbarSizer.Add(self.ButtonSavePSF , 1, flag=wx.ALL|wx.EXPAND, border=3)
        bbarSizer.Add(self.ButtonDisplayPSF , 1, flag=wx.ALL|wx.EXPAND, border=3)
        bbarSizer.Add(self.ButtonDisplayProf , 1, flag=wx.ALL|wx.EXPAND, border=3)
        bbarSizer.Add(self.ButtonOptions , 1, flag=wx.ALL|wx.EXPAND, border=3)
        bbarSizer.Add(self.ButtonQuit , 1, flag=wx.ALL|wx.EXPAND, border=3)
        bbar.SetSizerAndFit(bbarSizer)

        self.Bind(wx.EVT_BUTTON, self.ev_calcPSF, self.ButtonCompute)
        self.Bind(wx.EVT_BUTTON, self.ev_displayPSF, self.ButtonDisplayPSF)
        self.Bind(wx.EVT_BUTTON, self.ev_displayProfiles, self.ButtonDisplayProf)
        self.Bind(wx.EVT_BUTTON, self.ev_SaveAs, self.ButtonSavePSF)
        self.Bind(wx.EVT_BUTTON, self.ev_options, self.ButtonOptions)
        self.Bind(wx.EVT_BUTTON, self.OnClose, self.ButtonQuit)


        #==== Add items into top sizer of window
        topSizer.Add(sb1Sizer, 5, flag=wx.EXPAND|wx.ALL, border=10)
        topSizer.Add(nb, 10, flag=wx.EXPAND|wx.ALL, border=6)
        topSizer.Add(sb3Sizer, 9, flag=wx.EXPAND|wx.ALL, border=10)
        topSizer.Add(bbar, 0, flag=wx.EXPAND|wx.ALL, border=10)

        top_panel.SetSizerAndFit(topSizer)
        topSizer.Fit(self)
        self.SetSizeHints(self.GetSize().x,self.GetSize().y,-1,-1 ); #can get bigger but not smaller



        self.Bind(wx.EVT_CLOSE, self.OnClose)

        self.Center() 


    def _refresh_window(self):
        """ Force the window to refresh, and optionally to show itself if hidden (for recent matplotlibs)"""
        plt.draw()
        from distutils.version import StrictVersion
        if StrictVersion(matplotlib.__version__) >= StrictVersion('1.1'):
            plt.show(block=False)

    def log(self, messagestring):
        """ Display an informative string in both the log and the window status bar"""
        _log.info(messagestring)
        self.sb.SetStatusText(messagestring)
        self.Refresh()
        self.Update()
        wx.Yield()


    #--- Event Handlers for user-generated events
    def OnClose(self, event):
        dlg = wx.MessageDialog(self,
            "Do you really want to exit WebbPSF?",
            "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            if self.log_window is not None: self.log_window.Destroy()
            self.Destroy()

    def ev_SaveAs(self, event):
        "Event handler for Save As of output PSFs"
        dlg = wx.FileDialog(self, 'Choose Output Filename to Save PSF', defaultDir = os.getcwd(),
                defaultFile ='PSF_%s_%s.fits' %(self.iname, self.filter), 
                wildcard='*.fits',
                style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            filename = os.path.abspath(path)
            self.PSF_HDUlist.writeto(filename) 
            self.log("Saved to %s" % filename)
        else:
            self.log("User cancelled save.")

    def ev_SaveProfiles(self,event):
        """ Event handler for Save PSF Profiles"""
        raise NotImplementedError("Not implemented yet")

    def ev_Preferences(self,event):
        """ Event handler for WebbPSF overall preferences"""
        dlg = WebbPSFPreferencesDialog(self)
        results = dlg.ShowModal()
        dlg.Destroy()


    def ev_options(self, event):
        import copy
        import poppy
        oldoptions = copy.copy(self.advanced_options)

        dlg = WebbPSFOptionsDialog(self, input_options = self.advanced_options)
        results = dlg.ShowModal()

        iname = self.widgets['tabset'].GetCurrentPage().WebbPSFInstrumentName

        if dlg.results is not None: # none means the user hit 'cancel'
            self.advanced_options = dlg.results

            if self.advanced_options['monochromatic']:
                self.widgets["nlambda_pre_label"].SetLabel("Wavelength:    ")
                self.widgets["nlambda_post_label"].SetLabel(" [microns]")
            else:
                self.widgets["nlambda_pre_label"].SetLabel("# of wavelengths:")
                self.widgets["nlambda_post_label"].SetLabel("Leave blank for autoselect")

            if self.advanced_options['fov_in_arcsec']:
                self.widgets["FOV_post_label"].SetLabel("arcsec/side")
                if not oldoptions['fov_in_arcsec']: 
                    # we have to convert pixels to arcsec
                    self._setFOV( self._getFOV() * self.instrument[iname].pixelscale) 
            else:
                self.widgets["FOV_post_label"].SetLabel("pixels/side")
                if oldoptions['fov_in_arcsec']: 
                    # we have to convert arcsec to pixels
                    newfov = self._getFOV() / self.instrument[iname].pixelscale
                    if hasattr(newfov, '__len__') and len(newfov) > 1: 
                        newfov = np.asarray(newfov, dtype=int)
                    else:
                        newfov = int(newfov)
                    self._setFOV( newfov)
            poppy.conf.use_multiprocessing = self.advanced_options['parallelization']
 
        dlg.Destroy()


    def ev_showDocs(self,event):
        """ Event handler for show documentation"""
        raise NotImplementedError("Not implemented yet")

    def ev_AboutBox(self, event):
        about = AboutBox()
        about.ShowModal()


    def ev_plotspectrum(self, event):
        "Event handler for Plot Spectrum "
        self._updateFromGUI()
        self.log("Now calculating spectrum model...")

        self.log("Spectral type is "+self.sptype)
        self.log("Selected instrument tab is "+self.iname)
        self.log("Selected instrument filter is "+self.filter)


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
        #try:
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
        self.log("Spectrum displayed")

    def ev_calcPSF(self, event):
        "Event handler for PSF Calculations"
        self._updateFromGUI()
        self.log("Starting PSF calculation...")


        self._refresh_window() # pre-display window for calculation updates as it progresses...
        for w in ['Display PSF', 'Display profiles', 'Save PSF As...']:
           self.widgets[w].Enable(False)

        self.calcthread = PSFCalcThread()
        self.calcthread.runPSFCalc(self.inst, self) 
#        self.PSF_HDUlist = self.inst.calcPSF(source=source, 
#                detector_oversample= self.detector_oversampling,
#                fft_oversample=self.fft_oversampling,
#                fov_arcsec = self.FOV,  nlambda = self.nlambda, display=True)
        #self.PSF_HDUlist.display()


    def OnCalcDone(self, event):
        """ Called when worker thread returns results..."""
        if event.data is None:
            self.log("No result from calculation!")
            self.PSF_HDUlist = None
        else:
            self.PSF_HDUlist = event.data
            for w in ['Display PSF', 'Display profiles', 'Save PSF As...']:
               self.widgets[w].Enable(True)

        self._refresh_window()
        self.log("PSF calculation complete")
        # In either event, the worker is done
        self.calcthread = None


    def ev_displayPSF(self,event):
        "Event handler for Displaying the PSF"
        #self._updateFromGUI()
        #if self.PSF_HDUlist is not None:
        plt.clf()
        poppy.display_PSF(self.PSF_HDUlist, vmin = self.advanced_options['psf_vmin'], vmax = self.advanced_options['psf_vmax'], 
                scale = self.advanced_options['psf_scale'], cmap= self.advanced_options['psf_cmap'], normalize=self.advanced_options['psf_normalize'])
        self._refresh_window()
        self.log("PSF redisplayed")

    def ev_displayProfiles(self,event):
        "Event handler for Displaying the PSF"
        #self._updateFromGUI()
        poppy.display_profiles(self.PSF_HDUlist)        
        self._refresh_window()
        self.log("Radial profiles displayed.")

    def ev_displayOptics(self,event):
        "Event handler for Displaying the optical system"
        self._updateFromGUI()
        self.log("Selected OPD is "+str(self.opd_name))

        plt.clf()
        self.inst.display()
        self._refresh_window()
        self.log("Optical system elements displayed.")

    def ev_displayOPD(self,event):
        import poppy.utils

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
            poppy.utils.imshow_with_mouseover(masked_opd, cmap=cmap, interpolation='nearest', vmin=-0.5, vmax=0.5)
            plt.title("OPD from %s, #%d" %( os.path.basename(self.opd_name), self.opd_i))
            cb = plt.colorbar(orientation='vertical')
            cb.set_label('microns')

            f = plt.gcf()
            plt.text(0.4, 0.02, "OPD WFE = %6.2f nm RMS" % (masked_opd.std()*1000.), transform=f.transFigure)
        self.log("Optical Path Difference (OPD) now displayed.")

        self._refresh_window()

    def ev_launch_ITM_dialog(self, event):
        tkMessageBox.showwarning( message="ITM dialog box not yet implemented", title="Can't Display") 

    def ev_update_OPD_labels(self, event):
        "Update the descriptive text for all OPD files"
        for iname in self.instrument:
            self.ev_update_OPD_label(self.widgets[iname+"_opd"], self.widgets[iname+"_opd_label"], iname)

    def ev_update_OPD_label(self, widget_combobox, widget_label, iname):
        "Update the descriptive text displayed about one OPD file"
        from wx.lib.wordwrap import wordwrap
        showitm=False # default is do not show
        filename = os.path.join( self.instrument[iname]._datapath, 'OPD',  widget_combobox.GetValue() )

        if filename.endswith(".fits"):
            try:
                header_summary = fits.getheader(filename)['SUMMARY']
            except KeyError:
                header_summary = 'could not obtain a summary from FITS header'
            #self.widgets[iname+"_opd_i"]['state'] = 'readonly'
        else:  # Special options for non-FITS file inputs
            #self.widgets[iname+"_opd_i"]['state'] = 'disabled'
            val = widget_combobox.GetValue()
            if 'Zero' in val:
                header_summary = " 0 nm RMS"
            elif 'ITM' in val and self._enable_opdserver:
                header_summary= "Get OPD from ITM Server"
                showitm=True
            elif 'ITM' in val and not self._enable_opdserver:
                header_summary = "ITM Server is not running or otherwise unavailable."
            else: # other??
                header_summary = "   "

        widget_label.SetLabel( header_summary  )
        widget_label.Wrap(250)
        #widget_label.SetLabel( wordwrap( header_summary, 40, wx.ClientDC(self))  )


        #if showitm:
        #    self.widgets[iname+"_itm_coords"].grid() # re-show ITM options
        #else:
        #    self.widgets[iname+"_itm_coords"].grid_remove()  # hide ITM options

    def _updateFromGUI(self):
        """ Update the object's state with all the settings from the GUI
        """
        # get GUI values
        if _HAS_PYSYNPHOT:
            self.sptype = self.widgets['SpType'].GetValue()
        self.iname = self.widgets['tabset'].GetCurrentPage().WebbPSFInstrumentName


        if self.advanced_options['monochromatic']:
            # monochromatic mode
            self.nlambda=1
            self.monochromatic_wavelength = float(self.widgets['nlambda'].GetValue()) * 1e-6 # note that the GUI is in microns, so convert to meters here.

            if self.monochromatic_wavelength < 0.6e-6 or self.monochromatic_wavelength > 30e-6:
                _log.error("Invalid wavelength. Please enter a value between 0.6 - 30 microns")
                raise ValueError('Invalid wavelength outside of range 0.6 - 30 microns: {0:f}'.format(self.monochromatic_wavelength))
            try:
                pass
            except:
                _log.error("Could not obtain a wavelength. Please check the value of the wavelength field and try again.")
                raise ValueError('Could not obtain a wavelength from string "{0}"'.format(self.widgets['nlambda'].GetValue()))
        else:
            # normal broadband mode
            self.monochromatic_wavelength = None
            try:
                self.nlambda= int(self.widgets['nlambda'].GetValue())
            except:
                self.nlambda = None # invoke autoselect for nlambda
                



        self.FOV= self._getFOV() #float(self.widgets['FOV'].GetValue())
        self.fft_oversampling= int(self.widgets['fft_oversampling'].GetValue())
        self.detector_oversampling= int(self.widgets['detector_oversampling'].GetValue())

        self.output_type = self.widgets['output_format'].GetValue()

        options = {}
        options['rebin'] = not (self.output_type == 'Oversampled PSF only')  #was downsample, which seems wrong?
        options['mock_dms'] = (self.output_type == 'Mock full image from JWST DMS')
        jitterchoice = self.widgets['jitter'].GetValue()
        if jitterchoice == 'Just use OPDs':
            options['jitter'] = None
        elif jitterchoice == 'Gaussian jitter with 7 mas rms':
            options['jitter'] = 'gaussian'
            options['jitter_sigma'] = 0.007
        elif jitterchoice == 'Gaussian jitter with 30 mas rms':
            options['jitter'] = 'gaussian'
            options['jitter_sigma'] = 0.030
        else: 
            _log.error("Unknown value for jitter selection: "+jitterchoice)



        # and get the values that may have previously been set by the 'advanced options' dialog
        if self.advanced_options is not None:
            for a in self.advanced_options:
                options[a] = self.advanced_options[a]


        # configure the relevant instrument object
        self.inst = self.instrument[self.iname]
        self.filter = self.widgets[self.iname+"_filter"].GetValue() # save for use in default filenames, etc.
        self.inst.filter = self.filter
        _log.info("Selected filter: "+self.filter)

        self.opd_name = self.widgets[self.iname+"_opd"].GetValue()
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
            self.opd_name= self.widgets[self.iname+"_opd"].GetValue()
            self.opd_i= int(self.widgets[self.iname+"_opd_i"].GetValue())
            self.inst.pupilopd = (self.inst._datapath+os.sep+"OPD"+os.sep+self.opd_name,self.opd_i)  #filename, slice

        self.log("Selected OPD is "+str(self.opd_name))


        if self.iname+"_coron" in self.widgets:
            self.inst.image_mask = self.widgets[self.iname+"_coron"].GetValue()
            self.inst.pupil_mask = self.widgets[self.iname+"_pupil"].GetValue()
            # TODO read in mis-registration options here.


            options['source_offset_r'] = float(self.widgets["source_off_r"].GetValue())
            options['source_offset_theta'] = float(self.widgets["source_off_theta"].GetValue())
            options['pupil_shift_x'] = float(self.widgets[self.iname+"_pupilshift_x"].GetValue())/100. # convert from percent to fraction
            options['pupil_shift_y'] = float(self.widgets[self.iname+"_pupilshift_y"].GetValue())/100. # convert from percent to fraction

        self.inst.options = options


    def _getFOV(self):
        """ Get field of view, either as a scalar number for square or a 
        2-element ndarray for rectangular 
        
        Note that if it is a 2-element paid, we flip the order of the elements. 
        This is to facilitate a more intuitive "x,y" ordering for the user interface.
        """
        fovstr = self.widgets['FOV'].GetValue()

        try:
            if ',' in fovstr:
                parts = fovstr.split(',')
                return np.asarray( parts[0:2], dtype=float)[::-1]
            else:
                return float(self.widgets['FOV'].GetValue())
        except:
            _log.error("Invalid entry in FOV field. Please check it and try again")
            raise ValueError("Invalid entry in FOV field. Please check it and try again")


    def _setFOV(self, newvalue):
        """ Get field of view, either as a scalar number for square or a 
        2-element ndarray for rectangular 

        Note that if it is a 2-element paid, we flip the order of the elements. 
        This is to facilitate a more intuitive "x,y" ordering for the user interface.
        """
        if hasattr(newvalue, '__iter__'):
            newstring = ",".join( (str(newvalue[1]), str(newvalue[0]) ))
        else:
            newstring = str(newvalue)

        self.widgets['FOV'].SetValue(newstring)

#-------------------------------------------------------------------------
# Class to run the actual PSF calculation in a background thread, to keep the
# GUI still responsive
# This code based on examples at http://wiki.wxpython.org/LongRunningTasks         
# and http://www.blog.pythonlibrary.org/2010/05/22/wxpython-and-threads/

class PSFCalcThread(Thread):
    def __init__(self):
        """Init Worker Thread Class."""
        Thread.__init__(self)
        self.start()    # start the thread

    def runPSFCalc(self,  instrument, masterapp):

        if _HAS_PYSYNPHOT:
            source = poppy.specFromSpectralType(masterapp.sptype)
        else:
            source=None # generic flat spectrum
        print("starting calc in thread")

        if instrument.options['fov_in_arcsec']:
            fov_arcsec = masterapp.FOV
            fov_pixels = None
        else:
            fov_arcsec = None
            fov_pixels = masterapp.FOV
 
        PSF_HDUlist = instrument.calcPSF(source=source, 
                detector_oversample = masterapp.detector_oversampling,
                fft_oversample = masterapp.fft_oversampling,
                fov_arcsec = fov_arcsec, fov_pixels=fov_pixels,  
                nlambda = masterapp.nlambda, 
                monochromatic=masterapp.monochromatic_wavelength, 
                display = True)

        wx.PostEvent(masterapp, ResultEvent(PSF_HDUlist)) # send results back to master thread


# Define notification event for thread completion
EVT_RESULT_ID = wx.NewId()
def EVT_RESULT(win, func):
    """Define Result Event."""
    win.Connect(-1, -1, EVT_RESULT_ID, func)

class ResultEvent(wx.PyEvent):
    """Simple event to carry arbitrary result data."""
    def __init__(self, data):
        """Init Result Event."""
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_RESULT_ID)
        self.data = data



#-------------------------------------------------------------------------

class WebbPSFMenuBar(wx.MenuBar):
    def __init__(self,parent):
        wx.MenuBar.__init__(self)
        item_keys =['save_psf','save_profile', 'documentation', 'preferences', 'calc_options', 'calcPSF', 
                'display_spectrum', 'display_optics','display_opd', 'display_psf', 'display_profiles']

        self.ids = {}
        for key in item_keys:
            self.ids[key]=wx.NewId()

        # File menu
        filemenu = wx.Menu()
        self.SavePSF = filemenu.Append(self.ids['save_psf'], '&Save PSF as...\tCtrl+Shift+S')
        self.SaveProfile =filemenu.Append(self.ids['save_profile'], 'Save &profile data as...\tCtrl+Shift+P')
        filemenu.AppendSeparator()
        self.Preferences=filemenu.Append(self.ids['preferences'], 'Preferences...\tCtrl+,')
        filemenu.AppendSeparator()
        self.Exit = filemenu.Append(wx.ID_EXIT, 'E&xit', 'Exit this program')

        # these start out disabled since no PSF calculated yet:
        self.SavePSF.Enable(False)
        self.SaveProfile.Enable(False)

        # Edit menu
        editmenu = wx.Menu()
        editmenu.Append(wx.ID_CUT, 'Cut\tCtrl-X')
        editmenu.Append(wx.ID_COPY, 'Copy\tCtrl-C')
        editmenu.Append(wx.ID_PASTE, 'Paste\tCtrl-V')

        #Calculation Menu
        calcmenu = wx.Menu()
        calcmenu.Append(self.ids['calcPSF'], 'Compute PSF')
        self.Preferences=calcmenu.Append(self.ids['calc_options'], 'More Options...')
        calcmenu.AppendSeparator()
        calcmenu.Append(self.ids['display_spectrum'], 'Display Spectrum')
        calcmenu.Append(self.ids['display_optics'], 'Display Optics')
        calcmenu.Append(self.ids['display_opd'], 'Display OPD')
        calcmenu.Append(self.ids['display_psf'], 'Display PSF')
        calcmenu.Append(self.ids['display_profiles'], 'Display PSF Profiles')

        #Help Menu
        helpmenu = wx.Menu()
        self.Docs = helpmenu.Append(self.ids['documentation'], 'WebbPSF Documentation\tCtrl+d')
        self.About = helpmenu.Append(wx.ID_ABOUT, '&About WebbPSF', 'About this program')

        self.Append(filemenu, '&File')
        self.Append(editmenu, '&Edit')
        self.Append(calcmenu, '&Calculation')
        self.Append(helpmenu, '&Help')

#-------------------------------------------------------------------------
class WebbPSFDialog(wx.Dialog):
    """ Generic dialog box for WebbPSF 

    TODO: investigate wx.Validator to validate the text input fields
    """
    def __init__(self, parent=None, id=-1, title="Dialog", **kwargs):
        wx.Dialog.__init__(self,parent,id=id,title=title, **kwargs)

        self.parent=parent

        self.results = None # in case we cancel this gets returned
        self.widgets = {}
        self.values = {}


        #self._createWidgets()
 
    def _add_labeled_dropdown(self, name, parent, parentsizer, label="Entry:", choices=None, default=0, width=5, position=(0,0), columnspan=1, **kwargs):
        """convenient wrapper for adding a Combobox

        columnspan sets the span for the combobox itself
        """

        mylabel = wx.StaticText(parent, -1,label=label)
        parentsizer.Add( mylabel, position,  (1,1), wx.EXPAND)

        if choices is None:
            try:
                choices = self.values[name]
            except:
                choices = ['Option A','Option B']

        if isinstance(default,int):
            value = choices[default]
        else:
            value=default

        mycombo = wx.ComboBox(parent, -1,value=value, choices=choices, style=wx.CB_DROPDOWN|wx.CB_READONLY)
        parentsizer.Add( mycombo, (position[0],position[1]+1),  (1,columnspan), wx.EXPAND)

        self.widgets[name] = mycombo
 

    def _add_labeled_entry(self, name, parent, parentsizer, label="Entry:", value=None, format="%.2g", 
            width=5, position=(0,0), postlabel=None, **kwargs):
        "convenient wrapper for adding an Entry"
        mylabel = wx.StaticText(parent, -1,label=label)
        parentsizer.Add( mylabel, position,  (1,1), wx.EXPAND)

        if value is None:
            try:
                value=format % self.input_options[name]
            except:
                value=""
        else:
            try: 
                value = format % value
            except:
                pass

        mytext = wx.TextCtrl(parent, -1,value=value)
        parentsizer.Add( mytext, (position[0],position[1]+1),  (1,1), wx.EXPAND)

        self.widgets[name] = mytext

        if postlabel is not None:
            mylabel2 = wx.StaticText(parent, -1,label=postlabel)
            parentsizer.Add( mylabel2, (position[0],position[1]+2),  (1,1), wx.EXPAND)


    def _createWidgets(self):
        pass
        # subclass me


    def OnButtonOK(self, event):
        pass


class WebbPSFOptionsDialog(WebbPSFDialog):
    """ Dialog box for WebbPSF options 

    TODO: investigate wx.Validator to validate the text input fields
    """
    def __init__(self, parent=None, id=-1, title="WebbPSF Options", 
            input_options=_default_options()): 
        WebbPSFDialog.__init__(self, parent,id=id,title=title)
        self.input_options = input_options

        colortables = [
         ('Jet (blue to red)',matplotlib.cm.jet),
         ('Gray', matplotlib.cm.gray),
         ('Heat (black-red-yellow)', matplotlib.cm.gist_heat),
         ('Copper (black to tan)',matplotlib.cm.copper),
         ('Stern',matplotlib.cm.gist_stern),
         ('Prism (repeating rainbow)', matplotlib.cm.prism)]

        try:
            import collections
            self.colortables = collections.OrderedDict(colortables)
        except:
            self.colortables = dict(colortables)

        self.values['force_coron'] = ['regular propagation (MFT)', 'full coronagraphic propagation (FFT/SAM)']
        self.values['no_sam'] = ['semi-analytic method if possible', 'basic FFT method always']
        self.values['monochromatic'] = ['Broadband', 'Monochromatic']
        self.values['parallelization'] = ['Sequential', 'Parallelized']
        self.values['fov_in_arcsec'] = ['Arcseconds', 'Pixels']


        self._createWidgets()
 

    def _createWidgets(self):

        topSizer = wx.BoxSizer(wx.VERTICAL)

        topPanel = wx.Panel(self)
        topPanelSizer = wx.FlexGridSizer(rows=3,cols=1, hgap=5, vgap=10)
        #topPanelSizer = wx.BoxSizer(wx.VERTICAL)

        panel1 = wx.Panel(topPanel, style=wx.SIMPLE_BORDER|wx.EXPAND)
        sizer = wx.GridBagSizer()

        txt = wx.StaticText(panel1,-1,label="Propagation Calculation Options")
        sizer.Add(txt,(0,0),(1,3),wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)


        r=1
        self._add_labeled_dropdown("monochromatic", panel1,sizer, 
                label='    Broadband or monochromatic? ', 
                default = 1 if self.input_options['monochromatic'] else 0, position=(r,0))
        r+=1
        self._add_labeled_dropdown("parallelization", panel1,sizer, 
                label='    Parallelize Wavelengths? ', 
                default = 1 if self.input_options['parallelization'] else 0, position=(r,0))
        r+=1
        self._add_labeled_dropdown("force_coron", panel1,sizer, 
                label='    Direct imaging calculations use: ', 
                default = 1 if self.input_options['force_coron'] else 0, position=(r,0))

        r+=1
        self._add_labeled_dropdown("no_sam", panel1,sizer, 
                label='    Coronagraphic calculations use', 
                default= 1 if self.input_options['no_sam'] else 0, position=(r,0))
        r+=1
        self._add_labeled_dropdown("parity", panel1,sizer, 
                label='    Output pixel grid parity is', 
                choices=['odd', 'even', 'either'], default=self.input_options['parity'], position=(r,0))
        r+=1
        self._add_labeled_dropdown("fov_in_arcsec", panel1,sizer, 
                label='    Specify field of view in: ', 
                default = 0 if self.input_options['fov_in_arcsec'] else 1, position=(r,0))

 
        #sizer.AddGrowableCol(0)
        panel1.SetSizerAndFit(sizer)

        panel2 = wx.Panel(topPanel, style=wx.SIMPLE_BORDER|wx.EXPAND)
        sizer = wx.GridBagSizer()

        r=0
        txt = wx.StaticText(panel2,-1,label="PSF Display Options")
        sizer.Add(txt,(r,0),(1,3),wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)


        r+=1
        self._add_labeled_dropdown("psf_scale", panel2,sizer, label='    Display scale:', choices=['log','linear'],
                default=self.input_options['psf_scale'], position=(r,0))

        r+=1
        self._add_labeled_entry("psf_vmin", panel2,sizer, label='    Min scale value:', 
            format="%.2g", width=7, 
            position=(r,0))
        r+=1
        self._add_labeled_entry("psf_vmax", panel2,sizer, label='    Max scale value:', 
            format="%.2g", width=7, 
            position=(r,0))
        r+=1
        self._add_labeled_dropdown("psf_normalize", panel2,sizer,
                label='    Normalize PSF to:', choices=['Total', 'Peak'], 
                default=self.input_options['psf_normalize'], position=(r,0))
        r+=1
        self._add_labeled_dropdown("psf_cmap", panel2,sizer, label='    Color table:', 
                choices=[a for a in self.colortables],  
                default=self.input_options['psf_cmap_str'], 
                position=(r,0))
        panel2.SetSizerAndFit(sizer)


        bbar = self.CreateStdDialogButtonSizer( wx.OK | wx.CANCEL)
        self.Bind(wx.EVT_BUTTON, self.OnButtonOK, id=wx.ID_OK)
        #self.Bind(wx.EVT_BUTTON, self.OnButtonCancel, id=wx.ID_CANCEL)
 

        topPanelSizer.Add(panel1, 1,wx.EXPAND)
        topPanelSizer.Add(panel2, 1, wx.EXPAND)
        #topPanelSizer.Add(panel3,1, wx.EXPAND)
        topPanelSizer.Add(bbar,1, wx.EXPAND)
        #topPanel.AddGrowableCol(0)
        topPanel.SetSizerAndFit(topPanelSizer)
 
        topSizer.Add(topPanel, 1, flag=wx.EXPAND|wx.ALL, border=10)
        self.SetSizerAndFit(topSizer)
        self.Show(True)
    #def OnButtonCancel(self, event):
        #print("User pressed Cancel")
        #self.Close()
        #self.Destroy()

    def OnButtonOK(self, event):
        print("User pressed OK")
        try:
            results = {}
            results['force_coron'] =    self.widgets['force_coron'].GetValue() == 'full coronagraphic propagation (FFT/SAM)'
            results['no_sam'] =         self.widgets['no_sam'].GetValue() == 'basic FFT method always'
            results['monochromatic'] =  self.widgets['monochromatic'].GetValue() == 'Monochromatic'
            results['parallelization'] =  self.widgets['parallelization'].GetValue() == 'Parallelized'
            results['fov_in_arcsec'] =  self.widgets['fov_in_arcsec'].GetValue() == 'Arcseconds'
            results['parity'] =         self.widgets['parity'].GetValue() 
            results['psf_scale'] =      self.widgets['psf_scale'].GetValue() 
            results['psf_vmax'] = float(self.widgets['psf_vmax'].GetValue())
            results['psf_vmin'] = float(self.widgets['psf_vmin'].GetValue())
            results['psf_cmap_str'] =   self.widgets['psf_cmap'].GetValue()
            results['psf_cmap'] =       self.colortables[self.widgets['psf_cmap'].GetValue() ]
            results['psf_normalize'] =  self.widgets['psf_normalize'].GetValue()


            print(results)
            self.results = results # for access from calling routine

            self.Close()
            #self.Destroy() # return... If called as a modal dialog, should ShowModal and Destroy from calling routine?
        except:
            _log.error("Invalid entries in one or more fields. Please check values and re-enter!")

#-------------------------------------------------------------------------

class WebbPSFPreferencesDialog(WebbPSFDialog):
    """ Dialog box for WebbPSF options 

    TODO: investigate wx.Validator to validate the text input fields
    """
    def __init__(self, parent=None, id=-1, title="WebbPSF Preferences"): 
        WebbPSFDialog.__init__(self, parent,id=id,title=title, size=(800,400))

        self._createWidgets()

    def _createWidgets(self):

        topSizer = wx.BoxSizer(wx.VERTICAL)

        topPanel = wx.Panel(self)
        topPanelSizer = wx.FlexGridSizer(rows=4,cols=1, hgap=5, vgap=10)

        panel1 = wx.Panel(topPanel, style=wx.SIMPLE_BORDER|wx.EXPAND)
        sizer = wx.GridBagSizer()

        txt = wx.StaticText(panel1,-1,label="Directory Paths")
        sizer.Add(txt,(0,0),(1,3),wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)

        r=1
        self._add_labeled_entry("WEBBPSF_PATH", panel1,sizer, label='    WebbPSF Data Path:', 
            value = str(conf.WEBBPSF_PATH), 
            format="%50s", position=(r,0))

        self.ButtonBrowseDir = wx.Button(panel1, label='Browse...')
        sizer.Add(self.ButtonBrowseDir, (r,4), (1,1), wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.ev_browseDataDirectory, self.ButtonBrowseDir)
        panel1.SetSizerAndFit(sizer)



        panel2 = wx.Panel(topPanel, style=wx.SIMPLE_BORDER|wx.EXPAND)
        sizer = wx.GridBagSizer()

        txt = wx.StaticText(panel2,-1,label="Calculation Defaults")
        sizer.Add(txt,(0,0),(1,3),wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)

        r=1
        self._add_labeled_entry("default_oversampling", panel2,sizer, label='    Default Oversampling:', 
            value = conf.default_oversampling, 
            format="%.2g", width=7, position=(r,0))
        r+=1
        self._add_labeled_entry("default_fov_arcsec", panel2,sizer, label='    Default FOV [arcsec]:', 
            value = conf.default_fov_arcsec, 
            format="%.2g", width=7, position=(r,0))

        # update conf.default_oversampling , default_output_mode, default_fov_arcsec, WEBBPSF_PATH, 
        # use_multiprocessing= True/False 
        # n_processes= #
        # use_fftw = True/False

        panel2.SetSizerAndFit(sizer)

        panel3 = wx.Panel(topPanel, style=wx.SIMPLE_BORDER|wx.EXPAND)
        sizer = wx.GridBagSizer()

        r=0
        txt = wx.StaticText(panel3,-1,label="Fourier Transform Options")
        sizer.Add(txt,(r,0),(1,3),wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)

        r+=1
        self._add_labeled_dropdown("use_fftw", panel3,sizer, label='    Use FFTW for FFTs:', choices=['True','False'],
                default="True" if poppy.conf.use_fftw else "False", position=(r,0))
        r+=1
        self._add_labeled_dropdown("use_multiprocessing", panel3,sizer, label='    Use Multiprocessing for DFTs:', choices=['True','False'],
                default="True" if poppy.conf.use_multiprocessing else "False", position=(r,0))

        panel3.SetSizerAndFit(sizer)
#
#
#        r+=1
#        self._add_labeled_entry("psf_vmin", panel2,sizer, label='    Min scale value:', 
#            format="%.2g", width=7, position=(r,0))
#        r+=1
#        self._add_labeled_entry("psf_vmax", panel2,sizer, label='    Max scale value:', 
#            format="%.2g", width=7, 
#            position=(r,0))
#        r+=1
#        self._add_labeled_dropdown("psf_normalize", panel2,sizer,
#                label='    Normalize PSF to:', choices=['Total', 'Peak'], 
#                default=self.input_options['psf_normalize'], position=(r,0))
#        r+=1
#        self._add_labeled_dropdown("psf_cmap", panel2,sizer, label='    Color table:', 
#                choices=[a for a in self.colortables.keys()],  
#                default=self.input_options['psf_cmap_str'], 
#                position=(r,0))
#        panel2.SetSizerAndFit(sizer)
#

        bbar = self.CreateStdDialogButtonSizer( wx.OK | wx.CANCEL)
        self.Bind(wx.EVT_BUTTON, self.OnButtonOK, id=wx.ID_OK)

        topPanelSizer.Add(panel1, 1,wx.EXPAND)
        topPanelSizer.Add(panel2, 1, wx.EXPAND)
        topPanelSizer.Add(panel3, 1, wx.EXPAND)
        topPanelSizer.Add(bbar,1, wx.EXPAND)
        topPanel.SetSizerAndFit(topPanelSizer)
 
        topSizer.Add(topPanel, 1, flag=wx.EXPAND|wx.ALL, border=10)
        self.SetSize((800,300))
        self.SetSizerAndFit(topSizer)

        self.Show(True)

  
    def ev_browseDataDirectory(self, event):
        dlg = wx.DirDialog(self, 'Choose WebbPSF Data Directory', defaultPath=os.getcwd())

        if dlg.ShowModal() == wx.ID_OK:
            print("Directory Selected: "+dlg.GetPath())
        else:
            print("Dialog cancelled")

    def OnButtonOK(self, event):
        print("User pressed OK")
        try:
            conf.default_oversampling = int(self.widgets['default_oversampling'].GetValue())
            conf.default_fov_arcsec = float(self.widgets['default_fov_arcsec'].GetValue())

#            results['force_coron'] =    self.widgets['force_coron'].GetValue() == 'full coronagraphic propagation (FFT/SAM)'
#            results['no_sam'] =         self.widgets['no_sam'].GetValue() == 'basic FFT method always'
#            results['monochromatic'] =  self.widgets['monochromatic'].GetValue() == 'Monochromatic'
#            results['fov_in_arcsec'] =  self.widgets['fov_in_arcsec'].GetValue() == 'Arcseconds'
#            results['parity'] =         self.widgets['parity'].GetValue() 
#            results['psf_scale'] =      self.widgets['psf_scale'].GetValue() 
#            results['psf_vmax'] = float(self.widgets['psf_vmax'].GetValue())
#            results['psf_vmin'] = float(self.widgets['psf_vmin'].GetValue())
#            results['psf_cmap_str'] =   self.widgets['psf_cmap'].GetValue()
#            results['psf_cmap'] =       self.colortables[self.widgets['psf_cmap'].GetValue() ]
#            results['psf_normalize'] =  self.widgets['psf_normalize'].GetValue()
#
#
#            print(results)
#            self.results = results # for access from calling routine

            utils.save_config()
            self.Close()
            #self.Destroy() # return... If called as a modal dialog, should ShowModal and Destroy from calling routine?
        except:
            _log.error("Invalid entries in one or more fields. Please check values and re-enter!")


#-------------------------------------------------------------------------
#
# Classes for displaying log messages in a separate window

class WxLog(logging.Handler):
    def __init__(self, ctrl):
        logging.Handler.__init__(self)
        self.ctrl = ctrl
    def emit(self, record):
        self.ctrl.AppendText(self.format(record)+"\n")
        self.ctrl.Update()

class LogFrame(wx.Frame):
    def __init__(self, parent=None, id=-1, pos=None, size=(600,300)):
        wx.Frame.__init__(self, parent, id=id, title="WebbPSF tasks log", size=size,pos=pos)
        #self.level = 4
        log = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(log, 1, wx.EXPAND)
        self.SetSizer(sizer)


        hdlr = WxLog(log)
        hdlr.setFormatter(logging.Formatter('%(name)-8s %(levelname)-8s: %(message)s'))

        self.hdlr = hdlr # save for later deletion

        logging.getLogger('webbpsf').addHandler(hdlr)
        logging.getLogger('poppy').addHandler(hdlr)
        
    def __del__(self):
        logging.getLogger('webbpsf').removeHandler(self.hdlr)
        logging.getLogger('poppy').removeHandler(self.hdlr)

#-------------------------------------------------------------------------

class HtmlWindow(wx.html.HtmlWindow):
    def __init__(self, parent, id, size=(600,400)):
        wx.html.HtmlWindow.__init__(self,parent, id, size=size)
        if "gtk2" in wx.PlatformInfo:
            self.SetStandardFonts()

    def OnLinkClicked(self, link):
        wx.LaunchDefaultBrowser(link.GetHref())

class AboutBox(wx.Dialog):
    def __init__(self):
        wx.Dialog.__init__(self, None, wx.ID_ANY, "About WebbPSF",
            style=wx.DEFAULT_DIALOG_STYLE|wx.THICK_FRAME|wx.RESIZE_BORDER|
                wx.TAB_TRAVERSAL)


        aboutText = """<p>This is the <b>WebbPSF Point Spread Function Simulator</b> for the James Webb Space Telescope (JWST),
<p>
<center>
Version %(webbpsf)s <P>
See <a href="http://www.stsci.edu/jwst/software/webbpsf">the WebbPSF home page</a>
</center>
<p>
(c) 2010-2013 by Marshall Perrin, <a href="mailto:mperrin@stsci.edu">mperrin@stsci.edu</a>.
<br>
With contributions from: Anand Sivaramakrishnan, Remi Soummer, &amp; Klaus Pontoppidan
<p>
WebbPSF is running with the following software versions:
<ul>
<li><b>Python</b>: %(python)s 
<li><b>numpy</b>: %(numpy)s
<li><b>matplotlib</b>: %(matplotlib)s
<li><b>astropy</b>: %(astropy)s
<li><b>wxPython</b>: %(wxpy)s  
<li><b>pysynphot</b>: %(pysynphot)s  
<li><b>pyFFTW</b>: %(pyfftw)s  
<li><b>PyFFTW3</b>: %(pyfftw3)s  
</ul>
</p>"""

        import sys, astropy
        hwin = HtmlWindow(self, wx.ID_ANY, size=(400,200))
        vers = {}
        from . import _version
        vers["webbpsf"] = _version.__version__
        try:
            import pysynphot
            vers['pysynphot'] = pysynphot.__version__
        except:
            vers['pysynphot'] = "Not Found"
        try:
            import pyfftw
            vers['pyfftw'] = "Present"
        except:
            vers['pyfftw'] = "Not Found"
        try:
            import fftw3
            vers['pyfftw3'] = "Present"
        except:
            vers['pyfftw3'] = "Not Found"
            
        vers["python"] = sys.version.split()[0]
        vers["wxpy"] = wx.VERSION_STRING
        vers['numpy'] = np.__version__
        vers['matplotlib'] = matplotlib.__version__
        #vers['atpy'] = atpy.__version__
        vers['astropy'] = astropy.__version__
        hwin.SetPage(aboutText % vers)
        btn = hwin.FindWindowById(wx.ID_OK)
        irep = hwin.GetInternalRepresentation()
        hwin.SetSize((irep.GetWidth()+25, irep.GetHeight()+10))
        self.SetClientSize(hwin.GetSize())
        self.CentreOnParent(wx.BOTH)
        self.SetFocus()


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




def wxgui(fignum=1, showlog=True):
    import poppy
    # enable log message printout
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')


    # GUI does not play well with multiprocessing, so avoid that.
    if poppy.conf.use_multiprocessing:
        _log.error('Multiprocessing is not compatible with the GUI right now. Falling back to single-threaded.')
        poppy.conf.use_multiprocessing = False 

    # start the GUI
    app = wx.App()
    app.SetAppName("WebbPSF")
    gui = WebbPSF_GUI()
    gui.Show()

    if showlog: 
        # start it immediately below the main window
        gpos = gui.GetScreenPosition()
        gsize = gui.GetSize()

        logwin = LogFrame(pos=(gpos[0], gpos[1]+gsize[1]), size=(gsize[0], 200))
        logwin.Show()

        gui.log_window = logwin
    #bout = AboutBox()
    #about.ShowModal()
    #gui = WebbPSFOptionsDialog()
    #plt.figure(fignum)
    #plt.show(block=False)
    app.MainLoop()



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    wxgui()



