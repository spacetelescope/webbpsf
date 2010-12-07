#!/usr/bin/env python
import os
import numpy as N
import pylab as P
import matplotlib
import pyfits
import ttk
import Tkinter as tk
import tkMessageBox
import tkFileDialog
#from Tkinter import N,E,S,W

try:
    __IPYTHON__
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    pass


import jwopt


class JWPSF_GUI(object):
    """ A GUI for the PSF Simulator 

    Documentation TBD!

    """
    def __init__(self):
        # init the object and subobjects
        self.instrument = {}
        self.widgets = {}
        self.vars = {}
        insts = ['NIRCam', 'NIRSpec','MIRI', 'TFI', 'FGS']
        self.stars = jwopt.kurucz_stars()
        for i in insts:
            self.instrument[i] = jwopt.Instrument(i)




        #---- create the GUIs
        self.root = tk.Tk()
        self.root.geometry('+50+50')
        self.root.title("JWPSF GUI 4alpha")

        frame = ttk.Frame(self.root)
        #frame = ttk.Frame(self.root, padx=10,pady=10)

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
                label = 'IFU' if iname !='TFI' else 'TFI'
                ttk.Label(fr2, text='   %s wavelen: '%label , state='disabled').grid(row=0, column=0)
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
        self.widgets['nlambda'].insert(0,'1')
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



        self.root.update()

    def quit(self):
        " Quit the GUI"
        if tkMessageBox.askyesno( message='Are you sure you want to quit?', icon='question', title='Install') :
            self.root.destroy()

    def ev_SaveAs(self):
        "Event handler for Save As of output PSFs"
        filename = tkFileDialog.asksaveasfilename(
                initialfile='PSF_%s_%s.fits' %(self.iname, self.filter), 
                filetypes=[('FITS', '.fits')],
                defaultextension='.fits',
                mustexist=True,  # directory chosen must exist.
                parent=self.root)
        if len(filename) > 0:
            self.PSF_HDUlist.writeto(filename) 
            print "Saved to %s" % filename

    def ev_plotspectrum(self):
        "Event handler for Plot Spectrum "
        sptype = self.widgets['SpType'].get()
        iname = self.widgets[self.widgets['tabset'].select()]
        filter = self.widgets[iname+"_filter"].get()

        print "Spectral type is "+sptype
        print "Selected instrument tab is "+iname
        print "Selected instrument filter is "+filter


        P.clf()

        speclib = jwopt.kurucz_stars()

        spectrum = speclib.specFromSpectralType(sptype)
        P.loglog(spectrum['wavelength_um'],spectrum['flux']/spectrum['flux'].max(),label=sptype+" star")
        P.xlabel("Wavelength [$\mu$m]")
        P.ylabel("Flux [$erg s^{-1} cm^{-2} \AA^{-1} \\times$ arbitrary scale factor ]")

        filt = self.instrument[iname].getFilter(filter)
        P.plot(filt.WAVELENGTH, filt.THROUGHPUT+1e-9, "--" ,label=filter+" filter")
        P.gca().set_ybound(1e-6,10)
        P.legend(loc="upper left")

    def ev_calcPSF(self):
        "Event handler for PSF Calculations"
        self._updateFromGUI()
        self.PSF_HDUlist = self.inst.calcPSF(detector_oversample= self.detector_oversampling,
                calc_oversample=self.calc_oversampling,
                fov_arcsec = self.FOV,  nlambda = self.nlambda)
        self.widgets['SaveAs'].state(['!disabled'])

    def ev_displayPSF(self):
        "Event handler for Displaying the PSF"
        self._updateFromGUI()

    def ev_displayOptics(self):
        "Event handler for Displaying the optical system"
        self._updateFromGUI()
        P.clf()
        self.inst.display()

    def ev_displayOPD(self):
        self._updateFromGUI()
        if self.inst.pupilopd is None:
            tkMessageBox.showwarning( message="You currently have selected no OPD file (i.e. perfect telescope) so there's nothing to display.", title="Can't Display") 
        else:
            opd = pyfits.getdata(self.inst.pupilopd[0])

            opd = opd[self.opd_i,:,:] # grab correct slice

            masked_opd = N.ma.masked_equal(opd,  0) # mask out all pixels which are exactly 0, outside the aperture
            cmap = matplotlib.cm.jet
            cmap.set_bad('k', 0.8)
            P.clf()
            P.imshow(masked_opd, cmap=cmap, interpolation='nearest', vmin=-0.5, vmax=0.5)
            P.title("OPD from %s, #%d" %( os.path.basename(self.opd), self.opd_i))
            cb = P.colorbar(orientation='vertical')
            cb.set_label('microns')

            f = P.gcf()
            P.text(0.4, 0.02, "OPD WFE = %6.2f nm RMS" % (masked_opd.std()*1000.), transform=f.transFigure)

    def ev_update_OPD_labels(self):
        "Update the descriptive text for all OPD files"
        for iname in self.instrument.keys():
            self.ev_update_OPD_label(self.widgets[iname+"_opd"], self.widgets[iname+"_opd_label"], iname)

    def ev_update_OPD_label(self, widget_combobox, widget_label, iname):
        "Update the descriptive text displayed about one OPD file"
        #print "Here! for "+iname
        filename = self.instrument[iname]._datapath +os.sep+ 'OPD'+ os.sep+widget_combobox.get()
        if filename.endswith(".fits"):
            #print "read fits %s" % filename
            header_summary = pyfits.getheader(filename)['SUMMARY']
            self.widgets[iname+"_opd_i"]['state'] = 'readonly'
        else:
            #print "not fits: %s" % filename
            header_summary = " 0 nm RMS"
            self.widgets[iname+"_opd_i"]['state'] = 'disabled'

        widget_label.configure(text=header_summary, width=30)

    def _updateFromGUI(self):
        # get GUI values
        self.sptype = self.widgets['SpType'].get()
        self.iname = self.widgets[self.widgets['tabset'].select()]
        self.filter = self.widgets[self.iname+"_filter"].get()
        self.opd= self.widgets[self.iname+"_opd"].get()
        self.opd_i= int(self.widgets[self.iname+"_opd_i"].get())
        self.nlambda= int(self.widgets['nlambda'].get())
        self.FOV= float(self.widgets['FOV'].get())
        self.calc_oversampling= int(self.widgets['calc_oversampling'].get())
        self.detector_oversampling= int(self.widgets['detector_oversampling'].get())


        options = {}
        options['downsample'] = bool(self.vars['downsamp'])
        options['jitter'] = self.widgets['jitter'].get()
        #print "Downsamp value: ",  options['downsample']


        # configure the relevant instrument object
        self.inst = self.instrument[self.iname]
        self.inst.filter = self.filter
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


if __name__ == "__main__":

    gui = JWPSF_GUI()
    gui.mainloop()



