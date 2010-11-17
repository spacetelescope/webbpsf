#!/usr/bin/env python
#import os, sys
import numpy as N
import ttk
import Tkinter as tk
import tkMessageBox
from Tkinter import N,E,S,W

import popjw


class JWPSF_GUI(object):
    def __init__(self):

        # init the object and subobjects
        self.instrument = {}
        self.widgets = {}
        self.vars = {}
        insts = ['NIRCam', 'NIRSpec','MIRI', 'TFI', 'FGS']
        self.stars = popjw.kurucz_stars()
        for i in insts:
            self.instrument[i] = popjw.Instrument(i)





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
        self.widgets['SpType'] = ttk.Combobox(lf, textvariable = self.vars['SpType'],width=6)
        self.widgets['SpType'].grid(row=0, column=1)
        self.widgets['SpType']['values'] = self.stars.sptype_list
        self.widgets['SpType'].set('G0V')
        ttk.Label(lf, text='    ' ).grid(row=0, column=2)
        ttk.Button(lf, text='Plot spectrum').grid(row=0,column=3,sticky=E)
        lf.columnconfigure(2, weight=1)
        lf.grid(row=1, sticky=E+W, padx=10,pady=5)

        #-- instruments
        lf = ttk.LabelFrame(frame, text='Instrument Config')
        notebook = ttk.Notebook(lf)
        notebook.pack(fill='both')
        for i in insts:
            page = ttk.Frame(notebook)
            notebook.add(page,text=i)
            ttk.Label(page, text='    Configuration Options for '+i+"                      ").grid(row=0, columnspan=2)

            self.vars[i+"_filter"] = tk.StringVar()
            self.widgets[i+"_filter"] = ttk.Combobox(page,textvariable =self.vars[i+"_filter"], width=10)
            self.widgets[i+"_filter"]['values'] = self.instrument[i].filter_list
            self.widgets[i+"_filter"].set(self.widgets[i+"_filter"]['values'][0])
            #self.widgets[i+"_filter"]['readonly'] = True
            ttk.Label(page, text='    Filter: ' ).grid(row=1, column=0)
            self.widgets[i+"_filter"].grid(row=1, column=1)

            if len(self.instrument[i].image_mask_list) >0 :
                self.vars[i+"_coron"] = tk.StringVar()
                self.widgets[i+"_coron"] = ttk.Combobox(page,textvariable =self.vars[i+"_coron"], width=10)
                masks = self.instrument[i].image_mask_list
                masks.insert(0, "")
                self.widgets[i+"_coron"]['values'] = masks
                ttk.Label(page, text='    Coron: ' ).grid(row=2, column=0)
                self.widgets[i+"_coron"].set(self.widgets[i+"_coron"]['values'][0])
                self.widgets[i+"_coron"].grid(row=2, column=1)

                fr2 = ttk.Frame(page)
                self.vars[i+"_cor_off_r"] = tk.StringVar()
                self.vars[i+"_cor_off_theta"] = tk.StringVar()
                ttk.Label(fr2, text='target offset:  r=' ).grid(row=2, column=4)
                self.widgets[i+"_cor_off_r"] = ttk.Entry(fr2,textvariable =self.vars[i+"_cor_off_r"], width=5)
                self.widgets[i+"_cor_off_r"].insert(0,"0.0")
                self.widgets[i+"_cor_off_r"].grid(row=2, column=5)
                ttk.Label(fr2, text='arcsec,  PA=' ).grid(row=2, column=6)
                self.widgets[i+"_cor_off_theta"] = ttk.Entry(fr2,textvariable =self.vars[i+"_cor_off_theta"], width=3)
                self.widgets[i+"_cor_off_theta"].insert(0,"0")
                self.widgets[i+"_cor_off_theta"].grid(row=2, column=7)
                ttk.Label(fr2, text='deg' ).grid(row=2, column=8)
                fr2.grid(row=2,column=3, sticky=W)

 


            if len(self.instrument[i].image_mask_list) >0 :
                self.vars[i+"_pupil"] = tk.StringVar()
                self.widgets[i+"_pupil"] = ttk.Combobox(page,textvariable =self.vars[i+"_pupil"], width=10)
                masks = self.instrument[i].pupil_mask_list
                masks.insert(0, "")
                self.widgets[i+"_pupil"]['values'] = masks
                ttk.Label(page, text='    Lyot: ' ).grid(row=3, column=0)
                self.widgets[i+"_pupil"].set(self.widgets[i+"_pupil"]['values'][0])
                self.widgets[i+"_pupil"].grid(row=3, column=1)

                fr2 = ttk.Frame(page)
                self.vars[i+"_pupilshift_x"] = tk.StringVar()
                self.vars[i+"_pupilshift_y"] = tk.StringVar()
                ttk.Label(fr2, text=' pupil shift in X:' ).grid(row=3, column=4)
                self.widgets[i+"_pupilshift_x"] = ttk.Entry(fr2,textvariable =self.vars[i+"_pupilshift_x"], width=3)
                self.widgets[i+"_pupilshift_x"].insert(0,"0")
                self.widgets[i+"_pupilshift_x"].grid(row=3, column=5)
                ttk.Label(fr2, text='Y:' ).grid(row=3, column=6)
                self.widgets[i+"_pupilshift_y"] = ttk.Entry(fr2,textvariable =self.vars[i+"_pupilshift_y"], width=3)
                self.widgets[i+"_pupilshift_y"].insert(0,"0")
                self.widgets[i+"_pupilshift_y"].grid(row=3, column=7)
                ttk.Label(fr2, text='% of pupil' ).grid(row=3, column=8)
                fr2.grid(row=3,column=3, sticky=W)








        #notebook.tab('NIRCam').focus_set()
        lf.grid(row=2, sticky=E+W, padx=10, pady=5)

        lf = ttk.LabelFrame(frame, text='OTE Config')
        ttk.Label(lf, text='OTE OPD' ).grid(row=0)
        lf.grid(row=3, sticky=E+W, padx=10, pady=5)


        lf = ttk.LabelFrame(frame, text='Calculation Options')
        r =0
        ttk.Label(lf, text='Field of View:' ).grid(row=r, sticky=W)
        self.widgets['FOV'] = ttk.Entry(lf, width=5)
        self.widgets['FOV'].grid(row=r,column=1, sticky=E)
        self.widgets['FOV'].insert(0,'5')
        ttk.Label(lf, text='arcsec/side' ).grid(row=r, column=2, sticky=W)
        r+=1
        ttk.Label(lf, text='Oversampling:').grid(row=r, sticky=W)
        self.widgets['oversampling'] = ttk.Entry(lf, width=5)
        self.widgets['oversampling'].grid(row=r,column=1, sticky=E)
        self.widgets['oversampling'].insert(0,'2')
        ttk.Label(lf, text='x finer than instrument pixels' ).grid(row=r, column=2, sticky=W, columnspan=2)
        r=r+1
        ttk.Label(lf, text='Output filename:',  ).grid(row=r, sticky=W)
        self.widgets['outfile'] = ttk.Entry(lf, width=40)
        self.widgets['outfile'].grid(row=r,column=1, sticky=E)
        ttk.Label(lf, text='.fits' ).grid(row=r, column=2, sticky=W)
        ttk.Button(lf, text='Save As...', command=self.setSaveAsFilename ).grid(column=3, row=r)


        lf.grid(row=4, sticky=E+W, padx=10, pady=5)

        lf = ttk.Frame(frame)
        ttk.Button(lf, text='Compute PSF', command=self.calcPSF ).grid(column=0, row=0)
        ttk.Button(lf, text='Display PSF', command=self.calcPSF ).grid(column=1, row=0)
        ttk.Button(lf, text='Display Optics', command=self.calcPSF ).grid(column=3, row=0)
        ttk.Button(lf, text='Quit', command=self.quit).grid(column=5, row=0)
        lf.columnconfigure(2, weight=1)
        lf.columnconfigure(4, weight=1)
        lf.grid(row=5, sticky=E+W, padx=10, pady=15)



        #frame.pack(padx=10, pady=10)
        #frame.grid(row=0, sticky=N+E+S+W, padx=10, pady=10)
        frame.grid(row=0, sticky=N+E+S+W)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)



        self.root.update()

    def quit(self):
        if tkMessageBox.askyesno( message='Are you sure you want to quit?', icon='question', title='Install') :
            self.root.destroy()

    def calcPSF(self):
        pass

    def setSaveAsFilename(self):
        pass


    def mainloop(self):
        self.root.mainloop()


if __name__ == "__main__":

    gui = JWPSF_GUI()
    gui.mainloop()



