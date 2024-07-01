#!/usr/bin/env python

# This file is part of dsfrb Project.
#
# Copyright (C) 2023 dsfrb Corp.
#
# dsfrb Project is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# dsfrb Project is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with dsfrb Project; if not, write to the Free Software
# Please report any issues or bugs directly to the authors， xiejintao@zhejianglab.com

import sys
import time
import math
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, timezone


class FoldFitsArchive():

    def __init__(self,filename,apply_scl=False,apply_wts=False):
        
        self.filename = filename
        #if mode == 'r':
        self.is_psrfits(self.filename)
        self.hdu = fits.open(self.filename, mode='readonly')
        self.psrfits_version = 0.0
        # folded mode by default
        self.search_mode = False
        # on construction, the data have not been loaded from fits file
        self.loaded_from_fits = True
        #default fraction of the pulse period recorded (in turns)
        self.gate_duty_cycle = 1.0
        #  read header
        self.read_header()
        
        # load parameter from header
        self.load_header()
        # load the subint and data
        self.load_data()
        # load_sub_header
        self.load_sub_header()
        self.get_oreder()
        if apply_scl:
            self.apply_scl()
        
        if apply_wts:
            self.apply_wts()
        
        #elif mode == 'w':
        #    self.hdu = fits.HDUList([fits.PrimaryHDU()])
        
        #else:
        #    raise ValueError("Invalid mode: {}".format(mode))
            
    def read_header(self):
            
        ihdu = 0
        while True:
            try:
                if self.hdu[ihdu].name == 'PRIMARY':
                    self.header = self.hdu[ihdu].header
                    break
                else:
                    ihdu += 1
            except Exception as e:
                print("the file do not have PRIMARY part!")
                sys.exit(1)
        

    def is_psrfits(self,filename):
        """
        Determines whether the data is psrfits data
        """
        try:
            with fits.open(filename) as hdu:
                #  whether the data is psrfits data
                if "PSRFITS" in hdu[0].header["FITSTYPE"]:
                    print("\n")
                    print("Successfully read: {} ".format(filename))
                else:
                    print("Please the input file. This program only supports psrfits file.")
                    sys.exit(1)
        except Exception as e:
            print("Please the input file. This program only supports psrfits file.")
            sys.exit(1)


    def load_header(self):
        '''
        load the header of PRIMARY
        '''
        # file does conform to FITS standard
        self.simple = self.header['SIMPLE']
        # number of bits per data pixel
        self.bitpix = self.header['BITPIX']
        # number of data axes ,default 0
        self.naxis = self.header['NAXIS']
        #FITS dataset may contain extensions
        self.extend = self.header['EXTEND']
        self.comment = self.header['COMMENT']
        # Header version
        self.hdrver = self.header['HDRVER']
        # FITS definition for pulsar data files, PSRFITS
        self.fitstype = self.header['FITSTYPE']
        #  File creation date (YYYY-MM-DDThh:mm:ss UTC)
        self.date = self.header['DATE']
        # Observer name(s)
        self.observer = self.header['OBSERVER']
        # Project name
        self.projid = self.header['PROJID']
        # Telescope name
        self.telescope = self.header['TELESCOP']
        # [m] Antenna ITRF X-coordinate (D)
        self.ant_x = self.header['ANT_X']
        # [m] Antenna ITRF Y-coordinate (D)
        self.ant_y = self.header['ANT_Y']
        # [m] Antenna ITRF Z-coordinate (D)
        self.ant_z = self.header['ANT_Z']
        #  Rx and feed ID
        self.frontend = self.header['FRONTEND']
        # Beam ID for multibeam systems
        try:
            self.ibeamp = self.header['IBEAM']
        except Exception as e:
            self.ibeamp = 0
        # Number of receiver polarisation channels
        self.nrcvr = self.header['NRCVR']
        # LIN or CIRC
        self.fd_poln = self.header['FD_POLN']
        # +/- 1. +1 is LIN:A=X,B=Y, CIRC:A=L,B=R (I)
        self.fd_hand = self.header['FD_HAND']
        #  [deg] FA of E vect for equal sig in A&B (E)
        self.fd_sang = self.header['FD_SANG']
        # [deg] Phase of A^* B for injected cal (E)
        self.fd_xyph = self.header['FD_XYPH']
        #  Backend ID
        self.backend = self.header['BACKEND']
        #  Backend configuration file name
        self.beconfig = self.header['BECONFIG']
        # 0/+1/-1 BE cross-phase:0 unknown,+/-1 std/rev
        self.be_phase = self.header['BE_PHASE']
        # 0/1 BE downconversion conjugation
        self.be_dcc = self.header['BE_DCC']
        # [s] Backend propn delay from digitiser input
        self.be_delay = self.header['BE_DELAY']
        # [s] On-line cycle time (D)
        self.tcycle = self.header['TCYCLE']
        # observe mode: (PSR, CAL, SEARCH)
        self.obs_mode = self.header['OBS_MODE']
    
        if self.obs_mode == "PSR" or self.obs_mode == "LEVPSR":
            self.type = "Pulsar"
        elif self.obs_mode == "CAL" or self.obs_mode == "LEVCAL":
            self.type = "PolnCal"
        elif self.obs_mode == "FOF":
            self.type = "FluxCalOff"
        elif self.obs_mode == "FON":
            self.type = "FluxCalOn"
        elif self.obs_mode == "PCM":
            self.type = "Calibrator"
        elif self.obs_mode == "SEARCH" or self.obs_mode == "SRCH":
            self.search_mode = True
            self.type = "Unknown"
        else:
            self.type = "Unknown"
        
        # Date of observation (YYYY-MM-DDThh:mm:ss U
        self.date_obs = self.header['DATE-OBS']
        # [MHz] Centre frequency for observation
        self.obsfreq = self.header['OBSFREQ']
        # [MHz] Bandwidth for observation
        self.obsbw = self.header['OBSBW']
        # Number of frequency channels (original)
        self.obsnchan = self.header['OBSNCHAN']
        #  DM used to de -disperse each channel (pc/cm^3)
        self.chan_dm = self.header['CHAN_DM']
        #  Name or ID for pointing ctr (multibeam feeds)
        try:
            self.pnt_id = self.header['PNT_ID']
        except Exception as e:
            self.pnt_id = 0
        # Source or scan ID/name
        self.src_name = self.header['SRC_NAME']
        # Coordinate mode (J2000, GAL, ECLIP, etc.)
        self.coor_md = self.header['COORD_MD']
        #  Equinox of coords (e.g. 2000.0)
        self.equinox = self.header['EQUINOX']
        
        if self.coor_md == "EQUAT":
            self.coor_md = "J2000"
        elif  self.coor_md == "Gal" or  self.coor_md == "GALACTIC":
            self.coor_md = "GAL"
        elif self.coor_md == "UNSET":
            self.coor_md = "J2000"
        if self.coor_md == "J2000":
            # If RA exists, set the ra from the value of RA.
            # Otherwise, set the ra to the value of STT_CRD1
            try:
                # Right ascension (hh:mm:ss.ssss)
                self.ra = self.header['RA']
                # Declination (-dd:mm:ss.sss)
                self.dec = self.header['DEC']
            except Exception as e:
                self.ra = self.header['STT_CRD1']
                self.dec = self.header['STT_CRD2']
            # Start coord 1 (hh:mm:ss.sss or ddd.ddd)
            self.stt_crd1 = self.header['STT_CRD1']
            # Start coord 2 (-dd:mm:ss.sss or -dd.ddd)
            self.stt_crd2 = self.header['STT_CRD2']
            # Stop coord 1 (hh:mm:ss.sss or ddd.ddd)
            self.stp_crd1 = self.header['STP_CRD1']
            # Stop coord 2 (-dd:mm:ss.sss or -dd.ddd)
            self.stp_crd2 = self.header['STP_CRD2']
        
        else:
            print("FITSArchive:: not set other coordmode")
       
        # [deg] Beam major axis length
        self.bmaj = self.header['BMAJ']
        # [deg] Beam minor axis length
        self.bmin = self.header['BMIN']
        # [deg] Beam position angle
        self.bpa = self.header['BPA']
        # Track mode (TRACK, SCANGC, SCANLAT)
        self.trk_mode = self.header['TRK_MODE']
        # [s] Requested scan length (E)
        self.scanlen = self.header['SCANLEN']
        #  Feed track mode - FA, CPA, SPA, TPA
        self.fd_mode = self.header['FD_MODE']
        # [deg] Feed/Posn angle requested (E)
        self.fa_req = self.header['FA_REQ']
        # Cal mode (OFF, SYNC, EXT1, EXT2)
        self.cal_mode = self.header['CAL_MODE']
        # [Hz] Cal modulation frequency (E)
        self.cal_freq = self.header['CAL_FREQ']
        #  Cal duty cycle (E)
        self.cal_dcyc = self.header['CAL_DCYC']
        # Cal phase (wrt start time) (E)
        self.cal_phs = self.header['CAL_PHS']
        # Start MJD (UTC days) (J - long integer)
        self.stt_imjd = self.header['STT_IMJD']
        # [s] Start time (sec past UTC 00h) (J)
        self.stt_smjd = self.header['STT_SMJD']
        # [s] Start time offset (D)
        self.stt_offs = self.header['STT_OFFS']
        # [s] Start LST (D)
        self.stt_lst = self.header['STT_LST']

    def load_data(self):
        '''
        load data in search mode
        '''
        if self.search_mode == False:
            #print("\n")
            print("Read the data in Fold mode ...\n")
        else :
            print("\n Warning: this software only supports pulsar data in search mode")
            sys.exit(1)
            
        ihdu = 0
        while True:
            try:
                if self.hdu[ihdu].name == 'SUBINT':
                    self.subint_data = self.hdu[ihdu].data.copy()
                    self.hdu_subint = self.hdu[ihdu]
                    break
                else:
                    ihdu += 1
            except Exception as e:
                print("The {} do not have SUBINT data part, !".format(self.filename))
                sys.exit(1)
        
        self.obsdata = self.subint_data['DATA']
        self.dat_scl = self.subint_data['DAT_SCL']
        self.dat_offs = self.subint_data['DAT_OFFS']
        self.dat_wts = self.subint_data['DAT_WTS']
        self.dat_freq = self.subint_data['DAT_FREQ']
        self.indexval = self.subint_data['INDEXVAL']
        self.tsubint = self.subint_data['TSUBINT']
        self.offs_sub = self.subint_data['OFFS_SUB']
        self.axu_dm = self.subint_data['AUX_DM']
        self.axu_rm = self.subint_data['AUX_RM']
        self.period = self.subint_data['PERIOD']
        
        

 
    def load_sub_header(self):
        '''
        load the subint header
        '''
        self.sub_header = self.hdu_subint.header
        # data and header type: Subintegration data
        self.sub_xtension = self.sub_header['XTENSION']
        # the number of bit of sample
        self.sub_bitpix = self.sub_header['BITPIX']
        # the number of dimensional binary table
        self.sub_naxis = self.sub_header['NAXIS']
        # width of table in bytes
        self.sub_nbyte = self.sub_header['NAXIS1']
        # Number of rows in table (NSUBINT)
        self.sub_nrows = self.sub_header['NAXIS2']
        # size of special data area
        self.sub_pcount = self.sub_header['PCOUNT']
        #  one data group (required keyword)
        self.sub_gcount = self.sub_header['GCOUNT']
        # Number of fields per row
        self.sub_tfields = self.sub_header['TFIELDS']
        # Optionally used if INT_TYPE != TIME
        self.sub_ttype1 = self.sub_header['TTYPE1']
        # the data type of fields 1 :Double
        self.sub_tform1 = self.sub_header['TFORM1']
        #  Length of subintegration
        self.sub_ttype2 = self.sub_header['TTYPE2']
        #the data type of fields 2 :Double
        self.sub_tform2 = self.sub_header['TFORM2']
        # Offset from Start of subint centre
        self.sub_ttype3 = self.sub_header['TTYPE3']
        # the data type of fields 3: Double
        self.sub_tform3 = self.sub_header['TFORM3']
        # label for field
        self.sub_ttype4 = self.sub_header['TTYPE4']
        # the data type of fields 4: Double
        self.sub_tform4 = self.sub_header['TFORM4']
        # additional DM (ionosphere, corona, etc.)
        self.sub_ttype5 = self.sub_header['TTYPE5']
        # the data type of fields 5 : Double
        self.sub_tform5 = self.sub_header['TFORM5']
        # additional RM (ionosphere, corona, etc.)
        self.sub_ttype6 = self.sub_header['TTYPE6']
        # the data type of fields 6 : Double
        self.sub_tform6 = self.sub_header['TFORM6']
        # [MHz] Centre frequency for each channel
        self.sub_ttype7 = self.sub_header['TTYPE7']
        # the data type of fields 7 : Double
        self.sub_tform7 = self.sub_header['TFORM7']
        # Weights for each channel
        self.sub_ttype8 = self.sub_header['TTYPE8']
        # the data type of fields 8 : Float
        self.sub_tform8 = self.sub_header['TFORM8']
        #  Data offset for each channel
        self.sub_ttype9 = self.sub_header['TTYPE9']
        # the data type of fields 9 : Float
        self.sub_tform9 = self.sub_header['TFORM9']
        # Data scale factor (outval=dataval*scl + offs)
        self.sub_ttype10 = self.sub_header['TTYPE10']
        # the data type of fields 10 : Float
        self.sub_tform10 = self.sub_header['TFORM10']
        #  Subint data table
        self.sub_ttype11 = self.sub_header['TTYPE11']
        # the data type of fields 11 : Float
        self.sub_tform11 = self.sub_header['TFORM11']
        #(NBIN,NCHAN,NPOL) or (NCHAN,NPOL,NSBLK*NBITS/8)
        self.sub_tdim11 = self.sub_header['TDIM11']
        #Epoch convention (VALID, MIDTIME, STT_MJD)
        self.sub_epoch = self.sub_header['EPOCHS']
        # Time axis (TIME, BINPHSPERI, BINLNGASC, etc)
        self.sub_int_type = self.sub_header['INT_TYPE']
        # / Unit of time axis (SEC, PHS (0-1), DEG)
        self.sub_int_unit = self.sub_header['INT_UNIT']
        # / Intensity units (FluxDen/RefFlux/Jansky)
        self.sub_scale = self.sub_header['SCALE']
        # Nr of polarisations
        self.sub_npol = self.sub_header['NPOL']
        # Polarisation identifier (e.g., AABBCRCI, AA+BB)
        self.sub_pol_type = self.sub_header['POL_TYPE']
        # [s] Time per bin or sample
        self.sub_tbin = self.sub_header['TBIN']
        # Nr of bins (PSR/CAL mode; else 1)
        self.sub_nbin = self.sub_header['NBIN']
        # Nr of bins/pulse period (for gated data)
        self.sub_nbin_prd = self.sub_header['NBIN_PRD']
        # Phase offset of bin 0 for gated data
        self.sub_phs_offs = self.sub_header['PHS_OFFS']
        # Nr of bits/datum (SEARCH mode 'X' data, else 1)
        self.sub_nbits = self.sub_header['NBITS']
        # Subint offset (Contiguous SEARCH-mode files)
        self.sub_nsuboffs = self.sub_header['NSUBOFFS']
        #  Number of channels/sub-bands in this file
        self.sub_nchan = self.sub_header['NCHAN']
        # [MHz] Channel/sub-band width
        self.sub_chan_bw = self.sub_header['CHAN_BW']
        #  Channel/sub-band offset for split files
        self.sub_nchnoffs = self.sub_header['NCHNOFFS']
        # Samples/row (SEARCH mode, else 1)
        self.sub_nsblk = self.sub_header['NSBLK']
        #Zero offset for SEARCH-mode data
        self.sub_zero_off = self.sub_header['ZERO_OFF']
        # 0 / 1 for signed ints in SEARCH-mode data, else 0
        self.sub_signint = self.sub_header['SIGNINT']
        # [MHz] Reference frequency
        self.sub_reffreq = self.sub_header['REFFREQ']
        # [cm-3 pc] DM for post-detection dedisperion
        self.sub_dm = self.sub_header['DM']
        # [rad m-2] RM for post-detection deFaraday
        self.sub_rm = self.sub_header['RM']
        # Total number of samples (SEARCH mode, else 1)
        self.sub_nstot = self.sub_header['NSTOT']
        # name of this binary table extension
        self.sub_extname = self.sub_header['EXTNAME']
        # Units of field
        try:
            self.sub_tunit1 = self.sub_header['TUNIT1']
        except Exception as e:
            self.sub_tunit1 = 'None'
        # Units of field
        try:
            self.sub_tunit2 = self.sub_header['TUNIT2']
        except Exception as e:
            self.sub_tunit2 = 'None'
        # Units of field
        try:
            self.sub_tunit3 = self.sub_header['TUNIT3']
        except Exception as e:
            self.sub_tunit3 = 'None'
        # Units of field
        try:
            self.sub_tunit4 = self.sub_header['TUNIT4']
        except:
            self.sub_tunit4 = 'None'
        # Units of field
        try:
            self.sub_tunit5 = self.sub_header['TUNIT5']
        except:
            self.sub_tunit5 = 'None'
        # Units of field
        try:
            self.sub_tunit6 = self.sub_header['TUNIT6']
        except Exception as e:
            self.sub_tunit6 = 'None'
        # Units of field
        try:
            self.sub_tunit7 = self.sub_header['TUNIT7']
        except Exception as e:
            self.sub_tunit7 = 'None'
        # Units of field
        try:
            self.sub_tunit8 = self.sub_header['TUNIT8']
        except Exception as e:
            self.sub_tunit8 = 'None'
        # Units of field
        try:
            self.sub_tunit9 = self.sub_header['TUNIT9']
        except Exception as e:
            self.sub_tunit9 = 'None'
        # Units of field
        try:
            self.sub_tunit10 = self.sub_header['TUNIT10']
        except Exception as e:
            self.sub_tunit10 = 'None'
        # Units of field
        try:
            self.sub_tunit11 = self.sub_header['TUNIT11']
        except Exception as e:
            self.sub_tunit11 = 'None'
        # auto assigned by template parser
        self.sub_extver = self.sub_header['EXTVER']
        
    
    def get_oreder(self):
        '''
        gets the data dimension order
        '''
        OrderFPT = str((self.sub_nchan,self.sub_npol,self.sub_nbin))
        OrderTFP = str((self.sub_nbin,self.sub_nchan,self.sub_npol))
        OrderFPT = np.char.replace(OrderFPT, ' ', '')
        OrderTFP = np.char.replace(OrderTFP, ' ', '')
        #print(self.sub_tdim11)
        #print("Here",OrderFPT)
        if self.sub_tdim11 == OrderFPT:
            self.order = "OrderFPT"

        elif self.sub_tdim11 == OrderTFP:
            self.order = "OrderTFP"
        else:
            print("Please check the data dimension order")


    def apply_scl(self):
        '''
        sacle the data to float
        '''
        if self.order == "OrderTFP":
            self.obsdata = self.obsdata.astype(np.float64)
        elif self.order == "OrderFPT":
            self.obsdata = np.transpose(self.obsdata, (2, 0, 1)).astype(np.float64)
            #print("Transpose data...")

        for irow in range(self.sub_nrows):
            for ipol in range(self.sub_npol):
                for ichan in range(self.sub_nchan):
                    sc_of_idx = int(ipol*self.sub_nchan + ichan)
                    self.obsdata[irow][ipol][ichan] = self.obsdata[irow][ipol][ichan] \
                            * self.dat_scl[irow][sc_of_idx]


    
    def apply_wts(self):
        '''
        use the DAT_WTS to zapping the frequency channel
        '''
        
        if self.order == "OrderTFP":
            self.obsdata = self.obsdata
        elif self.order == "OrderFPT":
            self.obsdata = np.transpose(self.obsdata, (2, 0, 1))
       
        for ichan in range(self.sub_nchan):
            if self.dat_wts[:,ichan] ==0:
                self.obsdata[:,:,ichan] =  0#self.obsdata[:,:,ichan] * 0
            else:
                pass
                #print(ichan)
            
    

    def update_keyword(self, hdu_index,keyword, value, comment=None):
        '''
        update the value or comment of keyword
        '''
        self.hdu[hdu_index].header[keyword] = (value, comment)
        
        
        
        
if __name__ == '__main__':
    
    filename = "FRB180301_tracking-M01_0030.ar"
    ## psrfits archive
    arch = FoldFitsArchive(filename,apply_scl=True)
    print(arch.dat_freq)
    print(arch.order)
    print(arch.dat_wts.shape)
    print(arch.obsdata.shape)
    inten = np.sum(arch.obsdata[0,0:2,:,:],axis = 0)
    percentile_lower = 2
    # Upper quantile (e.g., upper 95%)
    percentile_upper = 98.5
    vmin = np.percentile(inten.flatten(), percentile_lower)
    vmax = np.percentile(inten.flatten(), percentile_upper)

    plt.imshow(inten,
        aspect='auto', origin='lower', cmap="viridis", vmin=vmin, vmax=vmax)
    
    plt.show()
    
