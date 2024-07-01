#!/usr/bin/env python

# This file is part of PSR/FRB Computing platform Project.
#
# Copyright (C) 2023 PSR/FRB Computing platform Corp.
#
# PSR/FRB Computing platform Project is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# PSR/FRB Computing platform Project is distributed in the hope that it will be useful,
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
from datetime import datetime, timezone
#from numba import jit

class FITSArchive:
    def __init__(self,filename,mode="Fold"):
        self.filename = filename
        #if mode == 'r':
        self.is_psrfits(self.filename)
        self.hdu = fits.open(self.filename, mode='readonly')
        self.psrfits_version = 0.0
        self.chanbw = 0.0
        self.scale_cross_products = False
        self.correct_P236_reference_epoch = False
        # use mode
        self.use_mode = mode
        # search mode by default
        self.search_mode = True
        # no auxiliary profiles
        self.naux_profile = 0
        self.aux_nsample = 0
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

    def read_header(self):
            
        ihdu = 0
        while True:
            try:
                if self.hdu[ihdu].name == 'PRIMARY':
                    self.header = self.hdu[ihdu].header
                    break
                else:
                    ihdu += 1
            except :
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
                    #pass
                    print("\n")
                    #print("Successfully read: {} ".format(filename))
                else:
                    print("Please check the input file. This program only supports psrfits file.")
                    sys.exit(1)
        except:
            print("Please check1 the input file. This program only supports psrfits file.")
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
        except:
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
            self.search_mode = False
        elif self.obs_mode == "CAL" or self.obs_mode == "LEVCAL":
            self.type = "PolnCal"
            self.search_mode = False
        elif self.obs_mode == "FOF":
            self.type = "FluxCalOff"
            self.search_mode = False
        elif self.obs_mode == "FON":
            self.type = "FluxCalOn"
            self.search_mode = False
        elif self.obs_mode == "PCM":
            self.type = "Calibrator"
            self.search_mode = False
        elif self.obs_mode == "SEARCH" or self.obs_mode == "SRCH":
            self.search_mode = True
            self.type = "Unknown"
        else:
            self.search_mode = False
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
        except:
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
            except:
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


    def unpack_1bit(self,data):
        '''
        Unpack 1-bit data that has been read in as bytes.
            Input:
                data: array of bits packed into an array of bytes.
            Output:
                outdata: unpacked array. The size of this array will
                    be eight times the size of the input data.
        '''
        val0 = np.bitwise_and(data >> 0x07, 0x01)
        val1 = np.bitwise_and(data >> 0x06, 0x01)
        val2 = np.bitwise_and(data >> 0x05, 0x01)
        val3 = np.bitwise_and(data >> 0x04, 0x01)
        val4 = np.bitwise_and(data >> 0x03, 0x01)
        val5 = np.bitwise_and(data >> 0x02, 0x01)
        val6 = np.bitwise_and(data >> 0x01, 0x01)
        val7 = np.bitwise_and(data, 0x01)
        
        return np.dstack([val0, val1, val2, val3, val4, val5, val6, val7]).flatten()

    def unpack_2bit(self, data):
        '''
            Unpack 2-bit data that has been read in as bytes.
                Input:
                    data: array of unsigned 2-bit ints packed into
                        an array of bytes.
                Output:
                    outdata: unpacked array. The size of this array will
                        be four times the size of the input data.
        '''
        val0 = np.bitwise_and(data >> 0x06, 0x03)
        val1 = np.bitwise_and(data >> 0x04, 0x03)
        val2 = np.bitwise_and(data >> 0x02, 0x03)
        val3 = np.bitwise_and(data, 0x03)
        return np.dstack([val0, val1, val2, val3]).flatten()

    def unpack_4bit(self, data):
        '''
            Unpack 4-bit data that has been read in as bytes.
            Input:
                data4bit: array of unsigned 4-bit ints packed into
                    an array of bytes.

            Output:
                outdata: unpacked array. The size of this array will
                    be twice the size of the input data.
        '''
        val0 = np.bitwise_and(data >> 0x04, 0x0F)
        val1 = np.bitwise_and(data, 0x0F)
        
        return np.dstack([val0, val1]).flatten()
    

    def load_data(self):
        '''
        load data in search mode
        '''
        if self.search_mode:
            pass
            #print("\n")
            #print("Read the data in search mode ...")
        #if self.search_mode ==False and self.type
        else :
            print("\n Read:  the data in {} mode".format(self.type))
            #sys.exit(1)
        ihdu = 0
        while True:
            try:
                if self.hdu[ihdu].name == 'SUBINT':
                    self.subint_data = self.hdu[ihdu].data
                    self.hdu_subint = self.hdu[ihdu]
                    break
                else:
                    ihdu += 1
            except :
                print("The {} do not have SUBINT data part, !".format(self.filename))
                sys.exit(1)
        
        if self.search_mode:
            bits_per_sample = self.hdu_subint.header['NBITS']
            nsamp_per_subint = self.hdu_subint.header['NSBLK']
        else:
            bits_per_sample = self.hdu_subint.header['BITPIX']
            nsamp_per_subint = self.hdu_subint.header['NBIN']
        
        nchan = self.hdu_subint.header['NCHAN']
        npol = self.hdu_subint.header['NPOL']
        nsubint = self.hdu_subint.header['NAXIS2']
        
        self.obsdata = self.subint_data['DATA']
        shp = self.obsdata.squeeze().shape
        
        
        #print("mode:",self.search_mode)
        #print("sub header: \n",repr(self.hdu_subint.header))
        #print("Data shape:",shp)
        #print("Bit num:",bits_per_sample)
        
        if bits_per_sample <8:
            # Unpack the bytes data
            if (shp[0] != nsamp_per_subint) and \
                    (shp[1] != nchan * npol * bits_per_sample // 8):
                self.obsdata = self.obsdata.reshape((nsubint,nsamp_per_subint,bits_per_sample*
                                        nchan * npol//8))
                                        
            new_data = np.zeros((nsubint,nsamp_per_subint*nchan * npol),dtype=np.uint8)
            print(f"Perform a {bits_per_sample}bit unpack...")
            for isub in range(nsubint):
            #print("step1:", self.obsdata.shape)
                if bits_per_sample == 4:
                    #print("Perform a 4bit unpack...")
                    new_data[isub] = self.unpack_4bit(self.obsdata[isub])
                    #print("unpack 4bit Done. \n")
                elif bits_per_sample == 2:
                    #print("Perform a 2bit unpack...")
                    new_data[isub] = self.unpack_2bit(self.obsdata[isub])
                    #print("unpack 2bit Done. \n")
                elif bits_per_sample == 1:
                    #print("Perform a 1bit unpack...")
                    new_data[isub] = self.unpack_1bit(self.obsdata[isub])
            print(f"unpack {bits_per_sample}bit Done. \n")
            
            self.obsdata = np.asarray(new_data, dtype=np.uint8)
            del new_data
        
        #print("step2:", self.obsdata.shape)
        # reshape the data
        
        if self.search_mode:
            if npol > 1:
                self.obsdata = self.obsdata.reshape((nsubint,nsamp_per_subint, npol, nchan))
            else:
                self.obsdata = self.obsdata.reshape((nsubint,nsamp_per_subint, nchan))
        
        #print("step3:", self.obsdata.shape)
        self.dat_scl = self.subint_data['DAT_SCL']
        self.dat_offs = self.subint_data['DAT_OFFS']
        self.dat_wts = self.subint_data['DAT_WTS']
        self.dat_freq = self.subint_data['DAT_FREQ']
        
        if self.dat_freq.shape[-1] > nchan:
            self.dat_freq = self.dat_freq[:,:nchan]
        else:
            pass
        try:
            self.tel_zen = self.subint_data['TEL_ZEN']
        except:
            pass
        try:
            self.tel_az = self.subint_data['TEL_AZ']
        except:
            pass
        try:
            self.par_ang = self.subint_data['PAR_ANG']
        except:
            pass
        try:
            self.pos_ang = self.subint_data['POS_ANG']
        except:
            pass
        try:
            self.fd_ang = self.subint_data['FD_ANG']
        except:
            pass
        try:
            self.glat_sub = self.subint_data['GLAT_SUB']
        except:
            pass
        try:
            self.glon_sub = self.subint_data['GLON_SUB']
        except:
            pass
        try:
            self.dec_sub = self.subint_data['DEC_SUB']
        except:
            pass
        try:
            self.ra_sub = self.subint_data['RA_SUB']
        except:
            pass
        try:
            self.lst_sub = self.subint_data['LST_SUB']
        except:
            pass
        self.offs_sub = self.subint_data['OFFS_SUB']
        self.tsubint = self.subint_data['TSUBINT']
            
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
        # Length of subintegration
        self.sub_ttype1 = self.sub_header['TTYPE1']
        # the data type of fields 1 :Double
        self.sub_tform1 = self.sub_header['TFORM1']
        # Offset from Start of subint centre
        self.sub_ttype2 = self.sub_header['TTYPE2']
        #the data type of fields 2 :Double
        self.sub_tform2 = self.sub_header['TFORM2']
        # tLST at subint centre
        self.sub_ttype3 = self.sub_header['TTYPE3']
        # the data type of fields 3: Double
        self.sub_tform3 = self.sub_header['TFORM3']
        # RA (J2000) at subint centre
        self.sub_ttype4 = self.sub_header['TTYPE4']
        # the data type of fields 4: Double
        self.sub_tform4 = self.sub_header['TFORM4']
        # Dec (J2000) at subint centre
        self.sub_ttype5 = self.sub_header['TTYPE5']
        # the data type of fields 5 : Double
        self.sub_tform5 = self.sub_header['TFORM5']
        # [deg] Gal longitude at subint centre
        self.sub_ttype6 = self.sub_header['TTYPE6']
        # the data type of fields 6 : Double
        self.sub_tform6 = self.sub_header['TFORM6']
        # [deg] Gal latitude at subint centre
        self.sub_ttype7 = self.sub_header['TTYPE7']
        # the data type of fields 7 : Double
        self.sub_tform7 = self.sub_header['TFORM7']
        # [deg] Feed angle at subint centre
        self.sub_ttype8 = self.sub_header['TTYPE8']
        # the data type of fields 8 : Float
        self.sub_tform8 = self.sub_header['TFORM8']
        #  [deg] Position angle of feed at subint centre
        self.sub_ttype9 = self.sub_header['TTYPE9']
        # the data type of fields 9 : Float
        self.sub_tform9 = self.sub_header['TFORM9']
        # [deg] Parallactic angle at subint centre
        self.sub_ttype10 = self.sub_header['TTYPE10']
        # the data type of fields 10 : Float
        self.sub_tform10 = self.sub_header['TFORM10']
        
        try:
            # [deg] Telescope azimuth at subint centre
            self.sub_ttype11 = self.sub_header['TTYPE11']
            # the data type of fields 11 : Float
            self.sub_tform11 = self.sub_header['TFORM11']
        except:
            pass
        try:
            #  [deg] Telescope zenith angle at subint
            self.sub_ttype12 = self.sub_header['TTYPE12']
            #  the data type of fields 12 : Float
            self.sub_tform12 = self.sub_header['TFORM12']
        except:
            pass
        try:
            # [MHz] Centre frequency for each channel
            self.sub_ttype13 = self.sub_header['TTYPE13']
            #  the data type of fields 13 : Nchan Float
            self.sub_tform13 = self.sub_header['TFORM13']
        except:
            pass
        try:
            # Weights for each channel
            self.sub_ttype14 = self.sub_header['TTYPE14']
            # the data type of fields 14 : Nchan Float
            self.sub_tform14 = self.sub_header['TFORM14']
        except:
            pass
        try:
            # Data offset for each channel
            self.sub_ttype15 = self.sub_header['TTYPE15']
            # the data type of fields 15: NCHAN*NPOL Floats
            self.sub_tform15 = self.sub_header['TFORM15']
        except:
            pass
        try:
            # Data scale factor for each channel
            self.sub_ttype16 = self.sub_header['TTYPE16']
            # the data type of fields 16: NCHAN*NPOL Floats
            self.sub_tform16 = self.sub_header['TFORM16']
        except:
            pass
        try:
            # Subint data table
            self.sub_ttype17 = self.sub_header['TTYPE17']
            # the data type of fields 17: NBIN*NCHAN*NPOL*NSBLK int, byte(B) or bit(X)
            self.sub_tform17 = self.sub_header['TFORM17']
        except:
            pass
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
        # name of this binary table extension
        self.sub_extname = self.sub_header['EXTNAME']
        # Units of field
        try:
            self.sub_tunit1 = self.sub_header['TUNIT1']
        except:
            self.sub_tunit1 = 'None'
        # Units of field
        try:
            self.sub_tunit2 = self.sub_header['TUNIT2']
        except:
            self.sub_tunit2 = 'None'
        # Units of field
        try:
            self.sub_tunit3 = self.sub_header['TUNIT3']
        except:
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
        except:
            self.sub_tunit6 = 'None'
        # Units of field
        try:
            self.sub_tunit7 = self.sub_header['TUNIT7']
        except:
            self.sub_tunit7 = 'None'
        # Units of field
        try:
            self.sub_tunit8 = self.sub_header['TUNIT8']
        except:
            self.sub_tunit8 = 'None'
        # Units of field
        try:
            self.sub_tunit9 = self.sub_header['TUNIT9']
        except:
            self.sub_tunit9 = 'None'
        # Units of field
        try:
            self.sub_tunit10 = self.sub_header['TUNIT10']
        except:
            self.sub_tunit10 = 'None'
        # Units of field
        try:
            self.sub_tunit11 = self.sub_header['TUNIT11']
        except:
            self.sub_tunit11 = 'None'
        # Units of field
        try:
            self.sub_tunit12 = self.sub_header['TUNIT12']
        except:
            self.sub_tunit12 = 'None'
        # Units of field
        try:
            self.sub_tunit13 = self.sub_header['TUNIT13']
        except:
            self.sub_tunit13 = 'None'
        # Dimensions (NBITS or NBIN,NCHAN,NPOL,NSBLK)
        try:
            self.sub_tdim17 = self.sub_header['TDIM17']
        except:
            self.sub_tdim17 = 'None'
        # Units of subint data
        try:
            self.sub_tunit17 = self.sub_header['TUNIT17']
        except:
            self.sub_tunit17 = 'None'
        # auto assigned by template parser
        try:
            self.sub_extver = self.sub_header['EXTVER']
        except:
            self.sub_extver = 'None'
    
    def get_oreder(self):
        '''
        gets the data dimension order
        '''
        
        if self.search_mode:
            nsamp_per_subint = self.hdu_subint.header['NSBLK']
        else:
            nsamp_per_subint = self.hdu_subint.header['NBIN']
        
        if self.sub_npol>1:
            
            if self.search_mode:
                OrderFPT = str((1,self.sub_nchan,self.sub_npol,nsamp_per_subint))
                OrderTFP = str((1,nsamp_per_subint,self.sub_nchan,self.sub_npol))
            else:
                OrderFPT = str((1,self.sub_nchan,self.sub_npol,nsamp_per_subint))
                OrderTFP = str((nsamp_per_subint,self.sub_nchan,self.sub_npol,self.sub_nsblk))
                
            OrderFPT = np.char.replace(OrderFPT, ' ', '')
            OrderTFP = np.char.replace(OrderTFP, ' ', '')
            OrderFPT2 = OrderFPT
            OrderTFP2 = OrderTFP
            
        else:
            if self.search_mode:
                OrderFPT = str((1,self.sub_nchan,nsamp_per_subint))
                OrderTFP = str((1,nsamp_per_subint,self.sub_nchan,self.sub_nsblk))
            else:
                OrderFPT = str((1,self.sub_nchan,nsamp_per_subint))
                OrderTFP = str((nsamp_per_subint,self.sub_nchan))
                
            OrderFPT = np.char.replace(OrderFPT, ' ', '')
            OrderTFP = np.char.replace(OrderTFP, ' ', '')
            
            OrderFPT2 = str((1,self.sub_nchan,self.sub_npol,nsamp_per_subint))
            OrderTFP2 = str((1,nsamp_per_subint,self.sub_nchan,self.sub_npol))
            OrderFPT2 = np.char.replace(OrderFPT2, ' ', '')
            OrderTFP2 = np.char.replace(OrderTFP2, ' ', '')
            
        if self.sub_tdim17 == 'None':
            if self.obsdata.shape != 1:
                datashape = list(self.obsdata.shape)
                datashape.append(1)
                datashape = tuple(datashape)
            self.sub_tdim17 = str(datashape[::-1][:-1])
            self.sub_tdim17 = np.char.replace(self.sub_tdim17, ' ', '')
        #print("Here",OrderFPT)
        if self.sub_tdim17 == OrderFPT or self.sub_tdim17 == OrderFPT2:
            self.order = "OrderFPT"
        elif self.sub_tdim17 == OrderTFP or self.sub_tdim17 == OrderTFP2:
            self.order = "OrderTFP"
        else:
            print("Please check the data dimension order",self.sub_tdim17,self.obsdata.shape,OrderTFP,OrderFPT)
            sys.exit(1)
            
            
        
    def update_keyword(self, hdu_index,keyword, value, comment=None):
        '''
        update the value or comment of keyword
        '''
        self.hdu[hdu_index].header[keyword] = (value, comment)
        
        
        
        
if __name__ == '__main__':
    
    filename = "FRB180301_tracking-M01_0030.fits"
    ## psrfits archive
    arch = FITSArchive(filename)
    print(arch.dat_freq)
    print(arch.order)
