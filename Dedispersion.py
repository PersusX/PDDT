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
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from PsrFitsArchive import FITSArchive
from FitsScale import Fits_scale as fisc


class channel_zap_median():
    def __init__(self,spectrum,rms_threshold = 4.0,
            window_size = 21):
        super(channel_zap_median,self).__init__()
        # Default constructor
        self.rms_threshold = rms_threshold
        self.madm_threshold = 0.0
        self.iqr_threshold = 0.0
        self.window_size = window_size
        self.bybin = False
        self.from_total = False
        self.paz_report = False
        self.single_mask = None
        self.statistic = None
        self.expression = ""
        self.spectrum = spectrum
        self.nchan = len(spectrum)
        self.mask = [False]*self.nchan
        #self.median_smooth(spectrum,self.window_size)
        self.rms_zap()

        
    def median_smooth(self,data,wsize):
        '''
        Median filter function
        '''
        # 检查输入数据的长度是否满足要求
        if len(data) < 4:
            print("The number of data points is too small, less than 4")
            sys.exit(1)
        # 检查滑动窗口大小是否满足要求
        if wsize < 3:
            print("Sliding window is too small, please set to > 3")
            sys.exit(1)
        # 检查输入数据长度是否大于窗口大小加一
        if len(data) < wsize + 1:
            print("Error: The number of data points is too small, less than wsize + 1")
            sys.exit()
        # 将窗口大小调整为奇数
        if wsize % 2 == 0:
            wsize += 1
        # 初始化变量
        middle = wsize // 2
        window = np.zeros(wsize)
        result = np.zeros(len(data))
        truncated = 0
        tmid = 0
        # 处理数据的前部分
        for ipt in range(middle):
            truncated = middle + 1 + ipt
            tmid = truncated // 2
            # 将窗口内的数据复制到临时数组 window
            window[:truncated] = data[:truncated]
            #找到窗口内的中值
            window[:truncated].sort()
            median = window[tmid]
            # 将中值赋值给结果数组
            result[ipt] = median
        # 处理数据的中间部分
        rsize = len(data) - wsize + 1
        for ipt in range(rsize):
            # 将窗口内的数据复制到临时数组 window
            window[:] = data[ipt : ipt + wsize]
            # 找到窗口内的中值
            window.sort()
            median = window[middle]
            # 将中值赋值给结果数组
            result[ipt + middle] = median
        # 处理数据的后部分
        for ipt in range(rsize, len(data) - middle):
            truncated = len(data) - ipt
            tmid = truncated // 2
            # 将窗口内的数据复制到临时数组 window
            window[:truncated] = data[ipt : ipt + truncated]
            # 找到窗口内的中值
            window[:truncated].sort()
            median = window[tmid]
            # 将中值赋值给结果数组
            result[ipt + middle] = median
        return result
 
 
    def rms_zap(self):
        # 执行中值平滑
        smoothed_spectrum = self.median_smooth(self.spectrum, self.window_size)
        spectrum = self.spectrum.copy()
        variance = 0.0
        total_chan = 0
        orig_total_chan = 0
        for ichan in range(self.nchan):
            spectrum[ichan] -= smoothed_spectrum[ichan]
            spectrum[ichan] *= spectrum[ichan]
            if not self.mask[ichan]:
                variance += spectrum[ichan]
                total_chan += 1
        variance /= (total_chan - 1)
        orig_total_chan = total_chan
        zapped = True
        round = 1
        while zapped:
            cutoff = self.rms_threshold * self.rms_threshold * variance
            zapped = False
            round += 1
            for ichan in range(self.nchan):
                if self.mask[ichan]:
                    continue
                if abs(spectrum[ichan]) > cutoff:
                    self.mask[ichan] = True
                if abs(spectrum[ichan] - spectrum[ichan-1]) > 2 * cutoff:
                    self.mask[ichan] = True
                if self.mask[ichan]:
                    variance -= spectrum[ichan] / (orig_total_chan - 1.)
                    total_chan -= 1
                    zapped = True
            variance *= (orig_total_chan - 1.) / (total_chan - 1.)
            orig_total_chan = total_chan









class Dedispersion():
    
    def __init__(self, filename, dm,threshold_snr=7, swt=0.3, savefig=False):
        super(Dedispersion, self).__init__()
        self.filename = filename
        self.basename = filename.split('.')[0]
        self.dm = dm
        self.thd_snr = threshold_snr
        self.searchWindowSize = swt # unit: s
        
        
        
        
        self.savefig = savefig
        
        self.arch = FITSArchive(self.filename)
        self.delay_t()
        self.dedispersion()        
        #self.hist_hdu = create_hdu_history().hdu
        #print(repr(self.hist_hdu.header))
        #self.arch.hdu.close()

    def delay_t(self):
        '''
        calculate the delay samples in each channel
        DM: Dispersion measure, in pc/cm^3
        freqs: observed freqency in each channel in MHz
        tbin: Sampling time in s
        '''
        # Reference frequency, using the highest frequency
        #freq_ref = freqs[len(freqs)//2]
        freq_ref = max(self.arch.dat_freq[0])
        
        print("the reference frequency:",freq_ref,self.arch.sub_tbin)
        #According to the given dm value, The time correction for each frequency channel is calculated
        self.t_delays = (4.148808 * self.dm * (self.arch.dat_freq[0] ** -2 - freq_ref ** -2) * 1e3 / self.arch.sub_tbin).astype(np.int64)
        # the delay time between the highest channel center frequency and the lowest center frequency
        self.max_delay = np.max(self.t_delays)
        
        
    def dedispersion(self):
        # Number of polarization
        npol = self.arch.sub_npol
        # Number of channel
        nchan = self.arch.sub_nchan
        # Number of sample in each row
        nsamp = self.arch.sub_nsblk
        # Number of rows
        nrows = self.arch.sub_nrows
        # the total number of sample in each channnel
        nsamp_ch = nrows * nsamp
        
        if self.arch.order == "OrderFPT":
            self.arch.obsdata =  self.arch.obsdata.reshape(
                nrows * nsamp, npol, nchan).T
        else:
            print("The program does not support data of other dimensional order at this time. Only TFP!")
            
        print("Start dedispersion ")
        new_data = np.zeros((nchan, npol, nsamp_ch - self.max_delay))
        sys.stdout.flush()
        oldpcnt = ""
        start_time = time.time()
        if nchan >=1:
            for ichan in range(nchan):
                # Refer to the delay time of the highest frequency
                ifreq_delay = self.t_delays[ichan]
                if ifreq_delay >= 0:
                    # Perform the dispersion process
                    new_data[ichan, :, :] = self.arch.obsdata[ichan, :,
                        ifreq_delay: ifreq_delay + nsamp_ch - self.max_delay]
                else:
                    print("Perform the dedispersion process error...")
                    #data[ichan,:,-ifreq_delay:] = data[ichan,:,:nsamp_ch+ifreq_delay]
                pcnt = "%d" %(100*(ichan+1)/nchan)
                #(100*(1/npol*(ichan+1)/nchan+ipol/npol))#(100/npol*ipol+ichan/nchan*npol)
                if pcnt != oldpcnt:
                    sys.stdout.write("% 4s%% complete\r" % pcnt)
                    sys.stdout.flush()
        else:
            print("There is channel wrong with {} file".format(self.filename))
        
        
        del self.arch.obsdata
        # set the fold mode shape
        self.arch.obsdata = new_data.reshape(1,nchan,npol,nsamp_ch-self.max_delay)
        #self.arch.obsdata = np.transpose(self.arch.obsdata, (0,2,1,3))
        
        del new_data
        
        end_time = time.time()
        spend_time = end_time- start_time
        
        print("Dedispersion done!")
        print(f"the dedispersion process spend time: {spend_time}")
    


        # search pulse 

        if self.arch.obsdata.shape[-1] > 4096:
            print("\n Search pulse position... \n")
            profile = np.sum(np.sum(self.arch.obsdata[0,:,0:2,:],axis=1),axis=0)
            
            self.arch.sub_nsblk = 4096
            pro_index = np.argmax(profile)
            self.arch.obsdata = self.arch.obsdata[:,:,:,pro_index-self.arch.sub_nsblk//2:pro_index+self.arch.sub_nsblk//2]
            #print(inten.shape)
            #print(self.arch.sub_tbin)
            self.arch.header["STT_SMJD"] += int((pro_index-self.arch.sub_nsblk//2)*self.arch.sub_tbin)
            stt_offs = ((pro_index-self.arch.sub_nsblk//2)*self.arch.sub_tbin
                    - int((pro_index-self.arch.sub_nsblk//2)*self.arch.sub_tbin))
            if stt_offs + self.arch.header["STT_OFFS"] >1:
                self.arch.header["STT_SMJD"] += 1
                self.arch.header["STT_OFFS"] += stt_offs + self.arch.header["STT_OFFS"] -1
            else :
                self.arch.header["STT_OFFS"] = stt_offs
                
            
        else:
            self.arch.sub_nsblk = nsamp_ch - self.max_delay
        # set the number of each channel
        
        self.sub_nrows = 1
        
        if self.savefig == True:
            self.save_figure()
        #scale the date to int16
        fits_scal =  fisc(self.arch.obsdata,self.sub_nrows,npol,nchan)
        self.arch.obsdata = fits_scal.data
        
        self.arch.obsdata = np.transpose(self.arch.obsdata, (0,2,1,3))
        # set DAT_SCL
        self.arch.dat_scl = fits_scal.scales
        # set DAT_OFFS
        self.arch.dat_offs = fits_scal.offsets
        self.arch.dat_wts = np.ones((self.sub_nrows,nchan),dtype=np.float32)
        self.arch.sub_header["DM"] = self.dm
        
        now_time  = datetime.now(timezone.utc)
        str_now_time = now_time.strftime('%Y-%m-%dT%H:%M:%S')
        self.arch.header["DATE"] = str_now_time
        
        
    def update_hdu(self):
        
        self.arch.header['COMMENT'] = 'FITS dataset may contain extensions', 'FITS (Flexible Image Transport System) format is defined in \'Astronomy and Astrophysics\', volume 376, page 359; bibcode: 2001A&A...376..359H FITS (Flexible Image Transport System) format defined in \'Astronomy and Astrophysics\' Supplement Series v44/p363, v44/p371, v73/p359, v73/p365. Contact the NASA Science Office of Standards and Technology for the FITS Definition document #100 and other FITS information.'

        # update the create file time
        now_time  = datetime.now(timezone.utc)
        str_now_time = now_time.strftime('%Y-%m-%dT%H:%M:%S')
        self.arch.header['DATE'] = str_now_time
        # update the observe mode
        self.arch.header['OBS_MODE'] = 'PSR'
        # update the header vision
        self.arch.header['HDRVER'] = '6.6'
        self.reset_columns()
        
        
    def multiplyList(myList):
      
        res = 1
        for i in myList:
            res = res * i
        return res

    def reset_columns(self):
        '''
        reset the colnmns
        '''
        
        columns = []
        # the number of subint
        nsub = self.obsdata.shape[0]
        #period = nsub * self.arch.sub_tbin * self.arch.sub_nsblk
        # the index of each sub
        indexval = np.zeros(nsub)
        columns.append(fits.Column(name = 'INDEXVAL',
            format = str(nsub)+"D", array =indexval))
        # the time length of subint
        tsubint = np.array([self.arch.sub_tbin*self.obsdata.shape[2]]*nsub)
        columns.append(fits.Column(name = 'TSUBINT',
            format = str(nsub)+"D", unit = 's',array = tsubint))
        # the time offset of subint
        offs_sub = np.zeros(nsub)
        columns.append(fits.Column(name = 'OFFS_SUB',
            format = str(nsub)+"D", unit = 's',array = offs_sub))
        # the period , defdault: the time length of subint
        period = np.array([self.arch.sub_tbin*self.obsdata.shape[2]]*nsub)
        columns.append(fits.Column(name = 'PERIOD',
            format = str(nsub)+"D", unit = 's',array = period))
        # Auxiliary radio observation parameter DM for frb or pulsar
        aux_dm = np.zeros(nsub)
        columns.append(fits.Column(name = 'AUX_DM',
            format = str(nsub)+"D", unit = 'CM-3',array = aux_dm))
        # Auxiliary radio observation parameter RM for frb or pulsar
        aux_rm = np.zeros(nsub)
        columns.append(fits.Column(name = 'AUX_RM',
            format = str(nsub)+"D", unit = 'RAD',array = aux_rm))
        # Observed frequency
        dat_freq = np.tile(self.arch.dat_freq[0], (nsub, 1))
        columns.append(fits.Column(name = 'DAT_FREQ',
            format = str(self.arch.sub_ncha)+"D", unit = 'MHz',array = dat_freq))
        # Weight value of each channel, the weights of each polarization is the same
        dat_wts = np.ones((nsub,self.arch.sub_nchan))
        columns.append(fits.Column(name = 'DAT_WTS',
            format = str(self.arch.sub_ncha)+"E", array = dat_wts))
        # Offset value of each channel
        dat_offs = np.zeros((nsub,self.arch.sub_npol,self.arch.sub_nchan))
        columns.append(fits.Column(name = 'DAT_OFFS',
            format = str(self.arch.sub_ncha*self.arch.sub_npol)+"E", array = dat_offs))
        # Scaling factor for the observed data
        dat_scl = np.ones((nsub, self.arch.sub_npol, self.arch.sub_nchan))
        columns.append(fits.Column(name = 'DAT_SCL',
            format = str(self.arch.sub_ncha*self.arch.sub_npol)+"E", array = dat_scl))
        # observed data, Data after Dedispersion
        
        dataDim = str(self.arch.sub_ncha*self.arch.sub_npol*self.arch.sub_nsblk)+'I'
        columns.append(fits.Column(name = 'DATA',
            format = dataDim, unit = 'Jy', array = self.obsdata))
        
        sellf.subint_columns = fits.ColDefs(columns)
        
        
    def save_figure(self):
        '''
        Output the image after the dispersion
        '''
        #if self.arch.sub_nrows * self.arch.sub_nsblk > 2048:
            
       #     profile = np.sum(np.sum(self.arch.obsdata[0,:,0:2,:],axis=1),axis=0)
       #     pro_index = np.argmax(profile)
            
       #     inten = np.sum(self.arch.obsdata[0,:,0:2,pro_index-1024:pro_index+1024],axis=1)
            #print(inten.shape)
            #print(self.arch.sub_tbin)
        #else:
        inten = np.sum(self.arch.obsdata[0,:,0:2,:],axis=1)
        # set x-axis and y-axis coordinate
        x_time, y_freq = np.meshgrid(np.arange(inten.shape[1]) * self.arch.sub_tbin *
            1e3, self.arch.dat_freq[0])
        
        
        vmin, vmax = np.sort(inten.flatten())[int(inten.shape[0]*
            inten.shape[1]/50)], np.sort(inten.flatten())[int
                (inten.shape[0]*inten.shape[1]/50*49)]
        
        print(inten.shape,x_time.shape,y_freq.shape)
        plt.pcolormesh(x_time, y_freq, inten, cmap='viridis', vmin=vmin, vmax=vmax)
        plt.xlabel('Time (ms)')
        plt.ylabel('Frequence (MHz)', labelpad=6)
        #plt.imshow(inten
        #    aspect='auto', origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        plt.savefig(self.basename+"dedis"+".png",dpi=300,format="png")
        plt.close()
        
    
if __name__ == '__main__':
        
    
    filename = "FRB180301_tracking-M01_0030.fits"
    Dedispersion("FRB180301_tracking-M01_0030.fits", 518, savefig = False)
    
    
    
    
    
    















