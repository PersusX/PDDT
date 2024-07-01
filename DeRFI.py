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
import os
import sys
import numpy as np

class channel_zap_median:
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
        # check if the input is of sufficient length
        if len(data) < 4:
            print("The number of data points is too small, less than 4")
            sys.exit(1)
        # check if the sliding window size is sufficient
        if wsize < 3:
            print("Sliding window is too small, please set to > 3")
            sys.exit(1)
        # check if input data length is greater than window size plus one
        if len(data) < wsize + 1:
            print("Error: The number of data points is too small, less than wsize + 1")
            sys.exit()
        # Resize the window to odd
        if wsize % 2 == 0:
            wsize += 1
        # Initialize variables
        middle = wsize // 2
        window = np.zeros(wsize)
        result = np.zeros(len(data))
        truncated = 0
        tmid = 0
        # Process the first part of the data
        for ipt in range(middle):
            truncated = middle + 1 + ipt
            tmid = truncated // 2
            # Copy the data inside the window to the temporary array window
            window[:truncated] = data[:truncated]
            # find the median value within the window
            window[:truncated].sort()
            median = window[tmid]
            # assign the median to the resulting array
            result[ipt] = median
        # process the middle part of the data
        rsize = len(data) - wsize + 1
        for ipt in range(rsize):
            # 将窗口内的数据复制到临时数组 window
            window[:] = data[ipt : ipt + wsize]
            # find the median value within the window
            window.sort()
            median = window[middle]
            # assign the median to the resulting array
            result[ipt + middle] = median
        # process the next part of the data
        for ipt in range(rsize, len(data) - middle):
            truncated = len(data) - ipt
            tmid = truncated // 2
            # Copy the data inside the window to the temporary array window
            window[:truncated] = data[ipt : ipt + truncated]
            # find the median value within the window
            window[:truncated].sort()
            median = window[tmid]
            # assign the median to the resulting array
            result[ipt + middle] = median
        return result
 
    def rms_zap(self):
        # perform median smoothing
        smoothed_spectrum = self.median_smooth(self.spectrum, self.window_size)
        spectrum = self.spectrum.copy()
        variance = 0.0
        total_chan = 0
        orig_total_chan = 0
        for ichan in range(self.nchan):
            spectrum[ichan] -= smoothed_spectrum[ichan]
            #print(spectrum[ichan],smoothed_spectrum[ichan])
            
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




class channel_zap_fft:
    def __init__(self,data,factor=0.05):
        # data order: [chan, sample]
        self.data = data
        # The scaling factor of the frequency band is deleted
        self.factor = factor
        self.nchan = self.data.shape[0]
        self.mask = [False]*self.nchan
        self.fft_zap()
        
    def get_indexes_of_max_percentage(self,data, percentage):
        # calculate the length of the section of the list to get
        part_length = int(len(data) * percentage)
        # sort the data and extract the largest part
        sorted_data = sorted(data, reverse=True)  # 降序排序以获取最大部分
        max_part = sorted_data[:part_length]
        # Get the index of the largest value in the original list
        indexes = [data.index(value) for value in max_part]
        return indexes

    def find_strongest_frequency(self,channel_data):
        # do FFT
        fft_result = np.fft.fft(channel_data)
        frequencies = np.fft.fftfreq(len(channel_data))
        # Find the index of the maximum intensity in the FFT result except 0 Hz
        fft_magnitudes = np.abs(fft_result)
        # Set the magnitude at the 0 Hz position to 0 and exclude the 0 Hz component
        fft_magnitudes[0] = 0
        max_intensity_index = np.argmax(fft_magnitudes)
        # Get the frequency corresponding to the maximum intensity
        max_intensity_frequency = frequencies[max_intensity_index]
        return [max_intensity_frequency,max(fft_magnitudes)]

    def fft_zap(self):
        Peri_sig_stre = []
        strongest_frequencies = []
        for channel_data in self.data:
            strongest_frequency_stre = self.find_strongest_frequency(channel_data)
            strongest_frequencies.append(strongest_frequency_stre[0])
            Peri_sig_stre.append(strongest_frequency_stre[1])
        max_idx = self.get_indexes_of_max_percentage(Peri_sig_stre, self.factor)
        for ichan in max_idx:
            self.mask[ichan] = True
        del self.data
