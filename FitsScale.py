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

import time
import math
import numpy as np

class Fits_scale():
    def __init__(self,data,nsub,npol,nchan,nbin=1):
        super(Fits_scale, self).__init__()
        
        print("\n Scale the date to int16 \n")
        self.data = data
        self.nsub = nsub
        self.npol = npol
        self.nchan = nchan
        self.nbin = nbin
        self.nsblk = data.shape[-1]
        #self.data = self.dedis.obsdata
        self.the_min = 1-math.pow(2,15)+16
        self.the_max = math.pow(2,15)-2-16
        self.FLT_MAX = 3.4e38
        self.FLT_MIN = 1.175494351e-38
        self.INT16_MIN = -32768
        self.INT16_MAX = 32767
    
        if nbin > 1 :
            print("Downsamp...")
            self.downsamp()
    
        self.measure_scale()
    
    
    def downsamp(self):
        '''
        Do downsampling on the data!
        '''
        self.downsamp_data = np.zeros((self.nsub,self.nchan,
            self.npol,self.nsblk//self.nbin))
        
        
        for isub in range(self.nsub):
            for ipol in range(self.npol):
                for ichan in range(self.nchan):
                
                    self.downsamp_data[isub][ichan][ipol] = np.mean(self.data[isub]
                        [ichan][ipol].reshape(-1, self.nbin),
                        axis=1)
        self.data = self.downsamp_data
        del self.downsamp_data
    
        
    def measure_scale(self):
        '''
        measure the scale and offset
        '''
        # waiting modify
        self.scales = np.zeros((self.nsub,self.npol*self.nchan))
        self.offsets = np.zeros((self.nsub,self.npol*self.nchan))
        for isub in range(self.nsub):
            profile = self.data[isub,:,:,:]
            inde = 0
            for ipol in range(self.npol):
                for ichan in range(self.nchan):
                    #downsamp_amp = np.mean(inten[ichan][ipol].reshape(-1, 2), axis=1)
                    max_val = np.max(profile[ichan][ipol])
                    min_val = np.min(profile[ichan][ipol])
                    # the scale
                    if np.fabs(max_val-min_val) > (100.0 * self.FLT_MIN):
                        iscale = (max_val-min_val) / (self.the_max - self.the_min)
                    else:
                        iscale = 1.0
                    # the offset
                    ioffset = (min_val*self.the_max - max_val*self.the_min)/(self.the_max - self.the_min)
                    
                    self.scales[isub][inde] = iscale
                    self.offsets[isub][inde] = ioffset
                    # scale the data to int16
                    flt_values = (profile[ichan][ipol] - ioffset)/iscale
                    # If any value exceeds the limit of int16, it is set to the limit value
                    flt_values = np.clip(flt_values, self.INT16_MIN, self.INT16_MAX).astype(np.int16)
                    #flt_values = np.int16(flt_values)
                    # Change data
                    self.data[isub,ichan,ipol] = flt_values
                    #print(flt_values.dtype)
                    inde += 1
           
        self.data = self.data.astype(np.int16)


if __name__ == '__main__':
    
    pass
    #Fits_scale()







