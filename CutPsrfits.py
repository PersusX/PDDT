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

import os
import re
import sys
import numpy as np
import multiprocessing
from astropy.io import fits
#sys.path.append(os.path.dirname(__file__))

#from PsrFitsArchive import FITSArchive


class CutPsrfits:
    def __init__(self, psrfitsfile,stt_subint=0,end_subint=10,progress_queue='',
                    outfilename=None):
   
        filename = psrfitsfile
        self.hdul = fits.open(filename)
        self.start_subint = stt_subint
        self.end_subint = end_subint
        
        if outfilename == None:
            out_basename = filename.split('.')[-2]
            directory_path = os.path.dirname(psrfitsfile)
            self.outfilename = out_basename + '_cut_sub_{}_{}.fits'.format(self.start_subint,self.end_subint)
        else:
            self.outfilename = outfilename
        
        self.cut_subints()

    def cut_subints(self):
        
        # 获取 SUBINT 表
        subint_hdu = self.hdul['SUBINT']
        subint_data = subint_hdu.data
        # 切割数据
        cut_data = subint_data[self.start_subint:self.end_subint]
       
       
        stt_imjd = self.hdul['PRIMARY'].header["STT_IMJD"]
        stt_smjd = self.hdul['PRIMARY'].header["STT_SMJD"]
        stt_offs = self.hdul['PRIMARY'].header["STT_OFFS"]
        
        nsamp_subint = self.hdul['SUBINT'].header["TBIN"]
        tsamp = self.hdul['SUBINT'].header["NSBLK"]
        # 计算开始时间的秒数偏移
        total_seconds = stt_smjd + stt_offs + self.start_subint* nsamp_subint * tsamp
        new_imjd = stt_imjd
        new_smjd = int(total_seconds)
        new_offs = total_seconds - new_smjd


        # 创建新的 HDU 列表
        new_hdul = fits.HDUList([self.hdul[0]])
        
        # 复制所有 HDUs，但替换 SUBINT
        for hdu in self.hdul[1:]:
            if hdu.header['EXTNAME'] == 'SUBINT':
                new_subint_hdu = fits.BinTableHDU(data=cut_data, header=hdu.header)
                new_hdul.append(new_subint_hdu)
            else:
                new_hdul.append(hdu)
        
        # 更新相关表头参数
        #new_hdul[0].header['NAXIS2'] = self.end_subint - self.start_subint
        new_hdul['SUBINT'].header['NAXIS2'] = self.end_subint - self.start_subint

        # 更新 MJD 参数
        new_hdul[0].header['STT_IMJD'] = new_imjd
        new_hdul[0].header['STT_SMJD'] = new_smjd
        new_hdul[0].header['STT_OFFS'] = new_offs

        # 保存新的 FITS 文件
        new_hdul.writeto(self.outfilename, overwrite=True)
    
        print("Save Cut data to {}.".format(self.outfilename))
    
    
    def close(self):
        self.hdul.close()




def Usage():

    print("\n")
    print("-"*10)
    print("\n")
    print("Usage: python CutPsrfits.py xxx.fits -sttsub 0 -endsub 10  \n")
    print("-sttsub :  get the start subint index value from input argv \n")
    print("-endsub :  get the end subint index value from input argv\n")
    print("-h : help information \n")
    print("-"*10)
    print("\n")


def is_psrfits(filename):
    """
    Determines whether the data is psrfits data
    """
    try:
        with fits.open(filename) as hdu:
            
            #  whether the data is psrfits data
            if "PSRFITS" in hdu[0].header["FITSTYPE"]:
                
                return True
            else:
                return False
    except:
        return False
        


def get_filename(arg):
    '''
    get the input file name
    '''
    
    file_list = []
    for iarg in arg:
        if is_psrfits(iarg) == True:
            file_list.append(iarg)
            
    if len(file_list) >= 1:
        return file_list
    else:
        print("\n")
        print("-"*10)
        print("\n")
        print("Please Check whether the input file is psrfits \n")
        print("-"*10)
        Usage()
        print("\n")
        
        sys.exit(1)



def get_stt_subint(arg):
    '''
    get the start subint value from input argv
    '''
    if "-sttsub" in arg:
        sttsub_index = arg.index("-sttsub") + 1
    
        try:
            sttsub = int(arg[sttsub_index])
            return sttsub
        except:
            print("Please input the correct STT_SUBINT  value")
            print("-"*10)
            print("\n")
            Usage()
        
    else:
        print("Please input STT_SUBINT, e.g.: -sttsub 8")
        Usage()
        sys.exit(1)


def get_end_subint(arg):
    '''
    get the end subint value from input argv
    '''
    if "-endsub" in arg:
        endsub_index = arg.index("-endsub") + 1
        try:
            endsub = int(arg[endsub_index])
            return endsub
        except:
            print("Please input the correct END SUBINT  value")
            print("-"*10)
            print("\n")
            Usage()
    else:
        print("Please input  END SUBINT, e.g.: -endsub 16")
        Usage()
        sys.exit(1)




if __name__ == '__main__':
    arg = sys.argv
    print("------> Processing data  <-----\n")
    # open a PSRFITS  file
    file_list = get_filename(arg)
    stt_subint = get_stt_subint(arg)
    end_subint = get_end_subint(arg)
    for file in file_list:
        CutPsrfits(file,stt_subint,end_subint)
    print("------> Done <-----\n")







