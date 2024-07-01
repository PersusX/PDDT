#!/usr/bin/env python
import re
import sys
import math
import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.signal import argrelextrema
from scipy.signal import correlate
from Dedispersion import Dedispersion as dds
from FoldDataArchive import FoldFitsArchive as fda
from PsrFitsArchive import FITSArchive as pfa
import pyqtgraph as pg
import rec_rc
from DeRFI import channel_zap_median
from DeRFI import channel_zap_fft
#from pyqtgraph.Qt import QtGui
#from PyQt5 import QtCore
#from pyqtgraph.Qt import QtWidgets
from PyQt5.QtWidgets import QWidget,QLabel,QVBoxLayout,\
        QApplication,QCalendarWidget,QComboBox,QDateEdit,\
        QPushButton,QToolButton,QGridLayout,QHBoxLayout,\
        QSpacerItem,QSizePolicy,QLineEdit,QFrame,QDialog,\
        QFileDialog,QTextEdit,QWidget,QMainWindow,QCheckBox,\
        QSpinBox
from PyQt5.QtGui import QMovie,QIcon,QTextCharFormat,QPixmap,\
        QRegExpValidator,QTextCursor, QTextBlockFormat,QFont,QPen,\
        QColor,QTransform
from PyQt5.QtCore import Qt,QDate,QSize,QRect,QRegExp
from PyQt5.QtCore import QThread, pyqtSignal
class CSpinBox(QWidget):
    subintSignal = pyqtSignal(int)  # 将输入框的数据发送处出去
    style = '''
        QSpinBox {
            padding-top: 2px;
            padding-bottom: 2px;
            padding-left: 5px;
            padding-right: 5px;
            border:0px solid rgb(20,20,20);
            border-radius: 15px;
            color: rgb(200,200,200);
            selection-background-color: darkgray;
            background-color:rgba(20,20,20,0);
            selection-background-color: rgb(83,121,180);
            font-size: 16pt;
        }
        QLabel {
            color: rgb(255,255,255);
            font-size:35;
        }

        QSpinBox:hover {
            background-color: rgb(59,59,59);
        }


        QSpinBox::up-button {
            background-color: rgba(200,200,200,0);
            border: none;
            width: 15px;
            height: 10px;
            color: white;
        }

        QSpinBox::down-button {
            background-color: rgba(200,200,200,0);
            border: none;
            width: 15px;
            height: 10px;
            color: white;
        }

        QSpinBox::up-arrow {
            image: url(:icon/Resource/arrow-up-filling.svg);
        }
    
        QSpinBox::down-arrow {
            image: url(:icon/Resource/arrow-down-filling.svg);
        }
 
        QSpinBox::up-arrow:hover {
            image: url(:icon/Resource/arrow-up-filling-hover.svg);
        }
 
        QSpinBox::down-arrow:hover {
            image: url(:icon/Resource/arrow-down-filling-hover.svg);
        }

    '''



    def __init__(self, defaultsubint = '1', width = 200, heigh = 30, parent=None):
        super(CSpinBox, self).__init__(parent)


        self.label = QLabel("Select Subint number:")

        self.spin_box = QSpinBox()
        self.spin_box.setMinimum(1)  # 设置最小值
        self.spin_box.setMaximum(120)  # 设置最大值
        self.spin_box.setValue(1)  # 设置初始值
        self.spin_box.valueChanged.connect(self.on_spin_box_value_changed)
        if sys.platform == 'darwin':
            self.spin_box.setAttribute(Qt.WA_MacShowFocusRect, 0)

        # 创建水平布局，将上面两个控件添加进来
        self.hLayout = QHBoxLayout()
        self.hLayout.setContentsMargins(0, 0, 10, 0)
        self.hLayout.addWidget(self.label)
        self.hLayout.addWidget(self.spin_box)


        self.setAttribute(Qt.WA_StyledBackground)  # 很重要, 没有这句QSS对QWeight的有些代码会无效
        self.setLayout(self.hLayout)
        self.setFixedSize(width, heigh)
        self.setStyleSheet(self.style)


    def on_spin_box_value_changed(self, value):
        self.subintSignal.emit(value)  # 将输入框的数据发送出去

        print("here")
        #self.label.setText(f"Selected  Subint number: {value}")


class CSearchEdit(QWidget):
    style = '''
        #cSearchWeight{
            border:0px solid rgb(20,20,20);
            border-radius:15px;
            background-color:rgb(200,200,200);
        }

        #lineEdit{
            border:1px solid rgba(41, 57, 85,0);
            color:rgb(20,20,20);
            border-radius:15px;
            padding-left:10px;
            selection-background-color: darkgray;
            background-color:rgba(20,20,20,0);
            font-family:Helvetica;
            font-size:16px;
        }


        #searchBtn{
            border:0;
            background-color:rgba(200,200,200,0); /*完全透明*/
            border-image:url(:icon/Resource/search-star.svg)
        }

        #searchBtn:hover{
            border:0;
            background-color:rgba(200,200,200,0); /*完全透明*/
            border-image:url(:icon/Resource/search-star-hover.svg)
        }
    '''

    searchSignal = pyqtSignal(str)  # 将输入框的数据发送处出去

    def __init__(self, placeHolderText = '请输入PSRJ', width = 200, heigh = 30, parent=None):
        super(CSearchEdit, self).__init__(parent)

        # 先创建一个QLineEdit
        self.lineEdit = QLineEdit()
        self.lineEdit.setFrame(False)
 
        
        self.lineEdit.setObjectName('lineEdit')
        self.lineEdit.setPlaceholderText(placeHolderText)
        self.lineEdit.setMinimumSize(150, heigh)
        
        ##如果是系统是MacOs 消除QLineEdit focus 边框
        #print("os name: ",sys.platform)
        if sys.platform == 'darwin':
            self.lineEdit.setAttribute(Qt.WA_MacShowFocusRect, 0)
        # 创建一个搜索按钮
        self.searchBtn = QPushButton()
        self.searchBtn.setObjectName('searchBtn')
        self.searchBtn.setFixedSize(int(heigh*3/5), int(heigh*3/5))  # 正方形搜索按钮

        # 创建水平布局，将上面两个控件添加进来
        self.hLayout = QHBoxLayout()
        self.hLayout.setContentsMargins(0, 0, 10, 0)
        self.hLayout.addWidget(self.lineEdit)
        self.hLayout.addWidget(self.searchBtn)

        self.setObjectName('cSearchWeight')
        self.setAttribute(Qt.WA_StyledBackground)  # 很重要, 没有这句QSS对QWeight的有些代码会无效
        self.setLayout(self.hLayout)
        self.setFixedSize(width, heigh)
        self.setStyleSheet(self.style)
        '''绑定槽函数'''
        self.searchBtn.clicked.connect(self.btnSearch)
        self.lineEdit.returnPressed.connect(self.btnSearch)

    def btnSearch(self):
        self.searchSignal.emit(self.lineEdit.text())  # 将输入框的数据发送出去

class CustomROI(pg.ROI):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # 设置handle的样式
        
        #self.handleBrush = pg.mkBrush(color='r',width=4)
        self.handlePen = pg.mkPen(color='teal', width=5)
        self.scaleHandlePen=pg.mkPen(color='b', width=2)
        self.handleHoverPen = pg.mkPen(color='teal', width=2)
        self.handlePen.setWidth(5)
        #handleSi
        #rect = self.boundingRect()
        #self.handleSize = max(rect.width(), rect.height()) * 0.05
        for handle in self.handles:
            handle['item'].setPen(width=12)  # 设置线条宽度为12个像素


    def paint(self, p, *args):
        # 重写绘制方法，在绘制handle前先设置样式
        p.setPen(self.handlePen)
        super().paint(p, *args)



class ViewPulse(QMainWindow):
    style = '''
        
        
         QMainWindow{
            background-color: rgba(10,12,15,0.9);
        }
        QLabel {
            color: rgb(235,125,168);
            font-size:35;
        }
        QCheckBox{
            color:rgb(200,200,200);
            font=size:15px;
        
        }

        QPushButton{
            background-color: rgba(110,112,115,0.6);
            border: 2px solid;
            border-color: rgba(45,45,45,1);
            border-radius: 4px;
            padding: 4px;
            font-size: 15px;
            color:rgba(200,200,200,1);
            line-height:48px;
        }
        
        QPushButton:hover {
            background-color: rgb(90, 85, 152);
        }
        QPushButton:pressed {
            background-color: rgb(189, 147, 249);
            color: rgb(255, 255, 255);
        }



        QCheckBox{
            color:rgb(200,200,200);
            font=size:15px;
        
        }
        QCheckBox::indicator {
            border: 3px solid rgb(52, 59, 72);
            width: 15px;
            color: rgb(200, 200, 200);
            height: 15px;
            font-size: 15px;
            border-radius: 10px;
            background: rgb(44, 49, 60);
        }
        QCheckBox::indicator:hover {
            border: 3px solid rgb(58, 66, 81);
        }
        QCheckBox::indicator:checked {
            background: 3px solid rgb(52, 59, 72);
            border: 3px solid rgb(52, 59, 72);
            background-image: url(:icon/Resource/check.svg);
        }
    '''

    def __init__(self,filename):
        super().__init__()
        
        self.setWindowTitle("View pulse images {}".format(filename))
        ## 设置窗口背景颜色
        #self.setStyleSheet("background-color: rgba(50, 50, 50,0.9);")
        self.resize(1200,1100)
        #self.setGeometry(100, 100, 400, 300)
        # 创建中央部件和布局
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)
         # 设置布局到中央部件
        self.central_widget.setLayout(self.layout)
        self.arch = pfa(filename)
        self.npol = self.arch.sub_npol
        self.nsubint = self.arch.sub_nrows
        self.nchan = self.arch.sub_nchan
        self.dat_freq = self.arch.dat_freq[0]
        self.sub_chan_bw = self.arch.sub_chan_bw
        self.nsampsubint = self.arch.sub_nsblk
        self.tsamp = self.arch.sub_tbin
        self.arch.obsdata = self.arch.obsdata.squeeze()
        #print(self.arch.obsdata.shape)
        if self.nsubint ==1:
            new_shape = (1,) + self.arch.obsdata.shape
            self.arch.obsdata = self.arch.obsdata.reshape((new_shape))
        if self.npol >1:
            print("data shape:",self.arch.obsdata.shape)
            self.arch.obsdata = np.transpose(self.arch.obsdata, (0,2,3,1))
        else:
            self.arch.obsdata = np.transpose(self.arch.obsdata, (0,2,1))
        #self.dat_wts = self.arch.dat_wts
        self.downchanrate = 1
        self.downsamprate = 1
        
        self.downChanMode = False
        self.downSampMode = False
        self.dedispersion_state = False
        self.removeRFIMode = False
        self.removeZeroDmRFIMode = False
        self.dm = 0
        self.dds_nsubint = 1
        
        self.init_data()
        self.setPulseData()
        self.setupUI()
       
    def init_data(self):
        '''
        
        '''
        if self.npol > 1:
            self.init_data = self.arch.obsdata
            del self.arch.obsdata
            self.init_data = np.mean(self.init_data[:,:2,:,:], axis=1)
            # change data dimensions
            self.init_data = np.transpose(self.init_data, (1,0,2)).astype(np.float64)
        else:
            self.init_data = self.arch.obsdata
            del self.arch.obsdata
            #self.data = np.sum(self.data[:,:,:,:], axis=1)
            self.init_data = np.transpose(self.init_data, (1,0,2)).astype(np.float64)
        

        # scale and subbase with subint
        for isub in range(self.nsubint):
            
            pulse_profile = np.mean(self.init_data[:,isub,:],axis=0)
            offset_stt,offset_end = self.offset_determine(pulse_profile,self.nsampsubint//8)
            
            #for ichan in range(self.nchan):
                #minval = np.min(self.init_data[ichan,isub,:])
                #maxval = np.max(np.abs(self.init_data[ichan,isub,:]))
                #if minval==maxval:
                #    self.init_data[ichan,isub,:] = (self.init_data[ichan,isub,:] - minval) * 2 - 1
                #else:
                #    self.init_data[ichan,isub,:] = ((self.init_data[ichan,isub,:] - minval) / (maxval - minval)) * 2 - 1
                    
                #ichan_offset = np.mean(self.init_data[ichan,isub,offset_stt:offset_end],dtype=np.float64)
                
                #self.init_data[ichan,isub,:]  -= ichan_offset
 
        #sys.exit(1)
 
 

    def downsample_chan(self, arr, down_rate):
        '''
        the input array is downsampled by the specified dimension, and the average is calculated every eight values.
        para:
        arr (numpy.ndarray)：The input array to downsample.
        dim_to_downsample (int)：The index of the dimension to downsample.
        return:
        numpy.ndarray：The downsampled array.
        '''
        # Gets the shape of the input array
        original_shape = arr.shape
        dim_to_downsample = 0
        # Calculate the shape after downsampling
        downsampled_shape = list(original_shape)
        downsampled_shape[dim_to_downsample] = original_shape[dim_to_downsample] // down_rate
        # Initialize the downsampled array
        downsampled_data = np.zeros(downsampled_shape)
        # Iterate over the original array and compute the average of every eight values
        for i in range(0, original_shape[dim_to_downsample], down_rate):
            slice_start = [slice(None)] * len(original_shape)
            slice_end = [slice(None)] * len(original_shape)
            # Slices on the selected dimension
            slice_start[dim_to_downsample] = slice(i, i + down_rate)  # slice using the slice object
            slice_end[dim_to_downsample] = slice(i, i + down_rate)
            # Calculating the average
            downsampled_data[i // down_rate] = np.mean(arr[tuple(slice_start)], axis=dim_to_downsample)
        return downsampled_data


    def downsamp_time(self, data, down_rate):
        '''
        Do downsampling on the data!
        '''
        # Remainder number between data length and
        dat_rem = -(data.shape[-1] % down_rate)
        if dat_rem !=0: data = data[:,:dat_rem]
        # Get shape of input array
        original_shape = data.shape
        # Calculate the shape after downsampling
        downsampled_shape = (*original_shape[:1], original_shape[1] // down_rate)
        # Initialize the downsampled array
        downsampled_data = np.zeros(downsampled_shape)
        # Iterate over the original array and compute the average of every eight values
        for j in range(0, original_shape[1], down_rate):
            slice_start = [slice(None)] * len(original_shape)
            slice_end = [slice(None)] * len(original_shape)
            # Slice on the second dimension
            slice_start[1] = j
            slice_end[1] = j + down_rate
            # Calculating the average
            downsampled_data[:, j // down_rate] = np.mean(data[:, slice_start[1]:slice_end[1]], axis=1)
        return downsampled_data


    def has_value_greater_than(self, arr, threshold):
        
        for element in arr:
            if element > threshold:
                return True
        return False

    def remove_zeroDMRFI(self,data,num_chan):
        '''
        remove the RFI of DM = 0
        '''
        sys.stdout.flush()
        oldpcnt = ""
        pulse_profile = np.mean(data,axis=0)
        data_mean = np.mean(data, axis=1)
        data_std = np.std(data, axis=1)
        data_shape = data.shape
        nsamp = len(pulse_profile)
        plot_result = False
        thd_para =  np.mean(pulse_profile)+ 1*np.std(pulse_profile)
        
        if self.has_value_greater_than(pulse_profile, thd_para):
            #pulse_profile[pulse_profile < thd_para] = 0
            #pulse_profile[pulse_profile > thd_para] = 1
            signal_rfi = np.where(pulse_profile < thd_para, 0, 1)
            for ichan in range(num_chan):
                pcnt = "%d" %(100*(ichan+1)/num_chan)
                if pcnt != oldpcnt:
                    sys.stdout.write("% 4s%% remove zero-dm RFI complete\r" % pcnt)
                    sys.stdout.flush()
                #thd_chan = np.mean(data[ichan]) + 2*np.std(data[ichan])
                thd_chan = data_mean[ichan] + 1 * data_std[ichan]
                if self.has_value_greater_than(data[ichan], thd_chan): pass
                chan_data = np.where(data[ichan] < thd_chan, 0, 1)
                mtp_chanpro_pro = np.multiply(chan_data,signal_rfi)
                mtp_chanpro_pro = mtp_chanpro_pro.astype(int)
                mean_val = np.mean(data[ichan]) + 1*np.random.normal(np.mean(data[ichan]),np.std(data[ichan]))
                mean_val = int(mean_val)
                data[ichan][mtp_chanpro_pro == 1] = mean_val
            
            return data
        else:
            return data


    def scale_data(self):
        '''
        test sacle data
        '''


    def setPulseData(self,pulse_index=0, averagemode=True):
        '''
        pulse_index: int ; pulse index, the first pulse index is 0
        average mode: if True , set average pulse
        '''
        if averagemode:
            #print(self.arch.obsdata.shape,np.mean(self.arch.obsdata,axis=0)[:2].shape)
            
            #if self.npol > 1:
            #    self.data = np.sum(np.mean(self.arch.obsdata,axis=0)[:2],axis=0)
            #else:
            self.data = np.mean(self.init_data,axis=1)
            
            
            self.pulse_profile = np.mean(self.data,axis=0)
            self.offset_start,self.offset_end = self.offset_determine(self.pulse_profile, self.nsampsubint//8)
            #print("offset region:",self.offset_start,self.offset_end)
            #print(self.dat_wts.shape)
            #self.dat_wts = np.prod(self.dat_wts,axis=0)
            self.dat_wts = np.mean(self.arch.dat_wts,axis=0)
            
            if self.removeRFIMode:
                # 中值滤波 如果有频带被标记为干扰，则进行消除
                spectrum = np.sum(self.data,axis=1)
                #set parameter for removing RFI
                rms_threshold = 3.0
                window_size =21
                mask = channel_zap_median(spectrum, rms_threshold, window_size).mask
                zap_chnum = [num for num in mask if num  == True]
                print("ChannelZapMedian::weight zapped {} channels".format(len(zap_chnum)))
            else:
                mask = [False]*self.nchan
            
            if self.removeZeroDmRFIMode:
                self.data = self.remove_zeroDMRFI(self.data,self.nchan)
            
            for ichan in range(self.nchan):
                if mask[ichan]:
                    self.data[ichan] =0
                    #print("zap:",i)
                else:
                    ichan_offset = np.mean(self.data[ichan][self.offset_start:self.offset_end])
                    self.data[ichan] -= ichan_offset
            try:
                if self.dedispersion_checkbox.isChecked():
                    self.dedispersion(self.dm)
                    self.dedispersion_state = True
                else:
                    self.dedispersion_state = False
            except:
                pass
           
           
            self.pulse_profile = np.mean(self.data,axis=0)
            self.set_SNR()
            
            try:
                self.setSnrLabel()
            except:
                pass
            
            try:
                self.img.setImage(self.data)
                ymin = np.min(self.dat_freq)
                ymax= np.max(self.dat_freq)
                self.p1.setLimits(xMin=0., yMin=ymin, xMax=1.0, yMax=ymax)
                #self.img.setTransform(self.tr.scale(1/self.nsampsubint,\
                        # self.sub_chan_bw).translate(0,1000/self.sub_chan_bw))
                percentile_lower = 2
                # Upper quantile (e.g., upper 95%)
                percentile_upper = 98.5
                vmin = np.percentile(self.data.flatten(), percentile_lower)
                vmax = np.percentile(self.data.flatten(), percentile_upper)
                #vmin, vmax = np.mean(self.data)-3*np.std(self.data), np.mean(self.data)+3*np.std(self.data)
                self.hist.setLevels(vmin, vmax)
                self.updatePlot()
                #print("Sucess set img")
                #self.img.setImage(self.data)
                #vmin, vmax = np.mean(self.data)-np.std(self.data), np.mean(self.data)+np.std(self.data)
                #self.p1.setLimits(xMin=0., yMin=ymin, xMax=1.0, yMax=ymax)
                #self.hist.setLevels(vmin, vmax)
                #self.updatePlot()
                #print("scess")
            except:
                #print("Something error")
                pass
            
        else:
            if not self.singlepulse_checkbox.isChecked():return
            #print("number nsubint for dds:",self.dds_nsubint)
            if self.dds_nsubint >1:
                # 设置数据范围
                pulse_end_index = pulse_index+self.dds_nsubint
                self.data = self.init_data[:,pulse_index:pulse_end_index,:]
                self.data = self.data.reshape(self.nchan,-1)
            else:
                self.data = self.init_data[:,pulse_index,:]
            
            #print("Deal data shape:",self.data.shape)
            self.data = self.data.astype(np.float64)

            if self.dedispersion_checkbox.isChecked():
                self.dedispersion(self.dm)
                self.dedispersion_state = True
                #print("Sucess dds!")
            else:
                self.dedispersion_state = False

            if self.downSampMode:
                self.data = self.downsamp_time(self.data, self.downsamprate)
                #print("downsamp", self.downsamprate)

            self.pulse_profile = np.mean(self.data,axis=0)
            #print("nsample",self.data.shape[1])
            self.offset_start,self.offset_end = self.offset_determine(self.pulse_profile, self.data.shape[1]//8)
            
            #print("Pulese index:", self.pulse_index)
            #print("offset region:",self.offset_start,self.offset_end)
            #print(self.dat_wts.shape)
            #self.dat_wts = np.prod(self.dat_wts,axis=0)
            self.dat_wts = self.arch.dat_wts[pulse_index]
            #print(self.dat_wts.shape)
            
            if self.removeRFIMode:
                # 中值滤波 如果有频带被标记为干扰，则进行消除
                spectrum = np.sum(self.data,axis=1)
                #set parameter for removing RFI
                #print("Spectrum val:",np.mean(spectrum),np.std(spectrum))
                #plt.plot(spectrum)
                #plt.show()
                rms_threshold = 2.5 # * (np.std(spectrum)/np.mean(spectrum))
                window_size = self.nchan//64
                print(f"rms_thd:{rms_threshold}, window_size:{window_size}")
                mask = channel_zap_median(spectrum, rms_threshold, window_size).mask
                zap_chnum = [num for num in mask if num  == True]
                print("ChannelZapMedian::weight zapped {} channels".format(len(zap_chnum)))
            else:
                mask = [False]*self.nchan

            for ichan in range(self.nchan):
                if  mask[ichan]:
                    self.data[ichan] = 0
                else:
                    minval = np.min(self.data[ichan])
                    maxval = np.max(np.abs(self.data[ichan]))
                    if minval==maxval:
                        self.data[ichan] = (self.data[ichan] - minval) * 2 - 1
                    else:
                        self.data[ichan] = ((self.data[ichan] - minval) / (maxval - minval)) * 2 - 1
                    
                    ichan_offset = np.mean(self.data[ichan][self.offset_start:self.offset_end])
                    self.data[ichan] -= ichan_offset
                    
            if self.removeZeroDmRFIMode:
                self.data = self.remove_zeroDMRFI(self.data,self.nchan)

            # set pulse profile
            self.pulse_profile = np.mean(self.data,axis=0)
            self.set_SNR()
            self.setSnrLabel()
            
            try:
                self.img.setImage(self.data)
                ymin = np.min(self.dat_freq)
                ymax= np.max(self.dat_freq)
                self.p1.setLimits(xMin=0., yMin=ymin, xMax=1.0*self.dds_nsubint/self.downsamprate, yMax=ymax)
                #self.img.setTransform(self.tr.scale(1/self.nsampsubint,\
                        # self.sub_chan_bw).translate(0,1000/self.sub_chan_bw))
                percentile_lower = 2
                # Upper quantile (e.g., upper 95%)
                percentile_upper = 98.5
                vmin = np.percentile(self.data.flatten(), percentile_lower)
                vmax = np.percentile(self.data.flatten(), percentile_upper)
                
                self.hist.setLevels(vmin, vmax)
                self.updatePlot()
                #print("Sucess set img")
            except:
                print("unscess")
                pass
    
    def setupUI(self):
        # set background color
        #pg.setConfigOption('background', 'w')

        hly = QHBoxLayout()
        self.averagepulse_checkbox = QCheckBox("Average Subint")
        self.averagepulse_checkbox.setChecked(True)  # 设置复选框默认为选中状态
        self.averagepulse_checkbox.stateChanged.connect(self.setAveragePulseMode)
        hly.addWidget(self.averagepulse_checkbox)

        self.singlepulse_checkbox = QCheckBox("Single Subint:")
        self.singlepulse_checkbox.setChecked(False)  # 设置复选框默认为选中状态
        self.singlepulse_checkbox.stateChanged.connect(self.setSingelePulseMode)


        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.singlepulse_checkbox)

        self.searchPulse = CSearchEdit('Input pulse index', 180)
        
        self.searchPulse.searchBtn.setEnabled(False)
        self.searchPulse.searchSignal.connect(self.setPulse)
        
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.searchPulse)
                
        # 创建QPushButton
        
        self.prevpulse_button = QPushButton('Prev')
        self.prevpulse_button.setFixedSize(100,30)
        
        self.prevpulse_button.clicked.connect(self.prevpulseClicked)
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.prevpulse_button)
        
        self.nextpulse_button = QPushButton('Next')
        self.nextpulse_button.setFixedSize(100,30)
        self.nextpulse_button.clicked.connect(self.nextpulseClicked)

        # default  average pulse mode, the next and previous pulse buttons cannot be clicked
        self.nextpulse_button.setEnabled(False)
        self.prevpulse_button.setEnabled(False)
        spacer = QSpacerItem(5, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.nextpulse_button)
        
        self.pulseshowlabel = QLabel()
        self.pulseshowlabel.setText("Current pulse: Average Pulse")
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.pulseshowlabel)

        #self.snr = np.sum(self.pulse_profile)
        #self.set_SNR()
        self.snrLabel = QLabel()
        self.snrLabel.setText("SNR: {:.2f} Peak SNR:{:.2f}".format(self.snr,self.peak_snr))
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.snrLabel)

        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        hly.addItem(spacerItem)

        # 将第一行控件添加到布局中
        self.layout.addLayout(hly)
        
        # 第二排控件
        hly = QHBoxLayout()
        
        self.dedispersion_checkbox = QCheckBox("De-Dispersion:")
        self.dedispersion_checkbox.setChecked(True)  # 设置复选框默认为选中状态
        #self.dedispersion_checkbox.setContentsMargins(300,0,30,0)
        self.dedispersion_checkbox.stateChanged.connect(self.setDedispersionMode)
        hly.addWidget(self.dedispersion_checkbox)
        
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        self.dedispersion_Edit = CSearchEdit('Input DM', 120)
        
        self.dedispersion_Edit.searchBtn.setEnabled(True)
        self.dedispersion_Edit.searchSignal.connect(self.setDisersionMeasure)
        #self.dedispersion_Edit.setContentsMargins(0,0,0,0)
        hly.addWidget(self.dedispersion_Edit)
  
        self.subint_Spinbox = CSpinBox(1, 230,50)
        self.subint_Spinbox.subintSignal.connect(self.setDdsNsubint)
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.subint_Spinbox)
  
        self.reRFI_checkbox = QCheckBox("Remove-ChannelRFI")
        self.reRFI_checkbox.setChecked(False)  # 设置复选框默认为非选中状态
        self.reRFI_checkbox.stateChanged.connect(self.setReRFIMode)
        
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.reRFI_checkbox)
  

        self.rezerodmRFI_checkbox = QCheckBox("Remove-ZeroDmRFI")
        self.rezerodmRFI_checkbox.setChecked(False)  # 设置复选框默认为非选中状态
        self.rezerodmRFI_checkbox.stateChanged.connect(self.setReZeroDmRFIMode)
        
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        hly.addWidget(self.rezerodmRFI_checkbox)


  
        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        hly.addItem(spacerItem)
        # 设置左对齐
        hly.setAlignment(Qt.AlignLeft)
        # 将第二行控件添加到布局中
        self.layout.addLayout(hly)
        
        # 第三排控件
        hly = QHBoxLayout()
        
        self.downchan_checkbox = QCheckBox("Down channel:")
        self.downchan_checkbox.setChecked(True)  # 设置复选框默认为选中状态
        #self.dedispersion_checkbox.setContentsMargins(300,0,30,0)
        self.downchan_checkbox.stateChanged.connect(self.setDownchanMode)
        hly.addWidget(self.downchan_checkbox)
        
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        self.downchan_Edit = CSearchEdit('DownRate ', 120)
        
        self.downchan_Edit.searchBtn.setEnabled(True)
        self.downchan_Edit.searchSignal.connect(self.setDownchanrate)
        #self.dedispersion_Edit.setContentsMargins(0,0,0,0)
        hly.addWidget(self.downchan_Edit)
        


        self.downsamp_checkbox = QCheckBox("Down time:")
        self.downsamp_checkbox.setChecked(False)  # 设置复选框默认为选中状态
        #self.dedispersion_checkbox.setContentsMargins(300,0,30,0)
        self.downsamp_checkbox.stateChanged.connect(self.setDownsampMode)
        hly.addWidget(self.downsamp_checkbox)
        
        spacer = QSpacerItem(20, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        hly.addItem(spacer)  # 添加间距
        self.downsamp_Edit = CSearchEdit('DownRate ', 120)
        
        self.downsamp_Edit.searchBtn.setEnabled(True)
        self.downsamp_Edit.searchSignal.connect(self.setDownsamprate)
        #self.dedispersion_Edit.setContentsMargins(0,0,0,0)
        hly.addWidget(self.downsamp_Edit)
        
        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        hly.addItem(spacerItem)
        # 设置左对齐
        hly.setAlignment(Qt.AlignLeft)
        # 将第二行控件添加到布局中
        self.layout.addLayout(hly)
        
        #win.addLayout(layout)
        plot_widget = pg.GraphicsLayoutWidget()
        self.layout.addWidget(plot_widget)

        # Interpret image data as row-major instead of col-major
        pg.setConfigOptions(imageAxisOrder='row-major')

        # A plot area (ViewBox + axes) for displaying the image
        self.p1 = plot_widget.addPlot(title="Phase vs Frequency")

        self.p1.showAxes( True, showValues=(True, False, False, True) )
        # 设置 动态铺图片与坐标轴贴合
        self.p1.setDefaultPadding(0.0)
        self.p1.getAxis('bottom').setLabel("Pulse Phase", units="")
        self.p1.getAxis('left').setLabel("Frequency (MHz)", units="")

        # Item for displaying image data
        self.img = pg.ImageItem()
        self.p1.addItem(self.img)

        yrange= np.max(self.dat_freq) - np.min(self.dat_freq)

        # Custom ROI for selecting an image region
        self.roi = CustomROI([0, np.min(self.dat_freq)], [1, yrange])

        # 创建自定义的handle样式
        handle_pen = QPen(QColor('red'))
        handle_pen.setWidth(5)

        # 在 四个边框 添加 控制 handle
        self.roi.addScaleHandle([0, 0.5], [1., 0.5])
        self.roi.addScaleHandle([1, 0.5], [0, 0.5])
        self.roi.addScaleHandle([0.5, 0.], [0.5, 1.])
        self.roi.addScaleHandle([0.5, 1.], [0.5, 0.])
        # 在 四个边角 添加 控制 handle
        self.roi.addScaleHandle([1, 1], [0., 0.])
        self.roi.addScaleHandle([0, 0], [1., 1.])
        self.roi.addScaleHandle([1, 0], [0., 1.])
        self.roi.addScaleHandle([0, 1], [1., 0.])

        roipen = pg.mkPen(color='#214761',width=2)
        self.roi.setPen(roipen)
        self.roi.hoverPen = pg.mkPen(color='teal',width=3)
        self.p1.addItem(self.roi)
        self.roi.setZValue(10)  # make sure ROI is drawn above image


        # Contrast/color control
        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img)
        plot_widget.addItem(self.hist)
        self.hist.vb.setMouseEnabled(y=False) # makes user interaction a little easier
        # Another plot area for displaying ROI data
        plot_widget.nextRow()
        self.p2 = plot_widget.addPlot(colspan=1)
        self.p2.setMaximumHeight(250)
        # 设置数据
        self.img.setImage(self.data)
        
        # 设置缩放与坐标轴范围
        ymin = np.min(self.dat_freq)
        ymax= np.max(self.dat_freq)
        
        self.p1.setLimits(xMin=0., yMin=ymin, xMax=1.0, yMax=ymax)
        self.img.setOpts(autoRange=False)
        self.img.setAutoDownsample(False)

        # 将回调函数连接到窗口的resize事件
        # 创建自定义颜色映射
        color_map = pg.ColorMap(
            pos=np.array([0.0, 0.5, 0.5, 1.0]),
            color=np.array([[255, 255, 0, 255], [255, 0, 0, 255], [255, 255, 0, 255], [255, 255, 255, 255]]))

        # 获取颜色映射的查找表
        self.lut = color_map.getLookupTable(nPts=256)

        percentile_lower = 2
        # Upper quantile (e.g., upper 95%)
        percentile_upper = 98.5
        vmin = np.percentile(self.data.flatten(), percentile_lower)
        vmax = np.percentile(self.data.flatten(), percentile_upper)
        #vmin, vmax = np.mean(self.data)-3*np.std(self.data), np.mean(self.data)+3*np.std(self.data)
        self.resizeEvent = self.onResize
        self.show()

        self.img.setLookupTable(self.lut)
        ## 设置初始 的 灰度范围值
        self.hist.setLevels(vmin, vmax)
        # set position and scale of image
        # 缩放，位置，以设置X、Y坐标系
        self.tr = QTransform()
        if self.sub_chan_bw<0:
            self.img.setTransform(self.tr.scale(1/self.nsampsubint, \
                    -self.sub_chan_bw).translate(0,min(self.dat_freq)/-self.sub_chan_bw))
        else:
             self.img.setTransform(self.tr.scale(1/self.nsampsubint, \
                    self.sub_chan_bw).translate(0,min(self.dat_freq)/self.sub_chan_bw))
        
        # 所选区域的脉冲轮廓
        # 创建一个PlotCurveItem对象，用于画 pulse profile
        self.curve_index = 0
        self.pulse_curve = pg.PlotCurveItem()
        pen = pg.mkPen(color='black', width=2)
        pen = pg.mkPen(color=(255, 0, 0, 245), width=2)  # 设置颜色的RGBA值，其中245表示50%的透明度
        self.pulse_curve.setPen(pen)
        self.p2.addItem(self.pulse_curve)

        # 区域更改 连接 plot 脉冲轮廓
        #self.roi.sigRegionChanged.connect(self.updatePlot)
        # 区域更改结束 连接 plot 脉冲轮廓
        self.roi.sigRegionChangeFinished.connect(self.updatePlot)
        
        self.updatePlot()

        self.img.hoverEvent = self.imageHoverEvent
        self.setStyleSheet(self.style)
        
        # 显示窗口
        #app.exec_()
        
        # 定义窗口大小变化的回调函数
    def onResize(self,event):
        # 创建自定义颜色映射
        self.img.setLookupTable(self.lut)

    
    def setSnrLabel(self):
        '''
        set SNR vale show in window
        '''
        self.snrLabel.setText("SNR: {:.2f} Peak SNR: {:.2f}".format(self.snr,self.peak_snr))
    

    def updatePlot(self):
    
        if self.curve_index != 0:
            self.pulse_curve.clear()
        
        selected = self.roi.getArrayRegion(self.data, self.img)
        # 设置 数据
        self.pulse_curve.setData(selected.mean(axis=0))
        self.curve_index += 1
    

    def imageHoverEvent(self,event):
        """
        Show the position, pixel, and value under the mouse cursor.
        """
        if event.isExit():
            self.p1.setTitle("")
            return
        pos = event.pos()
        i, j = pos.y(), pos.x()
        i = int(np.clip(i, 0, self.data.shape[0] - 1))
        j = int(np.clip(j, 0, self.data.shape[1] - 1))
        val = self.data[i, j]
        ppos = self.img.mapToParent(pos)
        x, y = ppos.x(), ppos.y()
        self.p1.setTitle("pos: (%0.1f, %0.1f)  pixel: (%d, %d)  value: %.3g" % (x, y, i, j, val))


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
        #print("the reference frequency:",freq_ref,self.arch.sub_tbin)
        #According to the given dm value, The time correction for each frequency channel is calculated
        self.t_delays = (4.148808 * self.dm * (self.arch.dat_freq[0] ** -2 - freq_ref ** -2) * 1e3 / self.arch.sub_tbin).astype(np.int64)
        # the delay time between the highest channel center frequency and the lowest center frequency
        self.max_delay = np.max(self.t_delays)


    def setDisersionMeasure(self,dm):
        '''
        Dedispersion is performed on the data
        '''
        # 设置 dm
        self.dm = float(dm)
        self.delay_t()
        
        # 设置数据
        if not self.dedispersion_checkbox.isChecked():return
        try:
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
        except:
            self.setPulseData(0, averagemode=True)
       
        #for ichan in range(self.nchan):
            # Refer to the delay time of the highest frequency
        #    ifreq_delay = self.t_delays[ichan]
        #    self.data[ichan] = np.roll(self.data[ichan],-ifreq_delay)
        #self.img.setImage(self.data)
        #vmin, vmax = np.sort(self.data.flatten())[int(self.data.shape[0]*
        #                self.data.shape[1]/50)], np.sort(self.data.flatten())[int
        #                    (self.data.shape[0]*self.data.shape[1]/50*49)]
        #vmin, vmax = np.mean(self.data)-3*np.std(self.data), np.mean(self.data)+3*np.std(self.data)
        #self.hist.setLevels(vmin, vmax)
        #self.updatePlot()
        

    def dedispersion(self,dm):
        '''
        Dedispersion is performed on the data
        '''
        # 设置 dm
        self.dm = float(dm)
        self.delay_t()
        for ichan in range(self.nchan):
            # Refer to the delay time of the highest frequency
            ifreq_delay = self.t_delays[ichan]
            self.data[ichan] = np.roll(self.data[ichan],-ifreq_delay)



    def setDedispersionMode(self):
        
        if self.dedispersion_checkbox.isChecked():
            if not self.dedispersion_state:
                if self.singlepulse_checkbox.isChecked():
                    self.setPulseData(self.pulse_index, averagemode=False)
                    self.dedispersion_state = True
                else:
                    self.setPulseData(0, averagemode=True)
                    self.dedispersion_state = True
        else:
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
                self.dedispersion_state = False
            else:
                self.setPulseData(0, averagemode=True)
                self.dedispersion_state = False
    
    
    def setReRFIMode(self):
        
        if self.reRFI_checkbox.isChecked():
            self.removeRFIMode = True
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
        else:
            self.removeRFIMode = False
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
            
    def setReZeroDmRFIMode(self):

        if self.rezerodmRFI_checkbox.isChecked():
            self.removeZeroDmRFIMode = True
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
        else:
            self.removeZeroDmRFIMode = False
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
        

    def setSingelePulseMode(self):
        '''
        set Pulse mode: single pulse mode
        '''
        if self.singlepulse_checkbox.isChecked():
            self.pulse_index = 0
            self.averagepulse_checkbox.setChecked(False)
            self.setPulseData(self.pulse_index,averagemode=False)
            self.searchPulse.searchBtn.setEnabled(True)
            self.nextpulse_button.setEnabled(True)
            self.prevpulse_button.setEnabled(True)
            
        else:
            self.averagepulse_checkbox.setChecked(True)
            self.setPulseData(0,averagemode=True)
            # if average pulse mode, the next and previous pulse buttons cannot be clicked
            self.searchPulse.searchBtn.setEnabled(False)
            self.nextpulse_button.setEnabled(False)
            self.prevpulse_button.setEnabled(False)
            

    def setAveragePulseMode(self):
        '''
        set Pulse mode: average pulse mode
        '''
        if self.averagepulse_checkbox.isChecked():
            self.singlepulse_checkbox.setChecked(False)
            self.setPulseData(0,averagemode=True)
            # if average pulse mode, the next and previous pulse buttons cannot be clicked
            self.searchPulse.searchBtn.setEnabled(False)
            self.nextpulse_button.setEnabled(False)
            self.prevpulse_button.setEnabled(False)
            
        else:
            self.singlepulse_checkbox.setChecked(True)
            self.searchPulse.searchBtn.setEnabled(True)
            self.nextpulse_button.setEnabled(True)
            self.prevpulse_button.setEnabled(True)


    def setDownchanMode(self):
        '''
        set down channel mode rate
        '''
        
        if self.downchan_checkbox.isChecked():
            self.downChanMode = True
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
        else:
            self.downChanMode = False
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
    
    def setDownchanrate(self,downrate):
        '''
        set down channel mode rate
        '''
        self.downchanrate = int(downrate)
        
        if not self.dedispersion_checkbox.isChecked():return
        try:
            if self.singlepulse_checkbox.isChecked():
                self.setPulseData(self.pulse_index, averagemode=False)
            else:
                self.setPulseData(0, averagemode=True)
        except:
            self.setPulseData(0, averagemode=True)



    def setDownsampMode(self):
        '''
        set down samp mode rate
        '''
        
        if self.downsamp_checkbox.isChecked():
            self.downSampMode = True
        else:
            self.downSampMode = False
           


    def setDownsamprate(self,downrate):
        '''
        set down samp mode rate
        '''
        self.downsamprate = int(downrate)



    def average_filter(self, signal, window_size):
        """
            平均滤波函数

        参数:
            signal (list): 输入信号列表
            window_size (int): 滤波窗口大小

        返回:
            filtered_signal (list): 平均滤波后的信号列表
        """
        filtered_signal = []
        half_window = window_size // 2

        for i in range(len(signal)):
            start = max(0, i - half_window)
            end = min(len(signal), i + half_window + 1)
            window = signal[start:end]
            average = sum(window) / len(window)
            filtered_signal.append(average)

        return filtered_signal


    def least_squares_fit(self , x, y):
        '''
        Fit a line using least squares
        '''
    
        n = len(x)
        # 计算x和y的均值
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        # 计算最小二乘法的系数
        numerator = np.sum((x - x_mean) * (y - y_mean))
        denominator = np.sum((x - x_mean) ** 2)
        slope = numerator / denominator
        intercept = y_mean - slope * x_mean
        return [slope, intercept]


    def fit_line(self, x, y):
        # 使用最小二乘法拟合直线
        params = self.least_squares_fit(x, y)
        return params

    def offset_determine(self, pulse_signal, window_size):
        
        offset = []
        #print(pulse_signal, window_size)
        num_window = len(pulse_signal) // window_size
        remainder = len(pulse_signal) % window_size
        if remainder != 0: num_window += 1
        slope_window =[]
        i = 0
        while i < num_window-1:
            start = i * window_size
            end = (i+1) * window_size
            if end > len(pulse_signal):end = len(pulse_signal)
            window = pulse_signal[start:end]
            x = np.arange(len(window))
            para = self.fit_line(x,window)
            slope_window.append(para[0])
            #print(para[0])
            phase = np.arange(start,end)
            i+=1
            #if para[0] > -0.005 and para[0] < 0.005:
            offset.extend(window)
        #print(np.fabs(slope_window))
        if len(offset) ==0:
            return start,end
        else:
            # The minimum slope window is taken as the offset.
            offset_index = np.argmin(np.fabs(slope_window))
            #print("offset index:",offset_index)
            start = offset_index * window_size
            end = (offset_index+1) * window_size
            if end > len(pulse_signal):end = len(pulse_signal)
            offset = pulse_signal[start:end]
        return start,end



    # 定义按钮的点击事件处理函数
    def nextpulseClicked(self):
        
        if self.pulse_index == self.nsubint - 1:
            return
        else:
            self.pulse_index += 1
        
        if self.pulse_index == self.nsubint - 1:
            self.pulseshowlabel.setText("Current Subint: %s last subint"%self.pulse_index)
        else:
            self.pulseshowlabel.setText("Current Subint: %s "%self.pulse_index)
        self.setPulseData(self.pulse_index,averagemode=False)
    

    def prevpulseClicked(self):
        
        if self.pulse_index == 0:
            return
        else:
            self.pulse_index -= 1
        
        if self.pulse_index == 0:
            self.pulseshowlabel.setText("Current Subint: %s (first subint)"%self.pulse_index)
        else:
            self.pulseshowlabel.setText("Current Subint: %s "%self.pulse_index)

        self.setPulseData(self.pulse_index,averagemode=False)


    def setDdsNsubint(self,dds_nsubint):
        self.dds_nsubint = dds_nsubint
        print("Dedispersion nsubint",self.dds_nsubint)
        

    def setPulse(self,index):
    
        index = int(index)

        if index >= 0 and index < self.nsubint:
            self.pulse_index = index
        else:
            return
        
        # set pulse data
        self.setPulseData(self.pulse_index,averagemode=False)
        
        # set pulse index label
        if self.pulse_index == self.nsubint - 1:
            self.pulseshowlabel.setText("Current Subint: %s last Subint"%self.pulse_index)
        else:
            self.pulseshowlabel.setText("Current Subint: %s "%self.pulse_index)

        if self.pulse_index == 0:
            self.pulseshowlabel.setText("Current Subint: %s (first Subint)"%self.pulse_index)
        else:
            self.pulseshowlabel.setText("Current Subint: %s "%self.pulse_index)
        
        #self.snr = np.sum(self.pulse_profile)/len(self.pulse_profile)
        self.set_SNR()
        # set snr label
        self.snrLabel.setText("SNR: {:.2f} Peak SNR:{:.2f}".format(self.snr,self.peak_snr))


    def pulse_template(self, amplitude,position,width,length):
        x = np.arange(length)
        template = amplitude * np.exp(-(x - position)**2 / (2 * width**2))
        return template
        
        

    def find_pulse_position(self,pulse_data, template,temp_poi):
        '''
        Using the cross correlation to obtaine pulse center position
        '''
        # 计算脉冲数据和模板的互相关
        correlation = np.correlate(pulse_data, template, mode='same')
        # 寻找互相关的峰值位置
        pulse_position =  temp_poi + np.argmax(correlation) - (len(template)//2)
        return pulse_position



    def find_pulse_range(self, noisy_pulse_data, pulse_template,temp_poi):
        '''
        Using the cross correlation to obtain the pulse range
        '''
        num_samp = len(noisy_pulse_data)
        t = np.arange(len(noisy_pulse_data))
        noisy_pulse_data = np.array(noisy_pulse_data)
        pulse_template = np.array(pulse_template)
        correlation = correlate(noisy_pulse_data, pulse_template, mode='same')
        # Find the index of maximum correlation
        max_correlation_index = np.argmax(correlation)
        # Calculate estimated pulse width (full width at half maximum, FWHM)
        half_max = correlation[max_correlation_index] / 2
    
        left_idx_candidates = np.where(correlation[:max_correlation_index] <= half_max)[0]
        if left_idx_candidates.size>0:
            left_idx = left_idx_candidates[-1]
        else:
            left_idx = temp_poi
      
        est_position = self.find_pulse_position(noisy_pulse_data, pulse_template,temp_poi)
        right_idx_candidates = np.where(correlation[max_correlation_index:] <= half_max)[0] + max_correlation_index
    
        if right_idx_candidates.size > 0:
            right_idx = right_idx_candidates[0]
        else:
            right_idx = est_position - num_samp//10
        #if left_idx_candidates.size == 0 or right_idx_candidates.size == 0:
        #    print("error:",left_idx,right_idx,est_position)
    
        if left_idx_candidates.size == 0 and right_idx_candidates.size == 0:
            left_idx = est_position - num_samp //10
            right_idx = est_position + num_samp //10
    
    
        elif left_idx_candidates.size == 0 and right_idx_candidates.size !=0:
            half_width = right_idx - est_position
            left_idx = est_position - half_width
        
        elif left_idx_candidates.size != 0 and right_idx_candidates.size ==0:
            half_width = est_position -left_idx
            right_idx = est_position + half_width
        else:
            pass
            #print("undo anything...")
        
    
        if left_idx < 0: left_idx=0
        if right_idx >= num_samp: right_idx= num_samp-1
        if (abs(right_idx-est_position) > 3 * abs(left_idx-est_position)) or \
                (abs(right_idx-est_position) <  abs(left_idx-est_position))/3:
            half_width = min(abs(right_idx-est_position), abs(left_idx-est_position))
            left_idx = est_position - half_width
            right_idx = est_position + half_width
        
        if left_idx == right_idx:
        
            if est_position > num_samp //20 and est_position < num_samp *19//20:
                left_idx = est_position - num_samp //20
                right_idx = est_position + num_samp //20
            else:
                est_position = np.argmax(noisy_pulse_data)
                left_idx = est_position - num_samp //20
                right_idx = est_position + num_samp //20
                
            
        if left_idx < 0: left_idx=0
        if right_idx >= num_samp: right_idx= num_samp-1
        
        estimated_width = t[right_idx] - t[left_idx]
        
        
        return left_idx,right_idx,estimated_width


    def calculate_signal(self,pulse_data, left_idx, right_idx, ctr_idx):
        '''
        calculate SNR for pulse
        '''
        noise = pulse_data[self.offset_start:self.offset_end]  #np.concatenate((pulse_data[:left_idx], pulse_data[right_idx:]))
        #noise = np.squeeze(noise)
        signal = pulse_data[left_idx:right_idx]
        #signal_energy = np.sum(signal ** 2)
        #noise_energy = np.sum(noise ** 2)
        #snr = 10 * np.log10(signal_energy / noise_energy)
        signal_power = np.var(signal)
        noise_power = np.var(noise)
        # 计算信噪比
        snr = signal_power / noise_power
        
        
        noise_rms = np.sqrt(np.mean(np.square(noise)))
        peak_snr = np.max(pulse_data)/noise_rms
        
        return snr,peak_snr


    def set_SNR(self):
        '''
        set pulse profile SNR
        '''
        
        # 平均滤波
        length = len(self.pulse_profile)
        if length >=128:
            filtered_signal = self.average_filter(self.pulse_profile,64)
        else:
            filtered_signal = self.average_filter(self.pulse_profile,length//4)
        # 互相关 获得脉冲位置
        temp_poi = length//2
        temp_width = length//20
        temp_amp = 12
        template = self.pulse_template(temp_amp,temp_poi,temp_width,length)
        estimated_position = self.find_pulse_position(self.pulse_profile,template,temp_poi)
        est_left_idx,est_right_idx,est_width = self.find_pulse_range(filtered_signal,template,temp_poi)
        #print(self.pulse_profile.shape)
        #print("est pulse range:",est_left_idx,est_right_idx)
        #print("est pulse poisition:",estimated_position)
        
        
        # 计算信噪比
        self.snr,self.peak_snr = self.calculate_signal(self.pulse_profile,est_left_idx,est_right_idx,estimated_position)
        
        
        #plt.plot(self.pulse_profile,color="k")
        #plt.plot(template,color="yellow")
        #plt.plot(filtered_signal,color = "teal")
        #plt.vlines(estimated_position,min(self.pulse_profile),max(self.pulse_profile),color="r")
        #plt.vlines(est_left_idx,min(self.pulse_profile),max(self.pulse_profile),color="orange")
        #plt.vlines(est_right_idx,min(self.pulse_profile),max(self.pulse_profile),color="tomato")
        #plt.xlabel('Time')
        #plt.ylabel('Amplitude')
        #plt.title('Simulated Pulse Data SNR %.2f'%self.snr)
        #plt.show()
    
            
    def buttonClicked(self):
        print('Button clicked!')


def Usage():

    print("\n")
    print("-"*10)
    print("\n")
    print("Usage: python view_pulse.py filename")
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
        #Usage()
        print("\n")




if __name__ == '__main__':

    arg = sys.argv
    if "-h" in arg:
        Usage()
        sys.exit(1)
    
    file_list = get_filename(arg)
    for filename in file_list:
        app = QApplication(sys.argv)
        window = ViewPulse(filename)
        window.show()
        sys.exit(app.exec_())
