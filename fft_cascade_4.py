#-------------------------------------------------------------------------------
# Name:        fft_cascade with original files.
# Purpose:
#
# Author:      FINOT_M
#
# Created:     28/06/2015
# Copyright:   (c) FINOT_M 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

''' 
idea:
- random sampling for SVD over time and year due to 
'''

# to run matplotlib without window
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import cPickle
from scipy import signal
import time
import ReformatDigFile1 as reformat

## example for wavelet analysis
def wavelet():
    t = np.linspace(-1, 1, 200, endpoint=False)
    sig  = np.cos(2 * np.pi * 7 * t) + signal.gausspulse(t - 0.4, fc=2)
    widths = np.arange(1, 31)
    cwtmatr = signal.cwt(sig, signal.ricker, widths)
    plt.imshow(cwtmatr, extent=[-1, 1, 1, 31], cmap='PRGn', aspect='auto', vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())
    plt.show()
    return


def svd_reduction(mat,k):
    mat_norm =  mat
    # svd decomposition   a =  u * np.diag(s) * v
    U,s,V = np.linalg.svd(mat_norm,full_matrices=False)
    # to reconstruct
    S = np.diag(s)
    #np.allclose(a, np.dot(U, np.dot(S, V)))
    R = np.dot(np.diag(s),V)
    z = np.dot(U.T,mat_norm)
    b = 0.0
    sumS = sum(s)
    s_sum = []
    for i in s:
        b=b+i/sumS
        s_sum.append(b)
    U_reduced = U[:,0:k]
    S_reduced = S[0:k]
    z_reduced = np.dot(U_reduced.T,mat_norm)
    mat_reduced = np.dot(U_reduced,np.dot(S_reduced, V))
    return mat_reduced, U_reduced, S_reduced, s_sum


def SVD_spectrum(a):
    ### compute the Signular value decompostion for a set of spectrum  a ~  u * np.diag(s) * v
    U, s, V = np.linalg.svd(a, full_matrices=False)
    return U, s, V


def attenuation(X,Y,freq_range, width_hertz = 2.0):
    ### remove peak frequencies ina frequency range using a gaussien filter
    for freq in freq_range:
        att = 1.0 - np.exp(-(X-freq)*(X-freq)/width_hertz)
        Y = Y  * att
    return Y

def gaussien_filter(X,Y,freq_range, width_hertz = 2.0):
    ### remove peak frequencies ina frequency range using a gaussien filter
    for freq in freq_range:
        att = np.exp(-(X-freq)*(X-freq)/width_hertz)
        Y = Y  * att
    return Y

def load_data(fullname):
    #datalength = 7199750
    f = open(fullname,'r')
    data = []
    date_start = f.readline()[:-1]
    for line in f:
        if len(line) <5:
            data.append(np.int(line))
    data_array = np.array(data)
    date_end = line
    f.close()
    return data_array, date_start, date_end

def fft_plot(period,Fs):
    n = period.shape[-1]
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(0,n/2)] # one side frequency range
    Y = np.fft.fft(period)/n # fft computing and normalization
    Y = Y[range(0,n/2)]

    plt.plot(frq,abs(Y),'r') # plotting the spectrum
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')
    return

def summary_digital(fullname,Fs = 400.0, timewindow = 2.0):
    data_array, date_start, date_end = load_data(fullname)
    size = data_array.shape[-1]
    width = int(Fs*timewindow)
    nt = range(0,width)
    nline = int(round(size / width))  # number of windows
    summary =[]
    for i in range(0,nline):
        period = data_array[i*width:(i+1)*width]
        summary.append(np.sum(period)/np.float(width))
    return summary

def image_series(filename,Fs = 4000.0, timewindow = 10.0 , plot_range = [0,1000,0,1], nb_grap = 10):
    data_array, date_start, date_end = load_data(filename)
    width = int(Fs*timewindow)
    size = data_array.shape[-1]
    nt = range(0,width/2)
    nline = int(round(size / width))  # number of windows
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')
    plt.axis(plot_range)

    for i in range(0,min(nline, nb_graph)):
        print i
        period = data_array[i*width:(i+1)*width]
        n = period.shape[-1]
        k = np.arange(n)
        T = n/Fs
        frq = k/T # two sides frequency range
        frq = frq[nt] # one side frequency range
        Y =  np.fft.fft(period)/n # fft computing and normalization
        Y = Y[nt]
        Y = attenuation(frq,Y,[60.0,120.0,240.0],width_hertz = 2.0)
        plt.plot(frq,abs(Y),'r') # plotting the spectrum
        fname = filename[:-4]+'%03d.png'%i
        print('Saving frame', fname)
        plt.savefig(fname)
    return


def countour_fft(array, Fs = 4000.0, timewindow = 10.0,fileformat = False ):
    ''' create an image with a given width in time window in second.
     will great an fft image
     @param: array  linear array of all the spectrum
     @param: Fs  frequency array
     @return:   img, frq, nline: set of spectrum, frequency array, number of spectrum
    '''
    if fileformat:
        filename  = array
        data_array, date_start, date_end = load_data(filename)
    else:
        data_array = array
    img = []
    size = data_array.shape[-1]
    width = Fs * timewindow  # number for fft window
    nline = int(round(size / width))  # number of windows
    nt = range(0,int(width/2))
    for i in range(0,nline):
        period = data_array[i*width:(i+1)*width]
        k = np.arange(width)
        T = timewindow
        frq = k/T # two sides frequency range
        frq = frq[nt] # one side frequency range
        Y =  np.fft.fft(period)/width # fft computing and normalization
        Y = abs(Y[nt])
        # poossible filtering of 60Hz frequency  harmonics
        #Y = attenuation(frq,Y,[60.0,120.0,240.0,360],width_hertz = 2.0)
        img = np.concatenate((img,Y))
    return img, frq, nline


def process_directory(directory):
    ''' return a dataframe with information related to the files
    @param directory name
    @return dataframe with information about all the files in the directory
    '''
    listfile = os.listdir(directory)
    location  = []
    letter = []
    date = []
    time = []
    typef = []
    tstamp = []
    filenamelist = []
    fullpath = []
    for filename in listfile:
        if filename[-3:]=='txt':
            name = filename.replace(' ','_')
            if 'Main_F' not in filename:
                name = name.replace('Main','Main_F')
            row = name.split('_')
            #print name
            try:
                location.append(row[0])
                letter.append(row[1])
                date.append(row[2])
                time.append(row[3])
                dtt = '2015 '+row[2]+' '+row[3]
                timestamp = dt.datetime.strptime(dtt,'%Y %b%d %H%M')
                if timestamp > dt.datetime(2015,4,1):
                    dtt = '2014 '+row[2]+' '+row[3]
                    timestamp = dt.datetime.strptime(dtt,'%Y %b%d %H%M')
                typef.append(row[4][:-4])
                tstamp.append(timestamp)
                filenamelist.append(filename)
                fullpath.append(os.path.join(directory,filename))
            except:
                print filename

    table_file = pd.DataFrame(location)
    table_file['location']=location
    table_file['letter']=letter
    table_file['date']=date
    table_file['time']=time
    table_file['typef']=typef
    table_file['filename']=filenamelist
    table_file['timestamp']=tstamp
    table_file['fullpath']=fullpath

    table_sort = table_file.sort(['timestamp','letter','typef'])
    return table_sort

def process_analog_0(table_sort, directory, target_directory =r'C:\Sandbox\coursera\cascade\data',letter='A', time_start = '2015-01-01',time_stop = '2015-01-15', plotting = True,Fs = 2000.0, timewindow = 4.0):
    ''' exedcute a series of operation on a dataframe of files defined by process_directory
    @param: table_sort  life of file 
    @param: directory 
    @param: target_directory to save data    
    
    '''

    table_analog = table_sort[(table_sort.typef=='analog') & (table_sort.letter ==letter) & (table_sort.timestamp < time_stop) &(table_sort.timestamp > time_start)]
    Ts = 1.0/Fs; # sampling interval
    t = np.arange(0,1,Ts) # time vector
    full_set = []
    for filename in table_analog.filename:
        #file_digital = filename[:-10] + 'digital.txt'
        fullname = os.path.join(directory,filename)
        #full_digital = os.path.join(directory,file_digital)
        #digital_info = summary_digital(fullname,Fs = 400.0, timewindow = 2.0)
        #image_series(fullname,Fs = 4000.0, timewindow = 10.0 )
        img, frq, nline = countour_fft(fullname,Fs = Fs, timewindow = timewindow)
        full_set = np.concatenate((full_set,img))
        width = Fs*timewindow
        b = round(len(img)/width)
        img2 = img[0:(b*width)].reshape((b,width))
        if plotting:
            plt.clf()
            plt.axis([0,1000,0,1000])
            plt.imshow(img2)
            fname = filename[:-4]+'.png'
            output = os.path.join(target_directory,fname)
            print('Saving frame', output)
            plt.savefig(output)
    return full_set

def process_analog(filelist):
    ''' execute a series of operation on a dataframe of files defined by process_directory
    Parameters:
    ----------
    filelist - dlist of files to process 
    figoutput - name of file for saving image
    return   combined_df: dataframe with sensor signal
    
    '''
    first = True

    for fullname in filelist:
        print 'processing:', fullname
        data_analog, start_date, stop_date = reformat.load_data(fullname)
        Start_Time = time.mktime(start_date.timetuple())
        Stop_Time=time.mktime(stop_date.timetuple())
    
        diff = stop_date-start_date
        time_diff_second = diff.seconds
        len_data= len(data_analog)
        analog_period = time_diff_second/np.float(len_data)
    
        # analog dataframe
        analog_df = pd.DataFrame(data_analog)
        analog_df.columns = ['sensor_amp']
        ts = np.arange(Start_Time,Stop_Time,analog_period)
        analog_df['timestamp'] = ts[:len_data]
        #normalization
        analog_df['sensor_amp'] = analog_df['sensor_amp'] - np.mean(analog_df['sensor_amp'])
        array_length1 = len(analog_df['timestamp'])
 
        Fs = int(np.float(len_data)/time_diff_second) # can be calculated based on time stamp.
    
        rng1 = pd.date_range(analog_df['timestamp'][0], periods = array_length1,freq = str(int(1000000/Fs))+'u')
        analog_df.index = rng1
        if first:
            combined_df = analog_df
            first = False
        else:
            combined_df =pd.concat([combined_df, analog_df]) 
        
    return combined_df, Fs
            
    
def reduce_analog(analog_df,Fs=2000.0, timewindow = 4.0, dimension = 20):    
    ''' reduce an analog signal to a set  of spectrum
    analog_df : dataframe with timestamp and signal amplitude
    '''    
    array_length1 = len(analog_df['timestamp'])
    ts1 = np.arange(0,array_length1)/(Fs*timewindow)
    # new time range
    width = Fs*timewindow
    img, frq, nline = countour_fft(analog_df['sensor_amp'].values,Fs = Fs, timewindow = timewindow)
    # refrmating image in 2D array
    rng = pd.date_range(analog_df['timestamp'][0], periods = nline,freq = str(int(timewindow))+'s')
    b = int(img.shape[0]/nline)
    img2 = img[0:b*nline].reshape((nline,b))
    mat_reduced, U_reduced, S_reduced, s_sum = svd_reduction(img2.T,dimension)
    z_reduced = np.dot(U_reduced.T,img2.T)
    
    return rng, z_reduced, img2,nline, frq, U_reduced
    
def process_digital(filelist):
    ''' execute a series of operation on a set digital files and conver them into a dataframe
    Parameters:
    ----------
    filelist - dlist of files to process 
    figoutput - name of file for saving image
    return   combined_df: dataframe with timestamp, flow and state of valve
    
    '''
    first = True
    for fullname in filelist:
        print 'processing:',fullname
        try:
            time_stamp, switch, flow = reformat.reformatdigfilex(fullname, save= False)
            digit_time_span = time_stamp[-1]-time_stamp[0]
    
        #digital datafrme
            digital_df = pd.DataFrame(time_stamp,columns = ['timestamp'])
            digital_df['state']= switch
            digital_df['flow'] = flow
            array_length2 = len(digital_df['timestamp'])
            
            Fs2 = int(array_length2/digit_time_span)   # can be calculated based on time stamp.
            rng2 = pd.date_range(digital_df['timestamp'][0], periods = array_length2,freq = str(int(1000000/Fs2))+'u')
            digital_df.index = rng2
            
            if first:
                combined_df = digital_df
                first = False
            else:
                combined_df =pd.concat([combined_df, digital_df]) 
        except:
            print 'issue with:',fullname
        
    return combined_df, Fs2

def compute_monthly_SVD():
    ''' compute the SVD of one day per month (first day of the month)
    save it as a cpickle file
    '''
    directory = r'C:\Users\home\Documents\Bob\Cascade Experiments'
    table_sort = process_directory(directory)
    target_directory =r'C:\Users\home\Documents\Bob\cascade\data'
    time_start_list = ['2014-08-01','2014-09-01','2014-08-01','2014-10-01','2014-11-01','2014-12-01']
    time_stop_list = ['2014-08-02','2014-09-02','2014-08-02','2014-10-02','2014-11-02','2014-12-02']
    for i in range(len(time_start_list)):
        time_start = time_start_list[i]
        time_stop = time_stop_list[i]
        for letter in ['A','B','C','D','E','F']:

            #img = process_digital(table_sort, directory, target_directory =r'C:\Users\home\Documents\Bob\cascade\data',letter = letter, time_start = '2015-02-01',time_stop = '2015-02-1' )
            full_set = process_analog(table_sort, directory, target_directory =r'C:\Users\home\Documents\Bob\cascade\data',letter = letter, time_start = time_start,time_stop = time_stop,Fs = 2000.0, timewindow = 4.0  )
            width = Fs*timewindow
            b = len(full_set)/width
            if b!=0:
                fft_array = full_set[0:(b*width)].reshape((b,width))
                U,s,V = np.linalg.svd(fft_array, full_matrices=False)
                plt.clf()
                fV = 'V_spectrum_'+letter +'_'+ time_start +'_'+time_stop+'.pickle'
                f = open(fV,'w')
                cPickle.dump(V,f)
                cPickle.dump(s,f)
                cPickle.dump(letter,f)
                cPickle.dump(time_start,f)
                cPickle.dump(time_stop,f)
                f.close()
                for i in range(0,6):
                    plt.plot(V[i])
                fname = 'analog_spectrum_'+letter +'_'+ time_start +'_'+time_stop+'.png'
                output = os.path.join(target_directory,fname)
                plt.savefig(output)
                plt.clf()
                plt.plot(s[0:20])
                fname = 'analog_s_coeff_'+letter +'_'+ time_start +'_'+time_stop+'.png'
                output = os.path.join(target_directory,fname)
                plt.savefig(output)
    return

def plot_profile():
    plt.subplot(2, 1, 1)
    plt.plot(rng, U_reduced[:,0])
    plt.xlabel('time')
    plt.subplot(2, 1, 2)
    plt.plot(mat_reduced[0])
    plt.xlabel('Frequency')
    return



def correlation_analog_digital(filelist, title = 'test'):
    ''' compute the correlation between the analog and the digital signal
    @param name: name of the file to analyze
    
    @todo different frequency for analog_periodg and digital - need to convert to same timescale
    '''
    
    #name = r'C:\Users\home\Documents\Bob\Cascade Experiments\CascadeKitchen_A_Nov01_1829'
    digital_df, Fs_digital =  process_digital(filelist)    
    analog_df, Fs =  process_analog(filelist) 
   
    reduced_analog_df, img2, nline, frq, U_reduced = reduce_analog(analog_df,Fs=2000.0, timewindow = 4.0, dimension = 20)
    reduced_digital_df = digital_df.resample(str(int(timewindow))+'s', how='mean')
    
    plt.close()
    plt.subplot(4,1,2)
    plt.xlabel('time')
    plt.ylabel('sensor')
#    plt.plot(dataframe1['timestamp'],dataframe1['sensor_amp'])
    plt.plot(analog_df['sensor_amp'])
    plt.subplot(4,1,3)
    plt.xlabel('time')
    plt.ylabel('flow')
    plt.plot(digital_df['flow'])

#    plt.subplot(3,1,1)
    
    z_reduced = np.dot(U_reduced.T,img2.T)
    
    for i in [0,1,2,3]:
        plt.subplot(4, 1, 1)
        plt.plot(frq, U_reduced[:,i])
        plt.ylabel('ampl')
        plt.subplot(4, 1, 4)
        plt.ylabel('spec. ampl.')
        plt.plot(z_reduced[i])
        plt.xlabel('time')
    output = r'C:\Users\home\Documents\Bob\steve\correlation_'+title+'.png'
    plt.savefig(output)
    #plt.show()
    plt.close()
    # resampling with the same period as the fft 

    
    x = reduced_digital_df['flow'].values[:nline]
    plt.ylabel('spect. amplitude')
    plt.xlabel('flow amplitude')
    xrg = np.arange(np.min(x),np.max(x),(np.max(x)-np.min(x))/30.0)
    for i in [0,1,2,3]:
        fit = np.polyfit(x,z_reduced[i],1)
        fit_fn = np.poly1d(fit)
        plt.plot(x,z_reduced[i],'.', xrg,fit_fn(xrg),'-')
    #plt.ylim(-0.03, 0.01)
    output = r'C:\Users\home\Documents\Bob\steve\correlation_flow_spectrum_quad_'+title+'.png'
    plt.savefig(output)
    return

def main():
    fullname = r'C:\Users\home\Documents\Bob\Cascade Experiments\CascadeBath_C_Jan01_0155_analog.txt'
    #image_series(fullname,Fs = 4000.0, timewindow = 40.0, plot_range = [0,1000,0,1] )
    directory = r'C:\Users\home\Documents\Bob\Cascade Experiments'

    table_sort = process_directory(directory)
    img, combined_df = process_digital(table_sort, directory, target_directory =r'C:\Users\home\Documents\Bob\cascade\data',Fs = 200.0, timewindow = 4.0 )

    #process_analog(table_sort,directory, target_directory =r'C:\Sandbox\coursera\cascade\data')
    #summary = summary_digital(fullname)
    #plt.plot(summary)
    pass

if __name__ == '__main__':
    date = ['Nov01_1929','Nov01_2029','Nov01_2129','Nov02_0729','Nov02_0829','Nov02_1829','Nov02_1929']
    name = r'C:\Users\home\Documents\Bob\Cascade Experiments\CascadeKitchen_A_Nov01_1829'
    for dd in date:
        name = r'C:\Users\home\Documents\Bob\Cascade Experiments\CascadeKitchen_A_'+dd
        #correlation_analog_digital(name)
        
    directory = r'C:\Users\home\Documents\Bob\Cascade Experiments'
    table_sort = process_directory(directory)
    target_directory =r'C:\Users\home\Documents\Bob\cascade\data'
    letter = 'A'
    time_start = '2014-11-01 8:00'
    time_stop = '2014-11-01  12:00'
    Fs = 2000.0
    timewindow = 4.0
#    full_set = process_analog(table_sort, directory, target_directory =r'C:\Users\home\Documents\Bob\cascade\data',letter = letter, time_start = time_start,time_stop = time_stop,Fs = Fs, timewindow = 4.0 )
    
    table_digital = table_sort[(table_sort.typef=='digital') & (table_sort.letter == letter) & (table_sort.timestamp < time_stop) & (table_sort.timestamp > time_start)]
    table_analog = table_sort[(table_sort.typef=='analog') & (table_sort.letter == letter) & (table_sort.timestamp < time_stop) & (table_sort.timestamp > time_start)]
    fname = 'digital_'+letter +'_'+ time_start +'_'+time_stop+'_'+str(int(timewindow))+'s'+'.png'
    figoutput = os.path.join(target_directory,fname)

    filelist_digital = table_digital['fullpath']
    filelist_analog = table_analog['fullpath']
    analog_df, Fs_analog =  process_analog(filelist_analog)
    rng, z_reduced, img2,nline, frq, U_reduced = reduce_analog(analog_df,Fs=2000.0, timewindow = 4.0, dimension = 20)
    #digital_df, Fs_digital = process_digital(filelist_digital)
    #digital_reduced  = digital_df.resample(str(int(timewindow))+'s', how='mean')

##timestamp = dt.datetime.strptime(dtt,'%b%d %H%M')