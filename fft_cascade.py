#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      FINOT_M
#
# Created:     28/06/2015
# Copyright:   (c) FINOT_M 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

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
    date_end = line[:-1]
    f.close()
    return data_array, date_start, date_end

def fft_plot(period,Fs):
    n = period.shape[-1]
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(1,n/2)] # one side frequency range
    Y = np.fft.fft(period)/n # fft computing and normalization
    Y = Y[range(1,n/2)]

    plt.plot(frq,abs(Y),'r') # plotting the spectrum
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')
    return

def summary_digital(fullname,Fs = 400.0, timewindow = 2.0):
    data_array, date_start, date_end = load_data(fullname)
    size = data_array.shape[-1]
    width = int(Fs*timewindow)
    nt = range(1,width/2)
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
    nt = range(1,width/2)
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


def countour_fft(filename, Fs = 4000.0, timewindow = 10.0 ):
    ### create an image with a given width in time window in second.
    ### will great an fft image
    data_array, date_start, date_end = load_data(filename)
    img = []
    size = data_array.shape[-1]
    Fs = 4000.0
    width = Fs * timewindow  # number for fft window
    nline = int(round(size / width))  # number of windows
    nt = range(1,int(width/2))
    for i in range(0,nline):
        period = data_array[i*width:(i+1)*width]
        k = np.arange(width)
        T = timewindow
        frq = k/T # two sides frequency range
        frq = frq[nt] # one side frequency range
        Y =  np.fft.fft(period)/width # fft computing and normalization
        Y = abs(Y[nt])
        Y = attenuation(frq,Y,[60.0,120.0,240.0,360],width_hertz = 2.0)
        img.append(Y)
    return img

##period = data_array[0:1000000]
##fft_plot(period)
##time.sleep(100)
##period = data_array[1000000:2000000]
##fft_plot(period)
##
####data = pd.read_table(fullname,dtype = np.int8, header = 0, nrows = 7199750)
####data2 = pd.read_table(fullname2,dtype = np.int8, header = 0, nrows = 719975)
####ts = pd.Series(data, index=pd.date_range('1/1/2000', periods=7199750, freq='250U')
##
##sp = np.fft.fft(period)

def process_directory(directory):
    ### return a dataframe with information related to the files
    listfile = os.listdir(directory)
    location  = []
    letter = []
    date = []
    time = []
    typef = []
    tstamp = []
    filenamelist = []
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
                typef.append(row[4][:-4])
                tstamp.append(timestamp)
                filenamelist.append(filename)
            except:
                print filename

    table_file = pd.DataFrame()
    table_file['location']=location
    table_file['letter']=letter
    table_file['date']=date
    table_file['time']=time
    table_file['typef']=typef
    table_file['filename']=filenamelist
    table_file['timestamp']=tstamp

    table_sort = table_file.sort(['timestamp','letter','typef'])
    return table_sort

def process_analog(table_sort, directory, target_directory =r'C:\Sandbox\coursera\cascade\data'):
    ### exedcute a series of operation on a dataframe of files defined by process_directory
    table_analog = table_sort[(table_sort.typef=='analog') & (table_sort.letter =='B') & (table_sort.timestamp < '2015-01-05') &(table_sort.timestamp > '2015-01-02')]
    Fs = 4000.0;  # sampling rate
    Ts = 1.0/Fs; # sampling interval
    t = np.arange(0,1,Ts) # time vector

    for filename in table_analog.filename:
        file_digital = filename[:-10] + 'digital.txt'
        fullname = os.path.join(directory,filename)
        full_digital = os.path.join(directory,file_digital)
        digital_info = summary_digital(fullname,Fs = 4000.0, timewindow = 2.0)
        #image_series(fullname,Fs = 4000.0, timewindow = 10.0 )
        img = countour_fft(fullname,Fs = 4000.0, timewindow = 2.0)
        plt.cla
        plt.axis([0,1000,0,1000])
        plt.imshow(img)
        fname = filename[:-4]+'.png'
        output = os.path.join(target_directory,fname)
        print('Saving frame', output)
        plt.savefig(output)
    return

def process_digital(table_sort, directory, target_directory =r'C:\Sandbox\coursera\cascade\data'):
    ### exedcute a series of operation on a dataframe of files defined by process_directory
    table_digital = table_sort[(table_sort.typef=='digital') & (table_sort.letter =='B') & (table_sort.timestamp < '2015-01-02') &(table_sort.timestamp > '2015-01-01')]
    img = []
    for filename in table_digital.filename:
        fullname = os.path.join(directory,filename)
        print filename
        summary = summary_digital(fullname,Fs = 400.0, timewindow = 2.0)
        #image_series(fullname,Fs = 4000.0, timewindow = 10.0 )
        img.append(summary)

    plt.cla
#    plt.axis([0,1000,0,1000])
    plt.imshow(img)
#    plt.show()
    fname = 'summary_test.png'
    output = os.path.join(target_directory,fname)
    print('Saving frame', output)
    plt.savefig(output)
    return img

def main():
    fullname = r'D:\Cascade Experiments\CascadeBath_C_Jan01_0155_analog.txt'
    fullname = r'D:\Cascade Experiments\CascadeBath_C_Jan01_0155_digital.txt'
    #image_series(fullname,Fs = 4000.0, timewindow = 40.0, plot_range = [0,1000,0,1] )
    directory = r'D:\Cascade Experiments'

    table_sort = process_directory(directory)
    img = process_digital(table_sort, directory, target_directory =r'C:\Sandbox\coursera\cascade\data')

    #process_analog(table_sort,directory, target_directory =r'C:\Sandbox\coursera\cascade\data')
    #summary = summary_digital(fullname)
    #plt.plot(summary)
    pass

if __name__ == '__main__':
    main()

##timestamp = dt.datetime.strptime(dtt,'%b%d %H%M')