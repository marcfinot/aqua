######################################################
# Reformats the flow file from
# timing and the original data stream. Puts it in a similar.
# format to the transducer file.
# Revision 1 writes in CSV, and truncates when required. It assumes there
# is a flow switch on Bit 0 and a flow meter on Bit 1. The flowmeter
# calibration is assumed to be 'ticks per liter', and is between either
# rising-to-rising or falling-to-falling edge.
#
#  Modifivcation  conversion of ReformatdigFile1.wak into python
#Marc Finot
# epoch time in second
# switch position: on: 1, off: 0
# 2015-8-13
# compute directly the interval, start and stop tie from the read file.
#
#  output: file with name filename_flow.txt
######################################################
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt


def load_data(fullname):
    ''' import a log file with date on first and last line
    parameters:
    fullname  - full path filename of the file to import
    output:
    data_array, date_start, date_end
    '''
    f = open(fullname,'r')
    data = []
    date_start = pd.to_datetime(f.readline()[:-1])
    for line in f:
        if len(line) <5:
            data.append(np.int(line))
    data_array = np.array(data)
    date_end = pd.to_datetime(line)
    f.close()
    return data_array, date_start, date_end

def reformatdigfilex(filename,Offset_Strt=0,Offset_Stp=0, save = True):
    # Time_inc = ARGV[2]
    # Start_Time = ARGV[3]
    # Stop_Time = ARGV[4]
    # Offset_Strt = ARGV[5]
    # Offset_Stp = ARGV[6]

    outputfile = filename[:-4] + '_flow.txt'
    TicksPerLiter = 200.
    SwitchPolarity = "POSITIVE"
    InstFlowRate = 0.0
    FlowMeterCount = 1.
    LastFlowMeterState = 0
    FlowRate = 0
    FlowRateThreshold = 0.01  #Threshold in liters per minute, below which flow is truncated to zero.

    AvgLength = 100
    CBHeadPtr = 1
    flow = []
    switch = []
    time_stamp = []
    #######################################
    #Set up and zero out a circular buffer,
    #used for calculating averages.
    #######################################

    CircBuff= np.zeros(AvgLength)

    if (Offset_Strt == 0) & (Offset_Stp == 0):
        OffsetAction = "NONE"
    else:
        OffsetAction = "YES"

    #loading the file and reading first and last time
    table, start_date, stop_date = load_data(filename)

    Start_Time = time.mktime(start_date.timetuple())
    Stop_Time=time.mktime(stop_date.timetuple())

    diff = stop_date-start_date
    time_diff_second = diff.seconds
    len_data= len(table)
    # compute time step based on number of points

    Time_inc = time_diff_second/np.float(len_data)
    Timestamp = Time_inc

    if save: g = open(outputfile,'w')
    for data in table:                   # lose the header
        ##########################################################################
        # Determine whether the flow switch is on or off, and what the flow rate
        # is. The flow rate is figured out by counting the amount of time
        # between successive ticks of the flow meter. Anything below the
        # <flow threshold> value is considered zero.
        ##########################################################################
        Bit0 = np.fmod(data,2)   #strip out bit 0 (switch).
        Bit1 = (data - Bit0)/2  #strip out bit 1 (flowmeter)
        FlowSwitch = 0
        if (SwitchPolarity == "POSITIVE") & ( Bit0 == 1 ):
            FlowSwitch = 1
        if (SwitchPolarity == "NEGATIVE") & ( Bit0 == 0 ):
            FlowSwitch = 1

        if(LastFlowMeterState == 0 ) & ( Bit1 == 1 ):            ##Have a rising edge on the flowmeter
            InstFlowRate = 60./(FlowMeterCount * Time_inc * TicksPerLiter)    ## Instantaneous Flow rate in liters per minute
            if( InstFlowRate < FlowRateThreshold):
                InstFlowRate = 0
            FlowMeterCount = 0

        FlowMeterCount=FlowMeterCount+1
        LastFlowMeterState = Bit1


        CircBuff[CBHeadPtr] = InstFlowRate
        CBHeadPtr=CBHeadPtr+1
        if(CBHeadPtr == AvgLength):
            CBHeadPtr = 0                ##Roll circular buffer over

        FlowRate = 0
        for i in range(0,AvgLength):
            FlowRate = FlowRate + CircBuff[i]

        FlowRate = FlowRate / AvgLength

        if FlowSwitch == 0:
            FlowRate = 0                    ##FLow meter has momentium. Flowswitch is off, no flow though.

        if OffsetAction == "NONE" :
            strf = str(Timestamp + Start_Time)+','+ str(FlowSwitch)+','+str(FlowRate)+'\n'
        #printf("%s, %d, %d, %d\n", OffsetAction, (Timestamp + Start_Time), Offset_Strt,Offset_Stp)
        if (OffsetAction == "YES") & ((Timestamp + Start_Time) >= Offset_Strt) & ((Timestamp + Start_Time) < Offset_Stp):
            strf = str(Timestamp + Start_Time)+','+ str(FlowSwitch)+','+str(FlowRate)+'\n'

        if save: g.write(strf)

        Timestamp = Timestamp + Time_inc
        flow.append(FlowRate)
        switch.append(FlowSwitch)
        time_stamp.append(Timestamp+Start_Time)
    if save: g.close()
    return time_stamp, switch, flow

def main():
    pass

if __name__ == '__main__':
    filename = r'C:\Users\home\Documents\Bob\Cascade Experiments\CascadeKitchen_A_Nov01_1829_digital.txt'
#    time_stamp, switch, flow = reformatdigfilex(filename, save= False)
#    plt.plot(time_stamp,flow)
#    plt.plot(time_stamp,switch)
#    plt.show()
    main()
