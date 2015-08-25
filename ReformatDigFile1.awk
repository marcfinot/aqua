######################################################
# Reformats the flow file from
# timing and the original data stream. Puts it in a similar. 
# format to the transducer file. The time is absolute,
# from teh start of 2014 in sec.
# Revision 1 writes in CSV, and truncates when required. It assumes there
# is a flow switch on Bit 0 and a flow meter on Bit 1. The flowmeter 
# calibration is assumed to be 'ticks per liter', and is between either 
# rising-to-rising or falling-to-falling edge.
#
#
# Calling: mawk -f reformatdigfilex <file> <time_inc> <Start_time> <Stop_Time> <offset_start> <offset_stop>
# <time_inc>: time increment between records, in sec
# <Start_time>: start time in secs since beginning 2014
# <Stop_time>: stop time in secs since beginning 2014
# <offset_start>: start time of subsection to convert, secs since beginning of 2014
# <offset_stop>: stop time of subsection to convert, secs since beginning of 2014
#
#
######################################################





BEGIN	{
	Time_inc = ARGV[2]
	Start_Time = ARGV[3]
	Stop_Time = ARGV[4]
	Offset_Strt = ARGV[5]
	Offset_Stp = ARGV[6]

	TicksPerLiter = 200
	SwitchPolarity = "POSITIVE"
  
	FlowMeterCount = 1
	LastFlowMeterState = 0
	FlowRate = 0
	FlowRateThreshold = 0.01  #Threshold in liters per minute, below which flow is truncated to zero.

	AvgLength = 100
	CBHeadPtr = 1

	#######################################
	#Set up and zero out a circular buffer,
	#used for calculating averages.
	#######################################
	for(i = 1; i < AvgLength; i++ ) {
	  CircBuff[i] = 0
	  }

	Timestamp = Time_inc
	ARGC = 2
	
	if( (Offset_Strt == 0) && (Offset_Stp == 0) ) OffsetAction = "NONE"
	else OffsetAction = "YES"
	}

	{
	if(NR > 1 ) {					# lose the header
##########################################################################
# Determine whether the flow switch is on or off, and what the flow rate
# is. The flow rate is figured out by counting the amount of time 
# between successive ticks of the flow meter. Anything below the
# <flow threshold> value is considered zero.
##########################################################################

	  Bit0 = 2*(($1/2)%1)	#strip out bit 0 (switch).
	  Bit1 = ($1 - Bit0)/2  #strip out bit 1 (flowmeter)
	  FlowSwitch = "OFF"
	  if((SwitchPolarity == "POSITIVE") && ( Bit0 == 1 )) FlowSwitch = "ON"
	  if((SwitchPolarity == "NEGATIVE") && ( Bit0 == 0 )) FlowSwitch = "ON"

	  if( (LastFlowMeterState == 0 ) && ( Bit1 == 1 ) ) {			##Have a rising edge on the flowmeter
	    InstFlowRate = 60/(FlowMeterCount * Time_inc * TicksPerLiter)	## Instantaneous Flow rate in liters per minute
	    if( InstFlowRate < FlowRateThreshold) InstFlowRate = 0
	    FlowMeterCount = 0
	    }
	
	  FlowMeterCount++
	  LastFlowMeterState = Bit1

	  CircBuff[CBHeadPtr++] = InstFlowRate
	  if(CBHeadPtr == AvgLength + 1) CBHeadPtr = 1				##Roll circular buffer over 
	  
	  FlowRate = 0
	  for(i = 1; i <= AvgLength; i++) {
	    FlowRate = FlowRate + CircBuff[i]
	    }
	  FlowRate = FlowRate / AvgLength


	  if(FlowSwitch == "OFF") FlowRate = 0					##FLow meter has momentium. Flowswitch is off, no flow though.

	  if( OffsetAction == "NONE" ) {
	    printf("%f,%s, %f\n", Timestamp + Start_Time, FlowSwitch, FlowRate)
	    }
#printf("%s, %d, %d, %d\n", OffsetAction, (Timestamp + Start_Time), Offset_Strt,Offset_Stp)
	  if( (OffsetAction == "YES") && ((Timestamp + Start_Time) >= Offset_Strt) && ((Timestamp + Start_Time) < Offset_Stp) ) {
	    printf("%f,%s, %f\n", Timestamp + Start_Time, FlowSwitch, FlowRate)
	    }
	  Timestamp = Timestamp + Time_inc
	  }
	}