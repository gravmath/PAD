import numpy as np
from copy import deepcopy 
###########################################################################
#    equation_of_time.py has a set of utility functions that perform basic 
#       time-related tasks
#       
#       
#       Most formulae and concepts are taken from Chapters 7, 22, & 25 of
#       'Astronomical Algorithms' by Jean Meeus - hereafter refered to as
#       Meeus or [1]
#       
#       The calculation of GMST is taken from 'Approximate Sidereal Time' '
#       from http://aa.usno.navy.mil/faq/docs/GAST.php - hereafter refered 
#       to as [2]
#          
#       The time standard is Gregorian calendar format with UTC time 
#       system
#
#        deg = pi/180
############################################################################
def DecimalDay(epoch_dict):

    day  = epoch_dict['day']
    hour = epoch_dict['hour']
    min  = epoch_dict['min']
    sec  = epoch_dict['sec']

    decimal_day = day + (hour + (min + sec/60.0)/60.0)/24.0

    return decimal_day

############################################################################
def JulianDates(epoch_dict):

    #unpack the epoch_dict for convenience (notation follows Eq. 4 of [1]
    Y = epoch_dict['year']
    M = epoch_dict['month']
    D = DecimalDay(epoch_dict)

    #make the appropriate adjustment to Y and M (year and month) based on 
    #the discussion on page 3 of [1]
    if( M == 1 or M == 2):
        M = M + 12
        Y = Y - 1
    A = np.floor(Y/100.0)
    B = 2.0 - A + np.floor(A/4.0)

    #calculate the Julian day (jul_day) using Eq. 4 of [1]
    jul_day =   np.floor( 365.25  * (Y+4716.0) ) + \
                np.floor( 30.6001 * (M+1.0)    ) + \
                D + B - 1524.5
          
    #calculate the elapsed Julian centuries since Jan 01 2000, 12:00:00 UTC
    T  = ( jul_day - 2451545.0 ) / 36525
    
    return jul_day, T

    
############################################################################
def calculate_obliquity(T):

    #Calculate the obliquity
    obliquity = 23.43929111111111         + \
                0.0130041666666667*T      + \
                1.638888888888889e-7*T**2 + \
                5.036111111111111e-7*T**3 
                
    return obliquity

    
############################################################################    
def calculate_GMST(epoch_dict):

    #Get current Julian date and centuries since J2000.0 Epoch
    JD,T = JulianDates(epoch_dict)

    #Calculate elapsed days since J2000.0 epoch   
    D    = JD - 2451545.0
    
    #Calculate the decimal hours into the current day
    frac = DecimalDay(epoch_dict) % 1.0
    H    = frac*24.0
    
    #Calculate elapsed days of previous midnight since J2000.0 epoch
    D0   = D - frac
    
    GMST = 6.697374558 + 0.06570982441908*D0 + 1.000273790935*H + 0.000026*T**2
    
    return (GMST % 24.0)*2*np.pi/24.0
