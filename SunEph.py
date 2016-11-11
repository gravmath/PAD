import numpy as np
import equation_of_time as equ_of_time
#import pdb
deg          = np.pi/180

def CalcSun_Low( epoch_dict ):
    '''CalcSunPos returns the sun's position with respect to the 
       earth for a specified epoch.  
       
       The formulae and concepts are taken from Chapters 7, 22, & 25 of
       'Astronomical Algorithms' by Jean Meeus - hereafter refered to as
       Meeus
   
       The time standard is Gregorian calendar format with UTC time standard
       deg = pi/180'''

    #1 - Get the Julian day (jul_day) and Julian century (T)
    (jul_day, T) = equ_of_time.JulianDates(epoch_dict)


    #2 - get sun distance, sun right ascension and declination
    (sun_dist, sun_RA, sun_dec, d_sun_dist, d_sun_RA, d_sun_dec) = CalcSunSpherical(T)

   
    #3 - construct the sun position wrt the center of the earth
    sun_pos = np.array([ np.cos(sun_RA*deg)*np.cos(sun_dec*deg), np.sin(sun_RA*deg)*np.cos(sun_dec*deg), np.sin(sun_dec*deg)])
    sun_pos = sun_dist*sun_pos
    
        
    sun_vel = [0,0,0]

    sun_vel = d_sun_dist*sun_pos/sun_dist + \
              sun_dist*deg*np.array([-d_sun_RA*np.sin(sun_RA*deg)*np.cos(sun_dec*deg) - d_sun_dec*np.cos(sun_RA*deg)*np.sin(sun_dec*deg), \
                             d_sun_RA*np.cos(sun_RA*deg)*np.cos(sun_dec*deg) - d_sun_dec*np.sin(sun_RA*deg)*np.sin(sun_dec*deg), \
                             d_sun_dec*np.cos(sun_dec*deg)])
                             
    return (sun_pos, sun_vel)
               
###########################################################################################
def CalcSunSpherical(T):
    
    #Calculate the mean longitude (L0) and the mean anomaly of the Sun (M0)
    L0 = 280.46646 + 36000.76983*T + 0.0003032*T**2
    M0 = 357.52911 + 35999.05029*T - 0.0001537*T**2

    #Calculate the 'center of the sun' 
    C0 = 1.914602 - 0.004817*T - 0.000014*T**2
    C1 = 0.019993 - 0.000101*T
    C2 = 0.000289
    C  = C0*np.sin(M0*deg) + C1*np.sin(2*M0*deg) + C2*np.sin(3*M0*deg)
    
    #Calculate the true longitude (true_long) and the true anomaly of the Sun (true_anom)
    true_long = L0 + C
    true_anom = M0 + C

    #Calculate the obliquity
    obliquity = equ_of_time.calculate_obliquity(T)

    #Calculate the sun's right ascension (sun_RA) and declination (sun_dec), both in
    #degrees, using equations 25.6 and 25.7 in 'Astronomical Algorithms' by Jean Meeus
    sun_RA  = np.arctan2(np.cos(obliquity*deg)*np.sin(true_long*deg),np.cos(true_long*deg))/deg
    sun_dec = np.arcsin(np.sin(obliquity*deg)*np.sin(true_long*deg))/deg

    #Calculate the eccentricity and radius of the Sun in its orbit
    e        = 0.016708634 - 0.000042037*T - 0.0000001267*T**2
    R        = (1.000001018*(1-e**2))/(1+e*np.cos(true_anom*deg))

    #Finally rescale distance to kilometers
    AU       = 149597871.0
    sun_dist = R*AU

    #Calculate the derivatives of the mean longitude (L0) and the mean anomaly of the Sun (M0)
    dL0_dT = 36000.76983 + 2*0.0003032*T
    dM0_dT = 35999.05029 - 2*0.0001537*T

    #Calculate the derivatives with the 'center of the sun' 
    dC0_dT = -0.004817 - 2*0.000014*T
    dC1_dT = -0.000101
    dC_dT  = dC0_dT*np.sin(M0*deg)          + dC1_dT*np.sin(2.0*M0*deg) + \
             dM0_dT*deg*( C0*np.cos(M0*deg) + 2.0*C1*np.cos(2.0*M0*deg) + 3.0*C2*np.cos(3.0*M0*deg) )

    #Calculate the derivatives of the true longitude (true_long) and the true anomaly of the Sun (true_anom)
    dtrue_long_dT = dL0_dT + dC_dT
    dtrue_anom_dT = dM0_dT + dC_dT

    #Calculate the derivatives of the obliquity
    dobliquity_dT =   0.0130041666666667      +\
                    2*1.638888888888889e-7*T  +\
                    3*5.036111111111111e-7*T**2

    #Calculate the derivatives of the sun's right ascension (sun_RA) and declination (sun_dec)
    coeff_obliquity = np.sin( obliquity*deg )*np.tan( true_long*deg )   *np.cos( sun_RA*deg )**2
    coeff_true_long = np.cos( obliquity*deg )*np.cos( sun_RA*deg )**2 /  np.cos( true_long*deg )**2
    dsun_RA_dT      = coeff_true_long * dtrue_long_dT - coeff_obliquity * dobliquity_dT
         
    dsun_dec_dT   = np.cos(obliquity*deg)*np.sin( true_long*deg )/np.cos( sun_dec*deg ) * dobliquity_dT  + \
                    np.sin(obliquity*deg)*np.cos( true_long*deg )/np.cos( sun_dec*deg ) * dtrue_long_dT 

    #Calculate the derivatives of the eccentricity and radius of the Sun in its orbit
    de_dT         = - 0.000042037 - 2*0.0000001267*T
    dR_dT         = 1.000001018*( -2*e*de_dT/(1+e*np.cos(true_anom*deg)) - \
                                (1-e**2)/(1+e*np.cos(true_anom*deg))**2*( de_dT*np.cos(true_anom*deg) - \
                                                                       e*dtrue_anom_dT*deg*np.sin(true_anom*deg)))

    #Rescale distance to kilometers and kilometers/Julian centuries
    dsun_dist_dT = dR_dT*AU

    #Pack the parameters up and rescale to secs. Note: angle-related parameters are still in deg & deg/sec
    dT_dt                        = 1.0/(36525.0*86400.0)
    d_sun_RA                     = dsun_RA_dT   * dT_dt
    d_sun_dec                    = dsun_dec_dT  * dT_dt
    d_sun_dist                   = dsun_dist_dT * dT_dt 

    return sun_dist, sun_RA, sun_dec, d_sun_dist, d_sun_RA, d_sun_dec

