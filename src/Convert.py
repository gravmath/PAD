# -*- coding: utf-8 -*-
import datetime           as dt
import matplotlib.dates   as mdates
import numpy              as np
import scipy.interpolate  as interp
from spacepy import pycdf

import SunEph 
import equation_of_time

eclip_pole_GCI = np.array([0.0,-0.39777715575399,0.917482062146321])
dipole_ECEF    = np.array([-0.05752247,0.17291892,-0.98325491])

###############################################################################
def convert_GCI_to_GSE(epoch_dict):
    S, V  = SunEph.CalcSun_Low(epoch_dict)
    X_hat = S/np.sqrt(S.dot(S))
    Y     = np.cross(eclip_pole_GCI,X_hat)
    Y_hat = Y/np.sqrt(Y.dot(Y))
    Z_hat = np.cross(X_hat,Y_hat)
    
    return np.vstack((X_hat,Y_hat,Z_hat))

###############################################################################    
def convert_ECEF_to_GCI(epoch_dict):
    gha   = equation_of_time.calculate_GMST(epoch_dict)
    cos_g = np.cos(gha)
    sin_g = np.sin(gha)
    
    return np.array([[cos_g,-sin_g,0],[sin_g,cos_g,0],[0,0,1]])

###############################################################################    
def convert_GSE_to_GSM(epoch_dict):
    A_gse_gci     = convert_GCI_to_GSE(epoch_dict)
    S, V          = SunEph.CalcSun_Low(epoch_dict)
    D             = convert_ECEF_to_GCI(epoch_dict).dot(dipole_ECEF)
    X_hat         = S/np.sqrt(S.dot(S))
    Y             = np.cross(S,D)
    Y_hat         = Y/np.sqrt(Y.dot(Y))
    Z_hat         = np.cross(X_hat,Y_hat)
    
    X_hat_GSE     = A_gse_gci[0,:]
    Y_hat_GSE     = A_gse_gci[1,:]
    Z_hat_GSE     = A_gse_gci[2,:]
    
    A00           = np.dot(X_hat_GSE,X_hat)
    A01           = np.dot(Y_hat_GSE,X_hat)
    A02           = np.dot(Z_hat_GSE,X_hat)
    A10           = np.dot(X_hat_GSE,Y_hat)
    A11           = np.dot(Y_hat_GSE,Y_hat)
    A12           = np.dot(Z_hat_GSE,Y_hat)
    A20           = np.dot(X_hat_GSE,Z_hat)
    A21           = np.dot(Y_hat_GSE,Z_hat)
    A22           = np.dot(Z_hat_GSE,Z_hat)

    return np.array([[A00,A01,A02],[A10,A11,A12],[A20,A21,A22]])

###############################################################################
def convert_GSM_to_ECEF(epoch_dict):
    T_GCI_ECEF = convert_ECEF_to_GCI(epoch_dict)
    T_GSE_GCI  = convert_GCI_to_GSE(epoch_dict)
    T_GSM_GSE  = convert_GSE_to_GSM(epoch_dict)
    
    T_GSM_ECEF = T_GSM_GSE.dot(T_GSE_GCI.dot(T_GCI_ECEF))
    
    return T_GSM_ECEF.transpose()

###############################################################################    
def convert_GCI_to_SM(epoch_dict):
    S, V  = SunEph.CalcSun_Low(epoch_dict)
    D     = convert_ECEF_to_GCI(epoch_dict).dot(dipole_ECEF)
    Z_hat = -D/np.sqrt(D.dot(D))
    Y     = np.cross(S,D)
    Y_hat = Y/np.sqrt(Y.dot(Y))
    X_hat = np.cross(Y_hat,Z_hat)
    return np.vstack((X_hat,Y_hat,Z_hat))    

###############################################################################
def convert_GCI_to_GSM(epoch_dict):
    return convert_GSE_to_GSM(epoch_dict).dot(convert_GCI_to_GSE(epoch_dict))
    
def calc_LT(V,x0_val=12):
    x, y, z = V
    LT = np.arctan2(y,x)*180/np.pi/15
    return (LT+x0_val) % 24    

###############################################################################
def convert_to_LM(B,U):
    B_norm = np.sqrt(B.dot(B))
    U_norm = np.sqrt(U.dot(U))
    z_LM   =  B/B_norm
    y      = np.cross(B,U)
    y_norm = np.sqrt(y.dot(y))
    y_LM   = y/y_norm
    x_LM   = np.cross(y_LM,z_LM)
    
    return np.vstack((x_LM,y_LM,z_LM))   

###############################################################################    
def construct_interpolants(cursor,fpi_prd1,obs,mode,descriptor,year,month,day):
    mquery = cursor.execute('Select ver, filename from mec_data where\
                                                obs        = "%s" and\
                                                mode       = "%s" and\
                                                descriptor = "%s" and\
                                                year       =  %s  and\
                                                month      =  %s  and\
                                                day        =  %s;' % \
                                         (obs,mode,descriptor,year,month,day))
    
    Re       = 6378.14
    mresults = mquery.fetchall()
    for mr in mresults:
        MEC_file = fpi_prd1+mr[1]
        MEC      = pycdf.CDF(MEC_file)
        mt       = mdates.date2num(MEC['Epoch'])
        mr_gsm   = np.asarray(MEC['%s_mec_r_gsm'%obs])/Re
        MEC.close()
        orbit_extent = np.max(np.abs(mr_gsm))
        if orbit_extent < 50.0*Re:
            x_gsm_spline = interp.splrep(mt,mr_gsm[:,0])
            y_gsm_spline = interp.splrep(mt,mr_gsm[:,1])
            z_gsm_spline = interp.splrep(mt,mr_gsm[:,2])
            return x_gsm_spline, y_gsm_spline, z_gsm_spline

    #if no MEC file works
    return False, False, False    
    
###############################################################################
def round_time(epoch, date_delta=dt.timedelta(minutes=1), to='average'):
    """
    Purpose:      Round a datetime object to a multiple of a timedelta
    epoch:        datetime.datetime object, no default
    dateDelta:    timedelta object, we round to a multiple of this, 
                  default 1 minute.
    Adapted from: http://stackoverflow.com/questions/3463930/
                  how-to-round-the-minute-of-a-datetime-object-python
    """
    round_to = date_delta.total_seconds()

    if epoch is None:
        epoch = dt.datetime.now()
    seconds = (epoch - epoch.min).seconds

    if to == 'up':
        # // is a floor division
        rounding = (seconds + round_to) // round_to * round_to
    elif to == 'down':
        rounding = seconds // round_to * round_to
    else:
        rounding = (seconds + round_to / 2) // round_to * round_to

    return epoch + dt.timedelta(0, rounding - seconds, -epoch.microsecond)
    