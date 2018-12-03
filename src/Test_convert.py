#A simple comment - that can be removed
import Convert
import datetime as dt
import numpy as np
#import pdb

epoch_dict          = {}
epoch_dict['year']  = 2015
epoch_dict['month'] = 12
epoch_dict['day']   = 1
epoch_dict['hour']  = 0
epoch_dict['min']   = 0
epoch_dict['sec']   = 0

epoch_obj = dt.datetime(2015,12,1,0,0,0)

#from the MEC file 
#'Z:/data/ftp/mms1/mec/srvy/l2/epht89d/2015/12/
#    mms1_mec_srvy_l2_epht89d_20151201_v2.1.0.cdf'
r_GCI = np.array([-39232.5166153,-38027.67374968,-19823.08526458])
v_GCI = np.array([-0.14577713,-1.80374737,-0.98742514])

r_GSM = np.array([54272.04901561,-20525.29361705,-3392.05804342])
v_GSM = np.array([1.9557809,0.626757,-0.17849566])

r_GSM_Convert = Convert.convert_GCI_to_GSM(epoch_obj).dot(r_GCI)
v_GSM_Convert = Convert.convert_GCI_to_GSM(epoch_obj).dot(v_GCI)

mag_r_GSM         = np.sqrt( r_GSM.dot(r_GSM) )
mag_v_GSM         = np.sqrt( v_GSM.dot(v_GSM) )
mag_r_GSM_Convert = np.sqrt( r_GSM_Convert.dot(r_GSM_Convert) )
mag_v_GSM_Convert = np.sqrt( v_GSM_Convert.dot(v_GSM_Convert) )

diff_r = r_GSM - r_GSM_Convert
diff_v = v_GSM - v_GSM_Convert

print "Looking at position conversion  "
print "The difference vector is:       ", diff_r
print "It's magnitude is:              ", np.sqrt( diff_r.dot(diff_r) )
print "It's percentage is:             ", np.sqrt( diff_r.dot(diff_r) )/mag_r_GSM*100.0
print "The difference in magnitudes is:", mag_r_GSM - mag_r_GSM_Convert
print "The angle between them is (deg):", np.arccos(r_GSM.dot(r_GSM_Convert)/mag_r_GSM/mag_r_GSM_Convert)*180/np.pi

print "Looking at velocity conversion  "
print "The difference vector is:       ", diff_v
print "It's magnitude is:              ", np.sqrt( diff_v.dot(diff_v) )
print "It's percentage is:             ", np.sqrt( diff_v.dot(diff_v) )/mag_v_GSM*100.0
print "The difference in magnitudes is:", mag_v_GSM - mag_v_GSM_Convert
print "The angle between them is (deg):", np.arccos(v_GSM.dot(v_GSM_Convert)/mag_v_GSM/mag_v_GSM_Convert)*180/np.pi

