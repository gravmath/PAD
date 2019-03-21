import datetime as dt
import dateutil.tz as tz
import h5py
import matplotlib.dates as mdates
import numpy    as np
from spacepy import pycdf

#configuration - time deltas (burst)
bpsd_f_delta     = 2.5
epsd_f_delta     = 2.5
dce_b_delta      = 0.00015
scpot_f_delta    = 0.001
fgm_b_delta      = 0.008
fgm_s_delta      = 1.0
edist_b_delta    = 0.031
emoms_b_delta    = 0.031
idist_b_delta    = 0.151
imoms_b_delta    = 0.151
edist_f_delta    = 5.0
emoms_f_delta    = 5.0
idist_f_delta    = 5.0
imoms_f_delta    = 5.0
mec_s_delta      = 45.0

#translations
# I believe that 'special' and 'eval' are no longer any different but the 
#naming has not been changed so as to not fix what ain't broke
#

#dsp (2) both fast
bpsd_f_translation = {'epochs':['Epoch','null'],
                      'freqs' :['"%s_b_freq" %(obs,)','special'],
                      'bpsd'  :['"%s_dsp_bpsd_omni_fast_l2" % (obs,)','eval']}

epsd_f_translation = {'epochs':['Epoch','null'],
                      'freqs' :['"%s_e_freq" %(obs,)','special'],
                      'epsd'  :['"%s_dsp_epsd_omni" % (obs,)','eval']}

#edp (2) 
dce_b_translation = {'epochs':['"%s_edp_epoch_brst_l2" % (obs,)','eval'],
                     'Egse'  :['"%s_edp_dce_gse_brst_l2" % (obs,)','eval']}

scpot_f_translation = {'epochs' : ['"%s_edp_epoch_fast_l2" % (obs,)','eval'],
                       'scpot'  : ['"%s_edp_scpot_fast_l2" % (obs,)','eval']}


#fgm (2) brst and fast
fgm_b_translation = {'epochs':['Epoch','null'],
                     'Bbcs'  :['"%s_fgm_b_bcs_brst_l2" % (obs,)','eval'],
                     'Bgse'  :['"%s_fgm_b_gse_brst_l2" % (obs,)','eval'],
                     'Bgsm'  :['"%s_fgm_b_gsm_brst_l2" % (obs,)','eval']}
                   
fgm_s_translation = {'epochs':['Epoch','null'],
                     'Bbcs'  :['"%s_fgm_b_bcs_srvy_l2" % (obs,)','eval'],
                     'Bgse'  :['"%s_fgm_b_gse_srvy_l2" % (obs,)','eval'],
                     'Bgsm'  :['"%s_fgm_b_gsm_srvy_l2" % (obs,)','eval']}                   

#fpi (8) e/i dist/moms brst/fast
edist_b_translation = {'epochs'    :['Epoch','null'],
                       'dist'      :['"%s_%s_dist_brst" % (obs,"des")','eval'],
                       'disterr'   :['"%s_%s_disterr_brst" % (obs,"des")','eval'],
                       'ergs'      :['"%s_%s_energy_brst" % (obs,"des")','eval'],
                       'parity'    :['"%s_%s_steptable_parity_brst" % (obs,"des")','eval'],
                       'phis'      :['"%s_%s_phi_brst" % (obs,"des")','eval'],
                       'start_dphi':['"%s_%s_startdelphi_count_brst" % (obs,"des")','eval'],
                       'thetas'    :['"%s_%s_theta_brst" % (obs,"des")','special']}
               
idist_b_translation = {'epochs'    :['Epoch','null'],
                       'dist'      :['"%s_%s_dist_brst" % (obs,"dis")','eval'],
                       'disterr'   :['"%s_%s_disterr_brst" % (obs,"dis")','eval'],
                       'ergs'      :['"%s_%s_energy_brst" % (obs,"dis")','eval'],
                       'phis'      :['"%s_%s_phi_brst" % (obs,"dis")','eval'],
                       'start_dphi':['"%s_%s_startdelphi_count_brst" % (obs,"dis")','eval'],
                       'thetas'    :['"%s_%s_theta_brst" % (obs,"dis")','special']}                   
               
emoms_b_translation = {'epochs' :['Epoch','null'],
                       'anti'   :['"%s_%s_energyspectr_anti_brst" % (obs,"des")','eval'],
                       'bulk_vs':['"%s_%s_bulkv_gse_brst" % (obs,"des")','eval'],
                       'ergs'   :['"%s_%s_energy_brst" % (obs,"des")','eval'],
                       'heats'  :['"%s_%s_heatq_gse_brst" % (obs,"des")','eval'],
                       'num_den':['"%s_%s_numberdensity_brst" % (obs,"des")','eval'],
                       'omnis'  :['"%s_%s_energyspectr_omni_brst" % (obs,"des")','eval'],
                       'par'    :['"%s_%s_energyspectr_par_brst" % (obs,"des")','eval'],
                       'perp'   :['"%s_%s_energyspectr_perp_brst" % (obs,"des")','eval'],                     
                       'pres_s' :['"%s_%s_prestensor_gse_brst" % (obs,"des")','eval'],
                       'T_par'  :['"%s_%s_temppara_brst" % (obs,"des")','eval'],
                       'T_perp' :['"%s_%s_tempperp_brst" % (obs,"des")','eval'],
                       'T_s'    :['"%s_%s_temptensor_gse_brst" % (obs,"des")','eval'] }

            
imoms_b_translation = {'epochs' :['Epoch','null'],
                       'bulk_vs':['"%s_%s_bulkv_gse_brst" % (obs,"dis")','eval'],
                       'ergs'   :['"%s_%s_energy_brst" % (obs,"dis")','eval'],
                       'heats'  :['"%s_%s_heatq_gse_brst" % (obs,"dis")','eval'],
                       'num_den':['"%s_%s_numberdensity_brst" % (obs,"dis")','eval'],
                       'omnis'  :['"%s_%s_energyspectr_omni_brst" % (obs,"dis")','eval'],
                       'pres_s' :['"%s_%s_prestensor_gse_brst" % (obs,"dis")','eval'],
                       'T_par'  :['"%s_%s_temppara_brst" % (obs,"dis")','eval'],
                       'T_perp' :['"%s_%s_tempperp_brst" % (obs,"dis")','eval'],
                       'T_s'    :['"%s_%s_temptensor_gse_brst" % (obs,"dis")','eval'] }  
                     
edist_f_translation = {'epochs'    :['Epoch','null'],
                       'dist'      :['"%s_%s_dist_fast" % (obs,"des")','eval'],
                       'disterr'   :['"%s_%s_disterr_fast" % (obs,"des")','eval'],
                       'ergs'      :['"%s_%s_energy_fast" % (obs,"des")','eval'],
                       'phis'      :['"%s_%s_phi_fast" % (obs,"des")','eval'],
                       'start_dphi':['"%s_%s_startdelphi_count_fast" % (obs,"des")','eval'],
                       'thetas'    :['"%s_%s_theta_fast" % (obs,"des")','special']}
               
idist_f_translation = {'epochs'    :['Epoch','null'],
                       'dist'      :['"%s_%s_dist_fast" % (obs,"dis")','eval'],
                       'disterr'   :['"%s_%s_disterr_fast" % (obs,"dis")','eval'],
                       'ergs'      :['"%s_%s_energy_fast" % (obs,"dis")','eval'],
                       'phis'      :['"%s_%s_phi_fast" % (obs,"dis")','eval'],
                       'start_dphi':['"%s_%s_startdelphi_count_fast" % (obs,"dis")','eval'],
                       'thetas'    :['"%s_%s_theta_fast" % (obs,"dis")','special']}                   
               
emoms_f_translation = {'epochs' :['Epoch','null'],
                       'anti'   :['"%s_%s_energyspectr_anti_fast" % (obs,"des")','eval'],
                       'bulk_vs':['"%s_%s_bulkv_gse_fast" % (obs,"des")','eval'],
                       'ergs'   :['"%s_%s_energy_fast" % (obs,"des")','eval'],
                       'heats'  :['"%s_%s_heatq_gse_fast" % (obs,"des")','eval'],
                       'num_den':['"%s_%s_numberdensity_fast" % (obs,"des")','eval'],
                       'omnis'  :['"%s_%s_energyspectr_omni_fast" % (obs,"des")','eval'],
                       'par'    :['"%s_%s_energyspectr_par_fast" % (obs,"des")','eval'],
                       'perp'   :['"%s_%s_energyspectr_perp_fast" % (obs,"des")','eval'],                     
                       'pres_s' :['"%s_%s_prestensor_gse_fast" % (obs,"des")','eval'],
                       'T_par'  :['"%s_%s_temppara_fast" % (obs,"des")','eval'],
                       'T_perp' :['"%s_%s_tempperp_fast" % (obs,"des")','eval'],
                       'T_s'    :['"%s_%s_temptensor_gse_fast" % (obs,"des")','eval'] }
            
imoms_f_translation = {'epochs' :['Epoch','null'],
                       'bulk_vs':['"%s_%s_bulkv_gse_fast" % (obs,"dis")','eval'],
                       'ergs'   :['"%s_%s_energy_fast" % (obs,"dis")','eval'],
                       'heats'  :['"%s_%s_heatq_gse_fast" % (obs,"dis")','eval'],
                       'num_den':['"%s_%s_numberdensity_fast" % (obs,"dis")','eval'],
                       'omnis'  :['"%s_%s_energyspectr_omni_fast" % (obs,"dis")','eval'],
                       'pres_s' :['"%s_%s_prestensor_gse_fast" % (obs,"dis")','eval'],
                       'T_par'  :['"%s_%s_temppara_fast" % (obs,"dis")','eval'],
                       'T_perp' :['"%s_%s_tempperp_fast" % (obs,"dis")','eval'],
                       'T_s'    :['"%s_%s_temptensor_gse_fast" % (obs,"dis")','eval'] }                       
                          
ecnts_translation = {'epochs'    :['Epoch','null'],
                     'cnts'      :['"%s_%s_brstSkyMap_cnts" % (obs,"des")','eval'] }                          
                         
#hpca (2)
#coming whenever
     
#mec (1) always fast
mec_s_translation = {'epochs'       : ['Epoch','null'],
                     'dipole_tilt'  : ['"%s_mec_dipole_tilt" % (obs,)','eval'],
                     'dst'          : ['"%s_mec_dst" % (obs,)','eval'],                   
                     'fieldline'    : ['"%s_mec_fieldline_type" % (obs,)','eval'],
                     'foot_n'       : ['"%s_mec_pfn_geod_latlon" % (obs,)','eval'],
                     'foot_s'       : ['"%s_mec_pfs_geod_latlon" % (obs,)','eval'],                   
                     'gse_pos'      : ['"%s_mec_r_gse" % (obs,)','eval'],
                     'gsm_pos'      : ['"%s_mec_r_gsm" % (obs,)','eval'],
                     'Kp'           : ['"%s_mec_kp" % (obs,)','eval'],
                     'Lshell'       : ['"%s_mec_l_dipole" %(obs,)','eval'],
                     'losscone_n'   : ['"%s_mec_loss_cone_angle_n" %(obs,)','eval'],
                     'losscone_s'   : ['"%s_mec_loss_cone_angle_s" %(obs,)','eval'],
                     'mag_model_ext': ['"%s_mec_ext_model" %(obs,)','eval'],
                     'mag_model_int': ['"%s_mec_int_model" %(obs,)','eval'],
                     'mlat'         : ['"%s_mec_mlat" % (obs,)','eval'],                   
                     'mlt'          : ['"%s_mec_mlt" % (obs,)','eval'],                   
                     'sm_pos'       : ['"%s_mec_r_sm" % (obs,)','eval']}
                         

###############################################################################	
#
#
# New function approach
# 
#
###############################################################################    
    
#############################################################################     
def make_data_dict_via_translation(name,translation):
    """Helper function for creating a custom data dictionary 
       by using a translation
       
       Arguments:
          name:        a string used to identify the name
          translation: translation structure used to structure the munge based 
                       on the CDF file

       Returns:
           The structured dictionary which is the atom for a munge
       
       Example use:  
           my_dict = make_data_dict_via_translation('des_moms',des_moms_translation)
              
       Note:  Should only be used from one of the core Munger functions
       """
       
    #make baseline structure for the munge
    A = {'name':name,
         'start':False,
         'stop':False,
         'num_segs':0}
    
    for k in translation.keys():
        A[k] = np.asarray([])
        
    return A    
    
#############################################################################
def make_munge_via_translation(obs,type,delta,file_list,translation,fn="None"):
    """Core function that uses an instrument-tailored hash to make a generic
      munge of all the data fed it in file_list
      
       Arguments:
          obs:         'mms1', 'mms2', 'mms3', or 'mms4'
          type:        string to identify and to name
          delta:       time delta in seconds to determine when one series ends
                       and the next begins
          file_list:   list of fully qualified filenames to be opened and looted 
                       of data
          translation: translation structure used to structure the munge based 
                       on the CDF file
          fn:          optional filename to redirect the output

       Returns:
           The structured munge - a list of dictionary data atoms
       
       Example use:  
           fgm_munge = Munger.make_munge_via_translation('mms1','fgm',Munger.fgm_s_delta,m1f['fgm_s'],Munger.fgm_s_translation) 
              
       Note:  The workhorse
      """
    #import pdb; pdb.set_trace()
    
    print "****************************************"
    print "Munging {type:} on obs {obs:} at {now:}".format(type=type,obs=obs,now=dt.datetime.now())
    if fn != "None":
        fh = open(fn,'a')
        fh.write("****************************************\n")
        fh.write("Munging %s on obs %s at %s\n" % (type,obs,dt.datetime.now()) )        
        
        
    B       = []
    counter = 0
    segment = 1
    A       = make_data_dict_via_translation('{type:}_stride{counter:02d}'.format(type=type,counter=counter),translation)


    for f in file_list:
        #open the cdf
        try:
            cdf  = pycdf.CDF(f)
        except:
            print "pycdf.CDF can't open {file:}!".format(file=f)
            if fn != "None":
                fh.write("pycdf.CDF can't open {file:}!\n".format(file=f))
            return False
        
        #make the temporary data dictionary and load
        temp = make_data_dict_via_translation('temp',translation)
        for k in translation:
            if translation[k][1] == 'null':
                temp[k] = np.asarray(cdf[translation[k][0]][:])
            if translation[k][1] == 'eval':
                temp[k] = np.asarray(cdf[eval(translation[k][0])][:])
            if translation[k][1] == 'special':
                temp[k] = np.asarray(cdf[eval(translation[k][0])][:])
        
        #close the cdf
        cdf.close()

        #report the time span of the first segment (note the temporary way of disabling the choice in printing to the screen)
        print 'segment %s - start: %s stop: %s' % (segment,temp['epochs'][0],temp['epochs'][-1])
        if fn != 'None':
            fh.write("segment %s - start: %s stop: %s\n" % (segment,temp['epochs'][0],temp['epochs'][-1]))
                    
        
        #pack it into the munge
        if A['num_segs'] == 0:
            print 'fresh segment - first stride'
            if fn != 'None':
               fh.write("fresh segment - first stride\n")
            A['start']      = temp['epochs'][0]
            A['stop']       = temp['epochs'][-1]
            for k in translation.keys():
                A[k] = temp[k]
            A['num_segs'] += 1
        elif (temp['epochs'][0] - A['stop']).total_seconds() < delta and (temp['epochs'][0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            print 'adjacency underway'
            if fn != 'None':
                fh.write("adjacency underway\n")
            A['stop']       = temp['epochs'][-1]
            for k in translation.keys():
                num_dim = len(temp[k].shape)
                if num_dim == 1:
                    A[k] = np.hstack((A[k],temp[k]))
                if num_dim > 1:
                    A[k] = np.vstack((A[k],temp[k]))
            A['num_segs'] += 1            
        elif (temp['epochs'][0] - A['stop']).total_seconds() < delta and (temp['epochs'][0] - A['stop']).total_seconds() < 0.0:
            #overlaps exist - assume that the coincident measurements are identical
            print 'overlap detected - dealing with it'
            if fn != 'None':
                fh.write("overlap detected - dealing with it\n")
            time_delta      = np.array([(temp['epochs'][j] - A['stop']).total_seconds() for j in range(len(temp['epochs']))])
            new_points      = np.where(time_delta > 0.0)
            if temp['epochs'][-1] > A['stop']:
                A['stop']   = temp['epochs'][-1]
            for k in translation.keys():
                num_dim = len(temp[k].shape)
                if num_dim == 1 and temp[k].shape == temp['epochs'].shape:
                    A[k] = np.hstack((A[k],temp[k][new_points]))
                elif num_dim == 1 and temp[k].shape != temp['epochs'].shape:
                    #special case where an oeverlap exists but the data
                    #are non-record varying
                    A[k] = np.hstack((A[k],temp[k]))
                if num_dim > 1:
                    A[k] = np.vstack((A[k],temp[k][new_points]))
            A['num_segs']  += 1            
        elif (temp['epochs'][0] - A['stop']).total_seconds() > delta:
            #break in the epochs so that this is a new segment
            print 'break in adjacency - new stride'
            if fn != 'None':
                fh.write("break in adjacency - new stride\n")
            B.append(A)
            counter += 1
            A          = make_data_dict_via_translation('{type:}_stride{counter:02d}'.format(type=type,counter=counter),translation)
            A['start'] = temp['epochs'][0]
            A['stop']  = temp['epochs'][-1]
            for k in translation.keys():
                A[k] = temp[k]    
            A['num_segs'] += 1
        segment += 1
            
    B.append(A)
    print 'Munged %s series for %s on %s!\n' % (len(B),type,obs)
    print "Finished at %s\n" % (dt.datetime.now(),)
    print "****************************************\n"

    if fn != 'None':
        fh.write("Munged %s series for %s on %s!\n" % (len(B),type,obs))
        fh.write("Finished at %s\n" % (dt.datetime.now(),))
        fh.write("****************************************\n")
        fh.close()
    return B
    
#############################################################################        
def interpolate_to_epoch(source_epoch,source_data,target_epoch):
    """Helper function that linearly interpolates a single time series 
      (source_epoch and source_data) to the time knots of another series
      (target_epoch)
      
       Arguments:
          source_epoch: a numpy.array of datetime objects
          source_data:  a numpy.array of data values
          target_epoch: a numpy.array of datatime object

       Returns:
           target_data from the interpolation appropriately masked
       
       Example use:  
           mimic_munge[j][k][:,m] = interpolate_to_epoch(source_epoch,source_data,target_epoch) 
              
       Note:       
           Assumes that the time-knots are python datetime objects so that it
           can use matplotlib's date2num function
      
           No error checking is performed
      
           Should only be called from within Munger
      """
       
    import matplotlib.dates  as mdates
    import scipy.interpolate as interp

    #construct the interpolant
    interpolant = interp.interp1d(mdates.date2num(source_epoch),source_data,kind='slinear')
    
    #check the epoch range
    good_indices              = np.where(np.logical_and(target_epoch > source_epoch[0],target_epoch < source_epoch[-1]))
    target_data               = np.ones(target_epoch.shape)*np.nan
    target_data[good_indices] = interpolant(mdates.date2num(target_epoch[good_indices]))
    
    return np.ma.masked_invalid(target_data)    
    
#############################################################################    
def make_mimic_munge(munge):
    """Helper function that makes a munge identical in 
       structure as an existing munge but lacking data
       
       Arguments:
          munge:  the munge to be mimicked 
          
       Returns:
           mimic: the mimicked munge
       
       Example use:  
           mimic_munge = make_mimic_munge(source_munge)
              
       Note:       
           Should only be called from within Munger
       """

       
    import copy
    
    #determine the number of strides
    num_strides = len(munge)
    
    mimic = []
    for j in range(num_strides):
        #determine the flavor
        if 'dce' in munge[0]['name']:
            print 'dce flavored'
            mimic.append(make_data_dict_via_translation('dce_mimic',dce_b_translation))
        if 'dist' in munge[0]['name']:
            print 'dist flavored'
            mimic.append(make_data_dict_via_translation('dist_mimic',edist_b_translation))
        if 'fgm' in munge[0]['name']:
            print 'fgm flavored'
            mimic.append(make_data_dict_via_translation('fgm_mimic',fgm_b_translation))
        if 'moms' in munge[0]['name']:
            print 'moms flavored'
            mimic.append(make_data_dict_via_translation('moms_mimic',emoms_b_translation))
        for k in munge[j].keys():
            if type(munge[j][k]) != type(np.array([])):
                mimic[j][k] = copy.deepcopy(munge[j][k])
 
    return mimic
    
#############################################################################   
def interpolate_munge_to_munge(source_munge,target_munge,fn="None"):
    """Core function designed to take two munges of the same form (i.e
       several burst segments from a given time period) and to linearly
       interpolate the data from the source one to the time-knots of the
       target.  
       
       In doing so, a mimicked source_munge is returned with time-knots 
       consistent with target_munge.
       
       The use case is to adapt fgm or dce data to des or dis time-knots
       for JdotE or curlometery or the like.
       
       Arguments:
          source_munge:  the munge with the raw data
          target_munge:  the munge with the epochs
          fn:            optional log file for workflow
          
       Returns:
           mimic_munge:  the mimicked munge with the source's data now resampled
                         to the target's epoch
       
       Example use:  
           efgm_munge = Munger.interpolate_munge_to_munge(fgm_munge,emoms_munge)
              
       Note:       
           The runner-up workhorse
       """
       
    import copy
    
    #import pdb; pdb.set_trace()

    print "****************************************"
    print "Interpolating %s to %s" % (source_munge[0]['name'],target_munge[0]['name'])
    #Open optional log file if fn is provided
    if fn != "None":
        fh = open(fn,'a')
        fh.write("****************************************\n")
        fh.write("Interpolating %s to %s \n" % (source_munge[0]['name'],target_munge[0]['name']))     
    
    if len(source_munge) != len(target_munge):
        print 'Direct adaptation of source_munge to target_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
        if fn != "None":
            fh.write("Direct adaptation of source_munge to target_munge not possible.")
            fh.write("Terminating with extreme prejudice!!!")
        return False
        
    #make a mimic of the source_munge
    mimic_munge = make_mimic_munge(source_munge)
    
    #determine the number of strides
    num_strides = len(source_munge)
    
    for j in range(num_strides):
        mimic_munge[j]['epochs'] = copy.deepcopy(target_munge[j]['epochs'])
        
    for j in range(num_strides):
        print "Interpolating stride %s of %s." % (j, num_strides)
        if fn != "None":
            fh.write("Interpolating stride %s of %s." % (j, num_strides))
        source_epoch     = source_munge[j]['epochs']
        target_epoch     = target_munge[j]['epochs']
        num_target_times = len(target_epoch) 
        for k in source_munge[j].keys():
            if k != 'epochs' and type(source_munge[j][k]) == type(np.array([])):
                #determine number of components (scale, vector, array, tensor...)
                #assumes that the ordering is always time (k) vs. components
                data_shape        = list(source_munge[j][k].shape)
                component_shape   = tuple(data_shape[1:])
                #if len(component_shape) == 0:
                #    component_shape = (1)
                data_shape[0]     = num_target_times
                data_shape        = tuple(data_shape)
                num_dim           = len(component_shape)
                mimic_munge[j][k] = np.zeros(data_shape)    
                if num_dim == 0: #scalar series
                    source_data = source_munge[j][k][:]
                    mimic_munge[j][k][:] = interpolate_to_epoch(source_epoch,source_data,target_epoch)
                if num_dim == 1: #vector series
                    for m in range(component_shape[0]):
                        source_data = source_munge[j][k][:,m]
                        mimic_munge[j][k][:,m] = interpolate_to_epoch(source_epoch,source_data,target_epoch)
                if num_dim == 2: #rank-2 tensor series
                    for m in range(component_shape[0]):
                        for n in range(component_shape[1]):
                            source_data = source_munge[j][k][:,m,n]
                            mimic_munge[j][k][:,m,n] = interpolate_to_epoch(source_epoch,source_data,target_epoch)
                if num_dim == 3: #distribution/skymap
                    for m in range(component_shape[0]):
                        for n in range(component_shape[1]):
                            for p in range(component_shape[2]):
                                source_data = source_munge[j][k][:,m,n,p]
                                mimic_munge[j][k][:,m,n,p] = interpolate_to_epoch(source_epoch,source_data,target_epoch)
                mimic_munge[j][k] = np.ma.masked_invalid(mimic_munge[j][k])
                
    print "All done with the interpolation of %s to %s at %s" % (source_munge[0]['name'],target_munge[0]['name'],dt.datetime.now())
    print "*************************************************"
    if fn != "None":
        fh.write("All done with the interpolation of %s to %s at %s" % (source_munge[0]['name'],target_munge[0]['name'],dt.datetime.now()))
        fh.write("*************************************************\n")
        fh.close()
                
    return mimic_munge

#############################################################################      
def average_time_series_to_time_knots(source_tknots,source_data,target_tknots):
    """Helper function designed to average a source time series to a target
       time knots.
       
       Arguments:
          source_tknots: times of the time series
          source_data:   corresponding data
          target_tknots: target times
          
       Returns:
           mimic_data:   the mimicked munge with the source's data now averaged
                         to the target's epoch
       
       Example use:  
           mimciked_efield = Munger.average_time_series_to_time_knots(E_times,E_vals,ion_times)
              
       Note:       
           Far too slow - how to speed it up????
       """
    num_tknots            = len(target_tknots)
    mimic_series_shape    = list(source_data.shape)
    mimic_series_shape[0] = num_tknots
    mimic_series_shape    = tuple(mimic_series_shape)
    mimic_value           = np.nan*np.ones(mimic_series_shape)
    start_k               = np.where(target_tknots > source_tknots[0] )[0][0]
    stop_k                = np.where(target_tknots < source_tknots[-1])[0][-2]
    curr_delta            = target_tknots[1] - target_tknots[0]

    #print 'Averaging %s tknots' % (num_tknots)
    for k in range(start_k,stop_k):
        #if k % 100 == 0:
            #print '*',
        curr_pts = np.logical_and(source_tknots >= target_tknots[k],\
                                  source_tknots <= target_tknots[k+1])
        mimic_value[k] = np.average(source_data[curr_pts],axis=0)
    #complete the final knot
    k = stop_k
    curr_pts = np.logical_and(source_tknots >= target_tknots[k],\
                              source_tknots <= target_tknots[k]+curr_delta)
    mimic_value[k] = np.average(source_data[curr_pts],axis=0)

    return mimic_value    
    
#############################################################################   
def adapt_munge_to_munge(source_munge,target_munge,fn="None"):
    """Core function designed to take two munges of the same form (i.e
       several burst segments from a given time period) and to average over the 
       time periods between the tknots of the target_munge.  
       
       In doing so, a mimicked source_munge is returned with time-knots 
       consistent with target_munge.
       
       The use case is to adapt fgm or dce data to des or dis time-knots
       for JdotE or curlometery or the like.
       
       Arguments:
          source_munge:  the munge with the raw data
          target_munge:  the munge with the epochs
          fn:            optional output log filename
          
       Returns:
           mimic_munge:  the mimicked munge with the source's data now resampled
                         to the target's epoch
       
       Example use:  
           efgm_munge = Munger.adapt_munge_to_munge(fgm_munge,emoms_munge)
              
       Note:       
           The runner-up workhorse
       """
       
    import copy
    
    #import pdb; pdb.set_trace()
    
    print "****************************************"
    print "Adapting %s to %s \n" % (source_munge[0]['name'],target_munge[0]['name'])
    #Open optional log file if fn is provided
    if fn != "None":
        fh = open(fn,'a')
        fh.write("****************************************\n")
        fh.write("Adapting %s to %s \n" % (source_munge[0]['name'],target_munge[0]['name']))     
    
    
    if len(source_munge) != len(target_munge):
        print 'Direct adaptation of source_munge to target_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
        if fn != "None":
            fh.write("Direct adaptation of source_munge to target_munge not possible.")
            fh.write("Terminating with extreme prejudice!!!")
        return False
        
    #make a mimic of the source_munge
    mimic_munge = make_mimic_munge(source_munge)
    
    #determine the number of strides
    num_strides = len(source_munge)
    
    for j in range(num_strides):
        mimic_munge[j]['epochs'] = copy.deepcopy(target_munge[j]['epochs'])
        
    for j in range(num_strides):
        source_epoch     = source_munge[j]['epochs']
        target_epoch     = target_munge[j]['epochs']
        num_target_times = len(target_epoch) 
        for k in source_munge[j].keys():
            if k != 'epochs' and type(source_munge[j][k]) == type(np.array([])):
                source_data = source_munge[j][k]
                mimic_munge[j][k] = average_time_series_to_time_knots(source_epoch,source_data,target_epoch)
                
    print "All done with the adaptation of %s to %s at %s" % (source_munge[0]['name'],target_munge[0]['name'],dt.datetime.now())
    print "*************************************************"
    if fn != "None":
        fh.write("All done with the adaptation of %s to %s at %s" % (source_munge[0]['name'],target_munge[0]['name'],dt.datetime.now()))
        fh.write("*************************************************\n")
        fh.close()
                
                
    return mimic_munge    
    
#############################################################################   
def adjust_epoch_by_delta(munge,delta):
    """Helper function designed to adjust the epochs in a given munge by
       a certain time delta; either to handle timing errors or to account
       for 'start time' versus 'center time' in the observations (e.g. FPI).
       
       Arguments:
          munge:    the munge to be adjusted
          delta:    the adjustment time (in seconds)
          
       Returns:
           nothing per se - the munge's epochs are now adjusted
       
       Example use:  
           adjust_epoch_by_delta(emoms_munge_f,4.5/2)
              
       Note:       
           Should be called before adaptation is done
       """    
       
    #create the time_delta
    my_time_delta = dt.timedelta(seconds=delta) 
    
    #determine the number of strides
    num_strides = len(munge)
    
    for N in range(num_strides):
        munge[N]['epochs'] = munge[N]['epochs'] + my_time_delta
        
        

############################################################################# 
def save_munge(munge,name,obs,h5file,fn="None"):
    """Main function designed to save the munges to an H5 file for later
    retrieval and modification.  
    
    possible:  This only works on memory-entrant data 
    from a set of CDFs.  IT IS NOT MEANT TO BE USED TO UPDATE THE H5.
    Proven true on 2/7/19 - the current workaround is to delete the data set to 
    be updated by hand and then rerun Munger.  But it would be better to 
    do this within the code.
       
       Arguments:
          munge:       the munge to be saved
          name:        name of the munge in the H5 file - need not be the 
                       same as the variable name
          obs:         string with 'mms1', 'mms2', ...
          h5file:      pointer to the h5file
          fn:          optional file for workflow output
          
       Returns:
           True if successful or a False and a printed error message 
           otherwise.
       
       Example use:  
           save_munge(munge,'skippy',my_h5)
              
       Note:       
           none
       """    
    #import pdb; pdb.set_trace()
    
    print "****************************************"
    print "Saving {name:} for obs {obs:} at {now:}".format(name=name,obs=obs,now=dt.datetime.now())
    #open the optional output file
    if fn != "None":
        fh = open(fn,'a')
        fh.write("****************************************\n")
        fh.write("Saving %s for obs %s at %s \n" % (name,obs,dt.datetime.now()) )        

    if type(munge) != list:
        print "Something has gone awry!  Given munge is not a list!"
        if fn != "None":
            fh.write("Something has gone awry!  Given munge is not a list!\n")
            fh.write("Done trying to save at {now:}".format(now=dt.datetime.now()))
            fh.write("****************************************\n")
            fh.close()
            return False            
        
    #determine the number of strides
    num_strides = len(munge)

    #loop over the strides and the keys to store the munge
    transaction_count = 0    
    for stride in range(num_strides):
        for k in munge[stride].keys():
            #increment the transaction count as a metric for data elements/series
            transaction_count += 1
            
            #first construct the dictionary structure of the portion of the munge
            #to be saved
            h5path     = '/{obs:}/{name:}/stride{stride:02d}/{key}'.format(obs=obs,name=name,stride=stride,key=k)
            
            #then determine the type of data being stored
            type_of_k = type(munge[stride][k])
            
            #finally, initially assume that the data can't be found in the file
            #and then try to find the data in the file and, if found, set the 
            #flag 'in_file' appropriately
            in_file   = False
            try:
                data  = h5file[h5path]
                in_file = True
            except:
                pass
            
            if in_file == False:
                #No corresponding data in the h5
                try:
                    #Try to save the data in whatever manner is appropriate to its type
                    
                    #handle the non-date metadata appropriately (currently the string name and the integer num_segs)
                    if type_of_k == str or type_of_k == int:
                        h5file.create_dataset(h5path,data=munge[stride][k])
                    #handle the date metadata appropriately (currently stride start and stop)
                    elif type_of_k == dt.datetime:
                        temp_epoch = mdates.date2num(munge[stride][k].replace(tzinfo=tz.tzutc()))
                        h5file.create_dataset(h5path,data=temp_epoch)
                    #if the data are contained in a numpy array, then store/update as appropriate 
                    else:
                        #It's an array of data
                        
                        #find the length of the array
                        stride_len = len(munge[stride][k])
                        
                        if munge[stride][k].dtype == 'object' or munge[stride][k].dtype == 'O':
                            #Its an array of epochs
                        
                            #convert epoch data to numbers via mdates.date2num before storing in the h5
                            epoch_data = mdates.date2num([ munge[stride][k][i].replace(tzinfo=tz.tzutc()) for i in range(stride_len)])
                            h5file.create_dataset(h5path,data=epoch_data)
                        else:
                            #Its numeric data of some kind
                            
                            #otherwise just store the data
                            h5file.create_dataset(h5path,data=(munge[stride][k]))
                    print 'Succeeded in saving stride{stride:02d}/{key:}'.format(stride=stride,key=k)
                    if fn != "None":
                       fh.write('Succeeded in saving stride{stride:02d}/{key:}\n'.format(stride=stride,key=k))
                except:
                    #The attempt to save the data failed miserably
                    #
                    
                    print 'Failed on saving stride{stride:02d}/(key:}'.format(stride=stride,key=k)
                    if fn != "None":
                        fh.write('Failed on saving stride{stride:02d}/{key:}\n'.format(stride=stride,key=k))
            
            elif in_file == True:
                #Data exists in the h5
                try:
                    #Try to update the data
                
                    if type_of_k == str or type_of_k == int:
                        data[...] = munge[stride][k]
                    elif type_of_k == dt.datetime:
                        temp_epoch = mdates.date2num(munge[stride][k].replace(tzinfo=tz.tzutc()))
                        data[...] = temp_epoch
                    else:
                        #It's an array of data 
                        
                        #find the length of the array
                        stride_len = len(munge[stride][k])
                        
                        if munge[stride][k].dtype == 'object' or munge[stride][k].dtype == 'O':
                            #Its an array of epochs
                        
                            #convert epoch data to numbers via mdates.date2num before storing in the h5
                            epoch_data = mdates.date2num([ munge[stride][k][i].replace(tzinfo=tz.tzutc()) for i in range(stride_len)])
                            data[...] = epoch_data
                        else:
                            #Its numeric data of some kind
                            
                            #otherwise just store the data
                            data[...] = munge[stride][k]
                    print 'Succeeded in updating stride{stride:02d}/{key:}'.format(stride=stride,key=k)
                    if fn != "None":
                        fh.write('Succeeded in updating stride{stride:02d}/{key:}\n'.format(stride=stride,key=k))
                except:
                    #Attempt to update the data failed miserably
                    if fn == "None":
                        print 'Failed on updating stride{stride:02d}/{key:}'.format(stride=stride,key=k)
                    else:
                        fh.write('Failed on updating stride{stride:02d}/{key:}\n'.format(stride=stride,key=k))

    print "Done with saving {num:} transactions at {now}".format(num=transaction_count,now=dt.datetime.now())
    print "****************************************"    
    if fn != "None":
        fh.write('Done with saving {num:} transactions at {now}\n'.format(num=transaction_count,now=dt.datetime.now()))
        fh.write("****************************************\n")
        fh.close()
        
    return True
    
#############################################################################     
def load_munge(name,obs,h5file,fn="None"):
    """Main function designed to load a munges from an H5 file for later
    saving and modification.
       
       Arguments:
          name:        the munge to be retrieved 
          obs:         string with 'mms1', 'mms2', ... 
          h5file:      pointer to the h5file
          fn:          optional file for workflow output
          
       Returns:
           munge if successful or a False and a printed error message 
           otherwise.
       
       Example use:  
           my_munge = load_munge('skippy','mms1',my_h5)
              
       Note:       
           none
       """ 
     
    #import pdb; pdb.set_trace()

    #open the optional output file
    print "****************************************"
    print "Loading %s for obs %s at %s" % (name,obs,dt.datetime.now())

    if fn != "None":
        fh = open(fn,'a')
        fh.write("****************************************\n")
        fh.write("Loading %s for obs %s at %s \n" % (name,obs,dt.datetime.now()) )    
    
    #check to make sure the name is in the file
    #by retrieving the number of strides
    try:
        groups      = h5file['/%s/%s'%(obs,name,)].items()
        num_strides = len(groups)
        print 'found %s strides' % (num_strides,)
        if fn != "None":
            fh.write("found %s strides\n" % (num_strides,))
    except:
        print "the h5 file doesn't have data under group: %s" % (name,)
        if fn != "None":
            fh.write("the h5 file doesn't have data under group: %s\n" % (name,))
        return False
        
    #now create the munge
    munge = []
    try:
        for stride in range(num_strides):
            munge_dict = {}
            h5_path    = '/{obs:}/{name:}/stride{stride_num:02d}/'.format(obs=obs,name=name,stride_num=stride)
            things     = h5file[h5_path].items()
            for thing in things:
                #the 0th component is needed since thing is a tuple 
                #giving the (<name>, <h5 type>) 
                munge_key = str(thing[0]) 
                h5_path   = '/{obs:}/{name:}/stride{stride_num:02d}/{thing}'.format(obs=obs,name=name,stride_num=stride,thing=thing[0])
                #perform the same unpacking as was done in save_munge
                if munge_key == 'name':
                    munge_dict[munge_key] = str(h5file[h5_path][...])
                elif munge_key == 'start' or munge_key == 'stop':
                    munge_dict[munge_key] = mdates.num2date(float(h5file[h5_path][...]))
                elif munge_key == 'num_segs':
                    munge_dict[munge_key] = int(h5file[h5_path][...])
                else:                  
                    munge_dict[munge_key] = h5file[h5_path][...]
            munge.append(munge_dict)
            
    except:
        print "Can't retrieve data!"
        if fn != "None":
            fh.write("Can't retrieve data!\n")
        return False

    convert_epochs(munge)
        
    print "Done with load of %s on %s at %s" % (name,obs,dt.datetime.now())
    print "****************************************\n"
    if fn != "None":
        fh.write("Done with load of %s on %s at %s\n" % (name,obs,dt.datetime.now()))
        fh.write("****************************************\n")
        fh.close()
        
    return munge
    
#############################################################################
def convert_epochs(munge):
    """helper function that changes an numpy array of num2date's to
       datetime objects for use in plotting.
       
       Arguments:
          munge:       the munge whose epochs are to be converted
          
       Returns:
          nothing per se
       
       Example use:  
           convert_epochs(munge)
              
       Note:       
           none
       """ 
    #determine the number of strides
    num_strides = len(munge)
    
    for num in range(num_strides):
        munge[num]['epochs'] = np.array(mdates.num2date(munge[num]['epochs']))
        
#############################################################################        
def report_epochs(munge):
    
    #determine the number of the strides in the munge
    num_strides = len(munge)
    
    #report the epochs
    for i in range(num_strides):
        print i, '\t', munge[i]['epochs'][0], '  ', munge[i]['epochs'][-1]         