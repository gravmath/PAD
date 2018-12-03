import datetime as dt
import dateutil.tz as tz
import h5py
import matplotlib.dates as mdates
import numpy    as np
from spacepy import pycdf

#configuration - time deltas (burst)
bpsd_delta     = 2.5
epsd_delta     = 2.5
scpot_delta    = 0.001
dce_delta      = 0.00015
fgm_delta      = 0.008
fgm_delta_srvy = 1.0
des_delta      = 0.031
dis_delta      = 0.151
des_delta_fast = 5.0
dis_delta_fast = 5.0
mec_delta      = 2.5

#translations
# I believe that 'special' and 'eval' are no longer any different but the 
#naming has not been changed so as to not fix what ain't broke
#

#dsp (2) both fast
bpsd_translation = {'epochs':['Epoch','null'],
                    'freqs' :['"%s_b_freq" %(obs,)','special'],
                    'bpsd'  :['"%s_dsp_bpsd_omni_fast_l2" % (obs,)','eval']}

epsd_translation = {'epochs':['Epoch','null'],
                    'freqs' :['"%s_e_freq" %(obs,)','special'],
                    'epsd'  :['"%s_dsp_epsd_omni" % (obs,)','eval']}

#edp (2) 
dce_translation = {'epochs':['"%s_edp_epoch_brst_l2" % (obs,)','eval'],
                   'Egse'  :['"%s_edp_dce_gse_brst_l2" % (obs,)','eval']}

scpot_translation = {'epochs' : ['"%s_edp_epoch_fast_l2" % (obs,)','eval'],
                     'scpot'  : ['"%s_edp_scpot_fast_l2" % (obs,)','eval']}


#fgm (2) brst and fast
fgm_translation = {'epochs':['Epoch','null'],
                   'Bbcs'  :['"%s_fgm_b_bcs_brst_l2" % (obs,)','eval'],
                   'Bgse'  :['"%s_fgm_b_gse_brst_l2" % (obs,)','eval'],
                   'Bgsm'  :['"%s_fgm_b_gsm_brst_l2" % (obs,)','eval']}
                   
fgm_translation_srvy = {'epochs':['Epoch','null'],
                        'Bbcs'  :['"%s_fgm_b_bcs_srvy_l2" % (obs,)','eval'],
                        'Bgse'  :['"%s_fgm_b_gse_srvy_l2" % (obs,)','eval'],
                        'Bgsm'  :['"%s_fgm_b_gsm_srvy_l2" % (obs,)','eval']}                   

#fpi (8) e/i dist/moms brst/fast
edist_translation = {'epochs'    :['Epoch','null'],
                     'dist'      :['"%s_%s_dist_brst" % (obs,"des")','eval'],
                     'disterr'   :['"%s_%s_disterr_brst" % (obs,"des")','eval'],
                     'ergs'      :['"%s_%s_energy_brst" % (obs,"des")','eval'],
                     'parity'    :['"%s_%s_steptable_parity_brst" % (obs,"des")','eval'],
                     'phis'      :['"%s_%s_phi_brst" % (obs,"des")','eval'],
                     'start_dphi':['"%s_%s_startdelphi_count_brst" % (obs,"des")','eval'],
                     'thetas'    :['"%s_%s_theta_brst" % (obs,"des")','special']}
               
idist_translation = {'epochs'    :['Epoch','null'],
                     'dist'      :['"%s_%s_dist_brst" % (obs,"dis")','eval'],
                     'disterr'   :['"%s_%s_disterr_brst" % (obs,"dis")','eval'],
                     'ergs'      :['"%s_%s_energy_brst" % (obs,"dis")','eval'],
                     'phis'      :['"%s_%s_phi_brst" % (obs,"dis")','eval'],
                     'start_dphi':['"%s_%s_startdelphi_count_brst" % (obs,"dis")','eval'],
                     'thetas'    :['"%s_%s_theta_brst" % (obs,"dis")','special']}                   
               
emoms_translation = {'epochs' :['Epoch','null'],
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

            
imoms_translation = {'epochs' :['Epoch','null'],
                     'bulk_vs':['"%s_%s_bulkv_gse_brst" % (obs,"dis")','eval'],
                     'ergs'   :['"%s_%s_energy_brst" % (obs,"dis")','eval'],
                     'heats'  :['"%s_%s_heatq_gse_brst" % (obs,"dis")','eval'],
                     'num_den':['"%s_%s_numberdensity_brst" % (obs,"dis")','eval'],
                     'omnis'  :['"%s_%s_energyspectr_omni_brst" % (obs,"dis")','eval'],
                     'pres_s' :['"%s_%s_prestensor_gse_brst" % (obs,"dis")','eval'],
                     'T_par'  :['"%s_%s_temppara_brst" % (obs,"dis")','eval'],
                     'T_perp' :['"%s_%s_tempperp_brst" % (obs,"dis")','eval'],
                     'T_s'    :['"%s_%s_temptensor_gse_brst" % (obs,"dis")','eval'] }  
                     
edist_translation_fast = {'epochs'    :['Epoch','null'],
                          'dist'      :['"%s_%s_dist_fast" % (obs,"des")','eval'],
                          'disterr'   :['"%s_%s_disterr_fast" % (obs,"des")','eval'],
                          'ergs'      :['"%s_%s_energy_fast" % (obs,"des")','eval'],
                          'phis'      :['"%s_%s_phi_fast" % (obs,"des")','eval'],
                          'start_dphi':['"%s_%s_startdelphi_count_fast" % (obs,"des")','eval'],
                          'thetas'    :['"%s_%s_theta_fast" % (obs,"des")','special']}
               
idist_translation_fast = {'epochs'    :['Epoch','null'],
                          'dist'      :['"%s_%s_dist_fast" % (obs,"dis")','eval'],
                          'disterr'   :['"%s_%s_disterr_fast" % (obs,"dis")','eval'],
                          'ergs'      :['"%s_%s_energy_fast" % (obs,"dis")','eval'],
                          'phis'      :['"%s_%s_phi_fast" % (obs,"dis")','eval'],
                          'start_dphi':['"%s_%s_startdelphi_count_fast" % (obs,"dis")','eval'],
                          'thetas'    :['"%s_%s_theta_fast" % (obs,"dis")','special']}                   
               
emoms_translation_fast = {'epochs' :['Epoch','null'],
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
            
imoms_translation_fast = {'epochs' :['Epoch','null'],
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
mec_translation = {'epochs'       : ['Epoch','null'],
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
           fgm_munge = Munger.make_munge_via_translation('mms1','fgm',Munger.fgm_delta_srvy,m1f['fgm_s'],Munger.fgm_translation_srvy) 
              
       Note:  The workhorse
      """
    #import pdb; pdb.set_trace()
    if fn != "None":
        fh = open(fn,'a')
        fh.write("****************************************\n")
        fh.write("On obs %s at %s\n" % (obs,dt.datetime.now()) )        
        
        
    B       = []
    counter = 1
    segment = 1
    A       = make_data_dict_via_translation('%s_stride%s' % (type,counter),translation)


    for f in file_list:
        #open the cdf
        cdf  = pycdf.CDF(f)
        
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

        #report the time span of the first segment
        if fn == 'None':
            print 'segment %s - start: %s stop %s' % (segment,temp['epochs'][0],temp['epochs'][-1])
        else:
            fh.write("segment %s - start: %s stop %s\n" % (segment,temp['epochs'][0],temp['epochs'][-1]))
                    
        
        #pack it into the munge
        if A['num_segs'] == 0:
            if fn == 'None':
                print 'fresh segment - first stride'
            else:
               fh.write("fresh segment - first stride\n")
            A['start']      = temp['epochs'][0]
            A['stop']       = temp['epochs'][-1]
            for k in translation.keys():
                A[k] = temp[k]
            A['num_segs'] += 1
        elif (temp['epochs'][0] - A['stop']).total_seconds() < delta and (temp['epochs'][0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            if fn == 'None':
                print 'adjacency underway'
            else:
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
            if fn == 'None':
                print 'overlap detected - dealing with it'
            else:
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
            if fn == 'None':
                print 'break in adjacency - new stride'
            else:
                fh.write("break in adjacency - new stride\n")
            B.append(A)
            counter += 1
            A          = make_data_dict_via_translation('%s_stride%s' % (type,counter),translation)
            A['start'] = temp['epochs'][0]
            A['stop']  = temp['epochs'][-1]
            for k in translation.keys():
                A[k] = temp[k]    
            A['num_segs'] += 1
        segment += 1
            
    B.append(A)
    if fn == 'None':
        print 'Munged %s series for %s on %s!' % (len(B),type,obs)
    else:
        fh.write("Munged %s series for %s on %s!\n" % (len(B),type,obs))
        fh.write("Finished at %s\n" % (dt.datetime.now(),))
        fh.write("****************************************\n\n")
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
            mimic.append(make_data_dict_via_translation('dce_mimic',dce_translation))
        if 'dist' in munge[0]['name']:
            print 'dist flavored'
            mimic.append(make_data_dict_via_translation('dist_mimic',edist_translation))
        if 'fgm' in munge[0]['name']:
            print 'fgm flavored'
            mimic.append(make_data_dict_via_translation('fgm_mimic',fgm_translation))
        if 'moms' in munge[0]['name']:
            print 'moms flavored'
            mimic.append(make_data_dict_via_translation('moms_mimic',emoms_translation))
        for k in munge[j].keys():
            if type(munge[j][k]) != type(np.array([])):
                mimic[j][k] = copy.deepcopy(munge[j][k])
 
    return mimic
    
#############################################################################   
def interpolate_munge_to_munge(source_munge,target_munge):
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
          
       Returns:
           mimic_munge:  the mimicked munge with the source's data now resampled
                         to the target's epoch
       
       Example use:  
           efgm_munge = Munger.interpolate_munge_to_munge(fgm_munge,emoms_munge)
              
       Note:       
           The runner-up workhorse
       """
       
    import copy
    
    if len(source_munge) != len(target_munge):
        print 'Direct adaptation of source_munge to target_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
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
def adapt_munge_to_munge(source_munge,target_munge):
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
          
       Returns:
           mimic_munge:  the mimicked munge with the source's data now resampled
                         to the target's epoch
       
       Example use:  
           efgm_munge = Munger.adapt_munge_to_munge(fgm_munge,emoms_munge)
              
       Note:       
           The runner-up workhorse
       """
       
    import copy
    
    if len(source_munge) != len(target_munge):
        print 'Direct adaptation of source_munge to target_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
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
def save_munge(munge,name,obs,h5file):
    """Main function designed to save the munges to an H5 file for later
    retrieval and modification.
       
       Arguments:
          munge:       the munge to be saved
          name:        name of the munge in the H5 file - need not be the 
                       same as the variable name
          obs:         string with 'mms1', 'mms2', ...
          h5filename:  pointer to the h5file
          
       Returns:
           True if successful or a False and a printed error message 
           otherwise.
       
       Example use:  
           save_munge(munge,'skippy',my_h5)
              
       Note:       
           none
       """    

    #determine the number of strides
    num_strides = len(munge)
    
    for stride in range(num_strides):
        for k in munge[stride].keys():
            h5path = '/%s/%s/stride%s/%s' % (obs,name,stride,k)
            type_of_k = type(munge[stride][k])
            if type_of_k == str or type_of_k == dt.datetime or type_of_k == int:
                pass
            else:
                try:
                    if munge[stride][k].dtype == 'object' or munge[stride][k].dtype == 'O':
                        epoch_data = mdates.date2num([ munge[stride][k][i].replace(tzinfo=tz.tzutc()) for i in range(len(munge[stride][k]))])
                        #h5file.create_dataset(h5path,data=(mdates.date2num(munge[stride][k])))
                        h5file.create_dataset(h5path,data=epoch_data)
                    else:
                        h5file.create_dataset(h5path,data=(munge[stride][k]))
                except:
                    print 'Failed on ', k

    return True

#############################################################################     
def load_munge(name,obs,h5file):
    """Main function designed to load a munges from an H5 file for later
    saving and modification.
       
       Arguments:
          name:        the munge to be retrieved 
          obs:         string with 'mms1', 'mms2', ... 
          h5filename:  pointer to the h5file
          
       Returns:
           munge if successful or a False and a printed error message 
           otherwise.
       
       Example use:  
           my_munge, h5file = load_munge('skippy',my_h5)
              
       Note:       
           none
       """ 
        
    #check to make sure the name is in the file
    #by retrieving the number of strides
    try:
        groups      = h5file['/%s/%s'%(obs,name,)].items()
        num_strides = len(groups)
        print 'found %s strides' % (num_strides,)
    except:
        print "the h5 file doesn't have data under group: %s" % (name,)
        return False
        
    #now create the munge
    munge = []
    try:
        for stride in range(num_strides):
            munge_dict = {}
            things     = h5file['/%s/%s/stride%s/'%(obs,name,stride)].items()
            for thing in things:
                munge_dict[str(thing[0])] = h5file['/%s/%s/stride%s/%s'%(obs,name,stride,thing[0],)]
                #print h5file['/%s/stride%s/%s'%(name,stride,thing[0],)].shape
            munge.append(munge_dict)
    except:
        print "Can't retrieve data!"
        return False
        
    convert_epochs(munge)
    
    return munge
    
#############################################################################
def convert_epochs(munge):
    #determine the number of strides
    num_strides = len(munge)
    
    for num in range(num_strides):
        munge[num]['epochs'] = np.array(mdates.num2date(munge[num]['epochs']))
        
        
def report_epochs(munge):
    
    #determine the number of the strides in the munge
    num_strides = len(munge)
    
    #report the epochs
    for i in range(num_strides):
        print i, '\t', munge[i]['epochs'][0], '  ', munge[i]['epochs'][-1]         