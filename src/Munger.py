import numpy as np
from spacepy import pycdf

#configuration - time deltas (burst)
bpsd_delta     = 2.5
epsd_delta     = 2.5
scpot_delta    = 0.001
dce_delta      = 0.00015
fgm_delta      = 0.008
fgm_delta_fast = 1.0
des_delta      = 0.031
dis_delta      = 0.151

#translations

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
                     'phis'      :['"%s_%s_phi_brst" % (obs,"des")','eval'],
                     'start_dphi':['"%s_%s_startdelphi_count_brst" % (obs,"des")','eval'],
                     'thetas'    :['"%s_%s_theta_brst" % (obs,"des")','eval']}
               
idist_translation = {'epochs'    :['Epoch','null'],
                     'dist'      :['"%s_%s_dist_brst" % (obs,"dis")','eval'],
                     'disterr'   :['"%s_%s_disterr_brst" % (obs,"dis")','eval'],
                     'ergs'      :['"%s_%s_energy_brst" % (obs,"dis")','eval'],
                     'phis'      :['"%s_%s_phi_brst" % (obs,"dis")','eval'],
                     'start_dphi':['"%s_%s_startdelphi_count_brst" % (obs,"dis")','eval'],
                     'thetas'    :['"%s_%s_theta_brst" % (obs,"dis")','eval']}                   
               
emoms_translation = {'epochs' :['Epoch','null'],
                     'bulk_vs':['"%s_%s_bulkv_gse_brst" % (obs,"des")','eval'],
                     'ergs'   :['"%s_%s_energy_brst" % (obs,"des")','eval'],
                     'heats'  :['"%s_%s_heatq_gse_brst" % (obs,"des")','eval'],
                     'num_den':['"%s_%s_numberdensity_brst" % (obs,"des")','eval'],
                     'omnis'  :['"%s_%s_energyspectr_omni_brst" % (obs,"des")','eval'],
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
                          'thetas'    :['"%s_%s_theta_fast" % (obs,"des")','eval']}
               
idist_translation_fast = {'epochs'    :['Epoch','null'],
                          'dist'      :['"%s_%s_dist_fast" % (obs,"dis")','eval'],
                          'disterr'   :['"%s_%s_disterr_fast" % (obs,"dis")','eval'],
                          'ergs'      :['"%s_%s_energy_fast" % (obs,"dis")','eval'],
                          'phis'      :['"%s_%s_phi_fast" % (obs,"dis")','eval'],
                          'start_dphi':['"%s_%s_startdelphi_count_fast" % (obs,"dis")','eval'],
                          'thetas'    :['"%s_%s_theta_fast" % (obs,"dis")','eval']}                   
               
emoms_translation_fast = {'epochs' :['Epoch','null'],
                          'bulk_vs':['"%s_%s_bulkv_gse_fast" % (obs,"des")','eval'],
                          'ergs'   :['"%s_%s_energy_fast" % (obs,"des")','eval'],
                          'heats'  :['"%s_%s_heatq_gse_fast" % (obs,"des")','eval'],
                          'num_den':['"%s_%s_numberdensity_fast" % (obs,"des")','eval'],
                          'omnis'  :['"%s_%s_energyspectr_omni_fast" % (obs,"des")','eval'],
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
                         
#hpca (2)
#coming whenever
     
#mec (1) always fast
mec_translation = {'epochs'    : ['Epoch','null'],
                   'fieldline' : ['"%s_mec_fieldline_type" % (obs,)','eval'],
                   'losscone_n': ['"%s_mec_loss_cone_angle_n" %(obs,)','eval'],
                   'losscone_s': ['"%s_mec_loss_cone_angle_s" %(obs,)','eval'],
                   'gse_pos'   : ['"%s_mec_r_gse" % (obs,)','eval'],
                   'gsm_pos'   : ['"%s_mec_r_gsm" % (obs,)','eval'],
                   'mlt'       : ['"%s_mec_mlt" % (obs,)','eval'],                   
                   'sm_pos'    : ['"%s_mec_r_sm" % (obs,)','eval']}
                         


                         

'''                         
#############################################################################
def make_dce_dict(name):
    """Helper function for making a dce data dictionary.
       
       Arguments:
          name:  string encoding the name of the dictionary (something
                 explanatory like dce_brst_mms1)

       Returns:
           Adce: data dictionary with standard keys:  'name', 'start', 'stop', 
                 'num_segs', 'epochs', as well as 'Egse' (low-frequency electric 
                 field in GSE)
       
       Example use:  
           my_dce_dict =  make_dce_dict('dce_brst_mms1')
              
       Note:  N/A
    """

    Adce = {'name'     :name,
            'start'    :False,
            'stop'     :False,
            'epochs'   :np.array([]),
            'Egse'     :np.array([]),
            'num_segs':0}

    return Adce

#############################################################################
def make_dist_dict(name):
    """Helper function for making a dist data dictionary.
       
       Arguments:
          name:  string encoding the name of the dictionary (something
                 explanatory like edist_brst_mms1)

       Returns:
           Adce: data dictionary with standard keys:  'name', 'start', 'stop', 
                 'num_segs', 'epochs', as well as 'dist', 'disterr', 'ergs',
                 'phis', 'startdelphi_count', and 'thetas'
       
       Example use:  
           my_edist_dict =  make_dist_dict('edist_brst_mms1')
              
       Note:  N/A
    """

    Adist = {'name'              :name,
             'start'             :False,
             'stop'              :False,
             'epochs'            :np.array([]),
             'dist'              :np.array([]),
             'disterr'           :np.array([]),
             'ergs'              :np.array([]),
             'phis'              :np.array([]),
             'startdelphi_count' :np.array([]),
             'thetas'            :np.array([]),
             'num_segs':0}

    return Adist    
    
#############################################################################
def make_fgm_dict(name):
    """Helper function for making a fgm data dictionary.
       
       Arguments:
          name:  string encoding the name of the dictionary (something
                 explanatory like fgm_brst_mms1)

       Returns:
           Afgm: data dictionary with standard keys:  'name', 'start', 'stop', 
                 'num_segs', 'epochs', as well as 'Bgse' and 'Bgsm' (low-
                 frequency magnetic field in GSE and GSM)
       
       Example use:  
           my_fgm_dict =  make_fgm_dict('fgm_fast_mms1')
              
       Note:  N/A
    """
    Afgm = {'name'     :name,
            'start'    :False,
            'stop'     :False,
            'epochs'   :np.array([]),
            'Bgse'     :np.array([]),
            'Bgsm'     :np.array([]),
            'num_segs':0}

    return Afgm

############################################################################
def make_moms_dict(name):
    """Helper function for making a moms data dictionary.
       
       Arguments:
          name:  string encoding the name of the dictionary (something
                 explanatory like emoms_brst_mms1)

       Returns:
           Amoms: data dictionary with standard keys:  'name', 'start', 'stop', 
                  'num_segs', 'epochs', as well as 'ergs', 'heats', 'num_den',
                  'omnis', 'pres_s', 'T_par', 'T_perp', and 'T_s' (ESA energies,
                  heat flux, number density, omni-directional spectrogram, 
                  pressure tensor, parallel temperature, perpendicular temperature,
                  and the temperature tensor - all rank-1 & greater tensors in GSE)
       
       Example use:  
           my_emoms_dict =  make_moms_dict('emoms_brst_mms1')
              
       Note:  The user specifies the species (eletrons or ion - DES or DIS), when 
              the data are placed into the dictionary
    """

    Amoms = {'name'      :name,
             'start'     :False,
             'stop'      :False,
             'bulk_vs'   :np.array([]),
             'epochs'    :np.array([]),
             'ergs'      :np.array([]),
             'heats'     :np.array([]),
             'num_den'   :np.array([]),          
             'omnis'     :np.array([]),
             'pres_s'    :np.array([]),
             'T_par'     :np.array([]),
             'T_perp'    :np.array([]),
             'T_s'       :np.array([]),
             'num_segs' :0}

    return Amoms

#############################################################################
def munge_dce(file_list,obs):
    """Core function for making a dce munged data dictionary.
       
       Arguments:
          file_list:  list of files to be munged
          obs:        'mms1', 'mms2', 'mms3', 'mms4'

       Returns:  munged dce data dictionary with standard keys:  'name', 
                 'start', 'stop', 'num_segs', 'epochs', as well as 'Egse' 
                 (low-frequency electric field in GSE)
                 
       Example use:  
           dce_brst_mms1 =  munge_dce(my_file_list)
              
       Note:  N/A
    """

    B       = []
    counter = 1
    A       = make_dce_dict('dce_stride%s' % (counter,))
    for file in file_list:
        print '*',
        dce        = pycdf.CDF(file) 
        temp_E_gse = np.asarray(dce['%s_edp_dce_gse_brst_l2' % (obs,)])
        temp_epoch = np.asarray(dce['%s_edp_epoch_brst_l2' % (obs,)])
        if A['num_segs'] == 0:
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['Egse']       = temp_E_gse
            A['epochs']     = temp_epoch
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < dce_delta and (temp_epoch[0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            A['stop']       = temp_epoch[-1]
            A['Egse']       = np.vstack((A['Egse'],temp_E_gse))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < dce_delta and (temp_epoch[0] - A['stop']).total_seconds() < 0.0:
            #overlaps exist - assume that the coincident measurements are identical
            print 'O',
            time_delta      = np.array([(temp_epoch[j] - A['stop']).total_seconds() for j in range(len(temp_epoch))])
            new_points      = np.where(time_delta > dce_delta)
            if temp_epoch[-1] > A['stop']:
                A['stop']   = temp_epoch[-1]
            A['epochs']     = np.hstack((A['epochs'],temp_epoch[new_points]))
            A['Egse']       = np.vstack((A['Egse'],temp_E_gse[new_points]))
            A['num_segs']  += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() > dce_delta:
            B.append(A)
            counter += 1
            A               = make_dce_dict('dce_stride%s' % (counter,))
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['Egse']       = temp_E_gse
            A['epochs']     = temp_epoch
            A['num_segs'] += 1
    B.append(A)
    return B
 
#############################################################################
def munge_dist(file_list,obs,species):
    """Core function for making a dist munged data dictionary.
       
       Arguments:
          file_list:  list of files to be munged
          obs:        'mms1', 'mms2', 'mms3', 'mms4'          
          species:    'des' or 'dis' for electrons or ions

       Returns:  munged moms data dictionary with standard keys: 'name', 'start', 
                 'stop', 'num_segs', 'epochs', as well as 'dist', 'disterr', 'ergs',
                 'phis', 'startdelphi_count', and 'thetas'
       
       Example use:  
           edist_brst_mms1 =  munge_dist(my_file_list,'des')
           
       Note:  The user specifies the species (eletrons or ion - DES or DIS), when 
              the data are placed into the dictionary                 
    """

    B       = []
  
    counter = 1
    A = make_dist_dict('%s_dist_stride%s' % (species,counter))
    for file in file_list:
        print '*',
        dist = pycdf.CDF(file)
        if species == 'des':
            dist_delta = des_delta
        if species == 'dis':
            dist_delta = dis_delta
        temp_epoch             = np.asarray(dist['Epoch'])
        temp_dist              = np.asarray(dist['%s_%s_dist_brst' % (obs,species)])
        temp_dist_err          = np.asarray(dist['%s_%s_disterr_brst' % (obs,species)])
        temp_erg               = np.asarray(dist['%s_%s_energy_brst' % (obs,species)])
        temp_phi               = np.asarray(dist['%s_%s_phi_brst' % (obs,species)])
        temp_startdelphi_count = np.asarray(dist['%s_%s_startdelphi_count_brst' % (obs,species)])
        temp_theta             = np.asarray(dist['%s_%s_theta_brst' % (obs,species)])
        if A['num_segs'] == 0:
            A['start']             = temp_epoch[0]
            A['stop']              = temp_epoch[-1]
            A['epochs']            = temp_epoch
            A['dist']              = temp_dist
            A['disterr']           = temp_dist_err
            A['ergs']              = temp_erg
            A['phis']              = temp_phi
            A['startdelphi_count'] = temp_startdelphi_count
            A['thetas']            = temp_theta
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < dist_delta and (temp_epoch[0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            A['stop']              = temp_epoch[-1]
            A['epochs']            = np.hstack((A['epochs'],temp_epoch))
            A['dist']              = np.vstack((A['dist'],temp_dist))
            A['disterr']           = np.vstack((A['disterr'],temp_dist_err))                        
            A['ergs']              = np.vstack((A['ergs'],temp_erg))
            A['phis']              = np.vstack((A['phis'],temp_phi))            
            A['startdelphi_count'] = np.hstack((A['startdelphi_count'],temp_startdelphi_count))
            A['thetas']            = np.vstack((A['thetas'],temp_omni))
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < dist_delta and (temp_epoch[0] - A['stop']).total_seconds() < 0.0:
            #overlaps exist - assume that the coincident measurements are identical
            print 'O',
            time_delta      = np.array([(temp_epoch[j] - A['stop']).total_seconds() for j in range(len(temp_epoch))])
            new_points      = np.where(time_delta > dce_delta)
            if temp_epoch[-1] > A['stop']:
                A['stop']   = temp_epoch[-1]
            A['epochs']            = np.hstack((A['epochs'],temp_epoch))
            A['dist']              = np.vstack((A['dist'],temp_dist))
            A['disterr']           = np.vstack((A['disterr'],temp_dist_err))                        
            A['ergs']              = np.vstack((A['ergs'],temp_erg))
            A['phis']              = np.vstack((A['phis'],temp_phi))            
            A['startdelphi_count'] = np.hstack((A['startdelphi_count'],temp_startdelphi_count))
            A['thetas']            = np.vstack((A['thetas'],temp_omni))
            A['num_segs']  += 1            
        elif (temp_epoch[0] - A['stop']).total_seconds() > dist_delta:
            B.append(A)
            counter += 1
            A               = make_dist_dict('%s_dist_stride%s' % (species,counter))
            A['start']             = temp_epoch[0]
            A['stop']              = temp_epoch[-1]
            A['epochs']            = temp_epoch
            A['dist']              = temp_dist
            A['disterr']           = temp_dist_err
            A['ergs']              = temp_erg
            A['phis']              = temp_phi
            A['startdelphi_count'] = temp_startdelphi_count
            A['thetas']            = temp_theta
            A['num_segs'] += 1
    B.append(A)
    return B
 
    
#############################################################################
def munge_fgm(file_list,obs):
    """Core function for making a fgm munged data dictionary.
       
       Arguments:
          file_list:  list of files to be munged
          obs:        'mms1', 'mms2', 'mms3', 'mms4'          

       Returns:  munged fgm data dictionary with standard keys:  'name', 
                 'start', 'stop', 'num_segs', 'epochs', as well as 'Bgse' 
                 and 'Bgsm' (low-frequency magnetic field in GSE and GSM)
                 
       Example use:  
           fgm_brst_mms1 =  munge_fgm(my_file_list)
              
       Note:  N/A
    """

    B       = []
    counter = 1
    A       = make_fgm_dict('fgm_stride%s' % (counter,))
    for file in file_list:
        print '*',
        fgm        = pycdf.CDF(file) 
        temp_B_gse = np.asarray(fgm['%s_fgm_b_gse_brst_l2' % (obs,)])
        temp_B_gsm = np.asarray(fgm['%s_fgm_b_gse_brst_l2' % (obs,)])
        temp_epoch  = np.asarray(fgm['Epoch'])
        if A['num_segs'] == 0:
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = temp_B_gse
            A['Bgsm']       = temp_B_gsm
            A['epochs']     = temp_epoch
            A['num_segs']  += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < fgm_delta and (temp_epoch[0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = np.vstack((A['Bgse'],temp_B_gse))            
            A['Bgsm']       = np.vstack((A['Bgsm'],temp_B_gsm))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
            A['num_segs']  += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < fgm_delta and (temp_epoch[0] - A['stop']).total_seconds() < 0.0:
            #overlaps exists - assume that the coincident measurements are identical
            print 'O',
            time_delta      = np.array([(temp_epoch[j] - A['stop']).total_seconds() for j in range(len(temp_epoch))])
            new_points      = np.where(time_delta > fgm_delta)
            if temp_epoch[-1] > A['stop']:
                A['stop']   = temp_epoch[-1]
            A['epochs']     = np.hstack((A['epochs'],temp_epoch[new_points]))
            A['Bgse']       = np.vstack((A['Bgse'],temp_B_gse[new_points]))
            A['Bgsm']       = np.vstack((A['Bgsm'],temp_B_gsm[new_points]))
            A['num_segs']  += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() > fgm_delta:
            B.append(A)
            counter += 1
            A               = make_fgm_dict('fgm_stride%s' % (counter,))
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = temp_B_gse
            A['Bgsm']       = temp_B_gsm
            A['epochs']     = temp_epoch
            A['num_segs'] += 1
    B.append(A)
    return B

	
#############################################################################
def munge_moms(file_list,obs,species):
    """Core function for making a moms munged data dictionary.
       
       Arguments:
          file_list:  list of files to be munged
          obs:        'mms1', 'mms2', 'mms3', 'mms4'          
          species:    'des' or 'dis' for electrons or ions

       Returns:  munged moms data dictionary with standard keys: 'name', 'start', 
                 'stop', 'num_segs', 'epochs', as well as 'ergs', 'heats', 'num_den',
                  'omnis', 'pres_s', 'T_par', 'T_perp', and 'T_s' (ESA energies,
                  heat flux, number density, omni-directional spectrogram, 
                  pressure tensor, parallel temperature, perpendicular temperature,
                  and the temperature tensor - all rank-1 & greater tensors in GSE)
       
       Example use:  
           emoms_brst_mms1 =  munge_moms(my_file_list,'des')
           
       Note:  The user specifies the species (eletrons or ion - DES or DIS), when 
              the data are placed into the dictionary                 
    """

    B       = []
  
    counter = 1
    A = make_moms_dict('%s_moms_stride%s' % (species,counter))
    for file in file_list:
        print '*',
        moms = pycdf.CDF(file)
        if species == 'des':
            moms_delta = des_delta
        if species == 'dis':
            moms_delta = dis_delta
        temp_bulk_v = np.asarray(moms['%s_%s_bulkv_gse_brst' % (obs,species)])
        temp_epoch  = np.asarray(moms['Epoch'])
        temp_erg    = np.asarray(moms['%s_%s_energy_brst' % (obs,species)])
        temp_heat   = np.asarray(moms['%s_%s_heatq_gse_brst' % (obs,species)])
        temp_num_e  = np.asarray(moms['%s_%s_numberdensity_brst' % (obs,species)])
        temp_omni   = np.asarray(moms['%s_%s_energyspectr_omni_brst' % (obs,species)])
        temp_pres   = np.asarray(moms['%s_%s_prestensor_gse_brst' % (obs,species)])
        temp_T_par  = np.asarray(moms['%s_%s_temppara_brst' % (obs,species)])
        temp_T_perp = np.asarray(moms['%s_%s_tempperp_brst' % (obs,species)])
        temp_T      = np.asarray(moms['%s_%s_temptensor_gse_brst' % (obs,species)])
        if A['num_segs'] == 0:
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['bulk_vs']    = temp_bulk_v
            A['epochs']     = temp_epoch
            A['ergs']       = temp_erg
            A['heats']      = temp_heat
            A['num_den']    = temp_num_e
            A['omnis']      = temp_omni
            A['pres_s']     = temp_pres
            A['T_par']      = temp_T_par
            A['T_perp']     = temp_T_perp
            A['T_s']        = temp_T            
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < moms_delta & (temp_epoch[0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            A['stop']       = temp_epoch[-1]
            A['bulk_vs']    = np.vstack((A['bulk_vs'],temp_bulk_v))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
            A['ergs']       = np.vstack((A['ergs'],temp_erg))
            A['heats']      = np.vstack((A['heats'],temp_heat))            
            A['num_den']    = np.hstack((A['num_den'],temp_num_e))
            A['omnis']      = np.vstack((A['omnis'],temp_omni))
            A['pres_s']     = np.vstack((A['pres_s'],temp_pres))
            A['T_par']      = np.hstack((A['T_par'],temp_T_par))
            A['T_perp']     = np.hstack((A['T_perp'],temp_T_perp))
            A['T_s']        = np.vstack((A['T_s'],temp_T))
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < moms_delta and (temp_epoch[0] - A['stop']).total_seconds() < 0.0:
            #overlaps exist - assume that the coincident measurements are identical
            print 'O',
            time_delta      = np.array([(temp_epoch[j] - A['stop']).total_seconds() for j in range(len(temp_epoch))])
            new_points      = np.where(time_delta > dce_delta)
            if temp_epoch[-1] > A['stop']:
                A['stop']   = temp_epoch[-1]
            A['epochs']     = np.hstack((A['epochs'],temp_epoch[new_points]))
            A['bulk_vs']    = np.vstack((A['bulk_vs'],temp_bulk_v[new_points]))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch[new_points]))
            A['ergs']       = np.vstack((A['ergs'],temp_erg[new_points]))
            A['heats']      = np.vstack((A['heats'],temp_heat[new_points]))            
            A['num_den']    = np.hstack((A['num_den'],temp_num_e[new_points]))
            A['omnis']      = np.vstack((A['omnis'],temp_omni[new_points]))
            A['pres_s']     = np.vstack((A['pres_s'],temp_pres[new_points]))
            A['T_par']      = np.hstack((A['T_par'],temp_T_par[new_points]))
            A['T_perp']     = np.hstack((A['T_perp'],temp_T_perp[new_points]))
            A['T_s']        = np.vstack((A['T_s'],temp_T[new_points]))
            A['num_segs']  += 1            
        elif (temp_epoch[0] - A['stop']).total_seconds() > moms_delta:
            B.append(A)
            counter += 1
            A               = make_moms_dict('%s_stride%s' % (species,counter))
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['bulk_vs']    = temp_bulk_v
            A['epochs']     = temp_epoch
            A['ergs']       = temp_erg
            A['heats']      = temp_heat
            A['num_den']    = temp_num_e
            A['omnis']      = temp_omni
            A['pres_s']     = temp_pres
            A['T_par']      = temp_T_par
            A['T_perp']     = temp_T_perp
            A['T_s']        = temp_T            
            A['num_segs'] += 1
    B.append(A)
    for i in range(len(B)):
        B[i]['trace_p'] = np.array([np.trace(B[i]['pres_s'][j]) for j in range(len(B[i]['epochs']))])
    return B
    
###############################################################################
def fetch_AFG_B_field(cdf_dict,species):
    B       = np.zeros( (len(cdf_dict['dist']['Epoch']),4) )
    B_str   = detect_DMPA_str(cdf_dict)

    counter = 0
    for e in cdf_dict['dist']['Epoch']:
        low_index, high_index = bisect_epochs(e,cdf_dict['AFG']['Epoch'])
        epoch_low    = cdf_dict['AFG']['Epoch'][low_index]
        epoch_high   = cdf_dict['AFG']['Epoch'][high_index]        
        h            = (epoch_high - epoch_low).total_seconds()
        dt           = (epoch_high - e).total_seconds()
        B_low        = cdf_dict['AFG'][B_str][low_index]
        B_high       = cdf_dict['AFG'][B_str][high_index]
        B_e          = B_low + (B_high-B_low)*(dt/h)
        B[counter,:] = B_e
        counter      = counter + 1
        
    B[:,3] = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)
    B[:,0] = B[:,0]/B[:,3]
    B[:,1] = B[:,1]/B[:,3]
    B[:,2] = B[:,2]/B[:,3]

    return B
    
#############################################################################
def bisect_epochs(target_epoch,epochs):
    low_index   = 0
    high_index  = len(epochs)

    old_trial_index = low_index
    trial_index     = high_index
    
    while( abs(old_trial_index - trial_index) > 1):
        old_trial_index = trial_index
        trial_index     = np.int(np.floor((high_index + low_index)/2.0 ))
        trial_epoch     = epochs[trial_index]
        if (trial_epoch - target_epoch).total_seconds() < 0:
            low_index       = trial_index
        else:
            high_index = trial_index

    #a test for what may happen for odd number of points in the 
    #epoch array
    if(high_index == low_index):
        #found out if the common epoch is above or below target
        if (epochs[high_index] - target_epoch).total_seconds() < 0:
            high_index = high_index + 1
        else:
            low_index  = low_index - 1

    return low_index, high_index
    
'''

###############################################################################	
#
#
# New function approach
# 
#
###############################################################################    
    
#############################################################################     
def make_data_dict_via_translation(name,translation):
    """Generizc helper function for creating a custom data dictionary 
       by using a translation"""
       
    #make baseline structure for the munge
    A = {'name':name,
         'start':False,
         'stop':False,
         'num_segs':0}
    
    for k in translation.keys():
        A[k] = np.asarray([])
        
    return A    
    
#############################################################################
def make_munge_via_translation(obs,type,delta,file_list,translation):
    """Core function that uses an instrument-tailored hash to make a generic
      munge of all the data fed it in file_list"""

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
                temp[k] = np.asarray(cdf[translation[k][0]])
            if translation[k][1] == 'eval':
                temp[k] = np.asarray(cdf[eval(translation[k][0])])
            if translation[k][1] == 'special':
                temp[k] = np.asarray(cdf[eval(translation[k][0])][:])
        
        #close the cdf
        cdf.close()

        #report the time span of the first segment
        print 'segment %s - start: %s stop %s' % (segment,temp['epochs'][0],temp['epochs'][-1])
        
        #pack it into the munge
        if A['num_segs'] == 0:
            print 'fresh segment - first stride'
            A['start']      = temp['epochs'][0]
            A['stop']       = temp['epochs'][-1]
            for k in translation.keys():
                A[k] = temp[k]
            A['num_segs'] += 1
        elif (temp['epochs'][0] - A['stop']).total_seconds() < delta and (temp['epochs'][0] - A['stop']).total_seconds() > 0.0:
            #segments are close enough to be considered adjacent but they don't overlap
            print 'adjacency underway'
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
    print 'Munged %s series for %s on %s!' % (len(B),type,obs)
    return B
    
#############################################################################        
def interpolate_to_epoch(source_epoch,source_data,target_epoch):
    """Helper function that linearly interpolates a single time series 
      (source_epoch and source_data) to the time knots of another series
      (target_epoch)
      
      Assumes that the time-knots are python datetime objects so that it
      can use matplotlib's date2num function
      
       No error checking is performed"""
       
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
       structure as an existing munge but lacking data"""
       
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
def adapt_munge_to_munge(source_munge,target_munge):
    """Core function designed to take two munges of the same form (i.e
       several burst segments from a given time period) and to linearly
       interpolate the data from the source one to the time-knots of the
       target.  
       
       In doing so, a mimicked source_munge is returned with time-knots 
       consistent with target_munge.
       
       The use case is to adapt fgm or dce data to des or dis time-knots
       for JdotE or curlometery or the like."""
       
    import copy
    
    if len(source_munge) != len(target_munge):
        print 'Direct adaptation of source_munge to target_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
        
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
                data_shape[0]     = num_target_times
                data_shape        = tuple(data_shape)
                num_dim           = len(component_shape)
                mimic_munge[j][k] = np.zeros(data_shape)                
                if num_dim == 1: #scalar or vector series
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