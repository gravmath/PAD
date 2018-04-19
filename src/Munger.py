import numpy as np
from spacepy import pycdf

#configuration - time deltas
dce_delta   = 0.001
des_delta   = 0.031
dis_delta   = 0.151
fgm_delta   = 0.008
scpot_delta = 0.001

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
        elif (temp_epoch[0] - A['stop']).total_seconds() < dce_delta:
            A['stop']       = temp_epoch[-1]
            A['Egse']       = np.vstack((A['Egse'],temp_E_gse))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
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
            A['num_segs'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < fgm_delta:
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = np.vstack((A['Bgse'],temp_B_gse))            
            A['Bgsm']       = np.vstack((A['Bgsm'],temp_B_gsm))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
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
    A = make_moms_dict('%s_stride%s' % (species,counter))
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
        elif (temp_epoch[0] - A['stop']).total_seconds() < moms_delta:
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
	

