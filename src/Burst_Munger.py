import numpy as np
from spacepy import pycdf

#configuration data- directories
obs_path      = ''
fast_path     = ''
brst_path     = ''
bpsd_dir      = ''
epsd_dir      = ''
fgm_dir       = ''
fpi_emoms_dir = ''
fpi_imoms_dir = ''

#configuration data - files
fgm_pre       = ''
fgm_post      = ''
fpi_emoms_pre = ''
fpi_imoms_pre = ''
fpi_post      = ''
scpot_pre     = ''
scpot_post    = ''

#configuration - time deltas
des_delta   = 0.031
dis_delta   = 0.151
fgm_delta   = 0.008
scpot_delta = 0.001


#############################################################################
def config_directories(basedir,obs,year,month,day):
    #configuration data- directories
    obs_path      = basedir+obs+'/'
    fast_path     = year+'/'+month+'/'
    brst_path     = year+'/'+month+'/'+day+'/'
    bpsd_dir      = obs_path+'dsp/fast/l2/bpsd/'+fast_path
    epsd_dir      = obs_path+'dsp/fast/l2/epsd/'+fast_path
    fgm_dir       = obs_path+'fgm/brst/l2/'+brst_path
    fpi_emoms_dir = obs_path+'fpi/brst/l2/des-moms/'+brst_path
    fpi_imoms_dir = obs_path+'fpi/brst/l2/dis-moms/'+brst_path
    
	#configuration data - files
    fgm_pre       = '%s_fgm_brst_l2_' % (obs,)
    fgm_post      = '_v5.92.0.cdf'
    fpi_emoms_pre = '%s_fpi_brst_l2_des-moms_' % (obs,)
    fpi_imoms_pre = '%s_fpi_brst_l2_dis-moms_' % (obs,)
    fpi_post      = '_v3.2.0.cdf'
    scpot_pre     = ''
    scpot_post    = ''

    
#############################################################################
def make_Amoms(name):
    Amoms = {'name'   :name,
             'start'  :False,
             'stop'   :False,
             'bulk_vs':np.array([]),
             'epochs' :np.array([]),
             'ergs'   :np.array([]),
             'heats'  :np.array([]),
             'num_es' :np.array([]),          
             'omnis'  :np.array([]),
             'pres_s' :np.array([]),
             'temp_s' :np.array([]),
             'num_brsts':0}

    return Amoms
	
#############################################################################
def munge_moms(epoch_strings,species):
    B       = []
  
    counter = 1
    A = make_Amoms('%s_stride%s' % (species,counter))
    for e in epoch_strings:
        print '*',
        if species == 'des':
            moms       = pycdf.CDF(fpi_emoms_dir+fpi_emoms_pre+e+fpi_post)
            moms_delta = des_delta
        if species == 'dis':
            moms       = pycdf.CDF(fpi_imoms_dir+fpi_imoms_pre+e+fpi_post) 
            moms_delta = dis_delta
        temp_bulk_v = np.asarray(moms['%s_%s_bulkv_gse_brst' % (obs,species)])
        temp_epoch  = np.asarray(moms['Epoch'])
        temp_erg    = np.asarray(moms['%s_%s_energy_brst' % (obs,species)])
        temp_heat   = np.asarray(moms['%s_%s_heatq_gse_brst' % (obs,species)])
        temp_num_e  = np.asarray(moms['%s_%s_numberdensity_brst' % (obs,species)])
        temp_omni   = np.asarray(moms['%s_%s_energyspectr_omni_brst' % (obs,species)])
        temp_pres   = np.asarray(moms['%s_%s_prestensor_gse_brst' % (obs,species)])
        temp_temp   = np.asarray(moms['%s_%s_temptensor_gse_brst' % (obs,species)])
        if A['num_brsts'] == 0:
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['bulk_vs']    = temp_bulk_v
            A['epochs']     = temp_epoch
            A['ergs']       = temp_erg
            A['heats']      = temp_heat
            A['num_es']     = temp_num_e
            A['omnis']      = temp_omni
            A['pres_s']     = temp_pres
            A['temp_s']     = temp_temp            
            A['num_brsts'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < moms_delta:
            A['stop']       = temp_epoch[-1]
            A['bulk_vs']    = np.vstack((A['bulk_vs'],temp_bulk_v))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
            A['ergs']       = np.vstack((A['ergs'],temp_erg))
            A['heats']      = np.vstack((A['heats'],temp_heat))            
            A['num_es']     = np.hstack((A['num_es'],temp_num_e))
            A['omnis']      = np.vstack((A['omnis'],temp_omni))
            A['pres_s']     = np.vstack((A['pres_s'],temp_pres))
            A['temp_s']     = np.vstack((A['temp_s'],temp_temp))
            A['num_brsts'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() > moms_delta:
            B.append(A)
            counter += 1
            A               = make_Amoms('%s_stride%s' % (species,counter))
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['bulk_vs']    = temp_bulk_v
            A['epochs']     = temp_epoch
            A['ergs']       = temp_erg
            A['heats']      = temp_heat
            A['num_es']     = temp_num_e
            A['omnis']      = temp_omni
            A['pres_s']     = temp_pres
            A['temp_s']     = temp_temp            
            A['num_brsts'] += 1
    B.append(A)
    return B
	
#############################################################################
def make_Afgm(name):
    Afgm = {'name'     :name,
            'start'    :False,
            'stop'     :False,
            'epochs'   :np.array([]),
            'Bgse'     :np.array([]),
            'Bgsm'     :np.array([]),
            'num_brsts':0}

    return Afgm
	
#############################################################################
def munge_fgm(epoch_strings):
    B       = []
    counter = 1
    A       = make_Afgm('fgm_stride%s' % (counter,))
    for e in epoch_strings:
        print '*',
        fgm        = pycdf.CDF(fgm_dir+fgm_pre+e+fgm_post) 
        temp_B_gse = np.asarray(fgm['%s_fgm_b_gse_brst_l2' % (obs,)])
        temp_B_gsm = np.asarray(fgm['%s_fgm_b_gse_brst_l2' % (obs,)])
        temp_epoch  = np.asarray(fgm['Epoch'])
        if A['num_brsts'] == 0:
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = temp_B_gse
            A['Bgsm']       = temp_B_gsm
            A['epochs']     = temp_epoch
            A['num_brsts'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() < fgm_delta:
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = np.vstack((A['Bgse'],temp_B_gse))            
            A['Bgsm']       = np.vstack((A['Bgsm'],temp_B_gsm))            
            A['epochs']     = np.hstack((A['epochs'],temp_epoch))
            A['num_brsts'] += 1
        elif (temp_epoch[0] - A['stop']).total_seconds() > fgm_delta:
            B.append(A)
            counter += 1
            A               = make_Afgm('fgm_stride%s' % (counter,))
            A['start']      = temp_epoch[0]
            A['stop']       = temp_epoch[-1]
            A['Bgse']       = temp_B_gse
            A['Bgsm']       = temp_B_gsm
            A['epochs']     = temp_epoch
            A['num_brsts'] += 1
    B.append(A)
    return B