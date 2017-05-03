# -*- coding: utf-8 -*-
import numpy                   as np
import scipy                   as sp

num_polar = 16
num_az    = 32
###############################################################################
###############################################################################
#  DISTRIBUTION DATA
###############################################################################
###############################################################################

############################################################################### 
#
# unpack_dist_cdf takes a file handle for a distribution file
# (Fast survey or burst - specifed via type) and creates a set of  
# data structures to move the data around.  
# These structures are:
#
# dist:      dictionary with the following keys
#            Dist  - phase space distribution
#            Err   - phase space errors
#            Epoch - name says it all
#            Flag  - goodness flag  
#
# params:    dictionary with the following keys
#            Erg   - the fast survey energies
#            Phi   - instrument look direction azimuths
#            Theta - instrument look direction polar angles
#
###############################################################################
def unpack_dist_cdf(cdf_dict,obs,mode,species,ver,corrections_on,correction_override=0.0):
    dist_cdf       = cdf_dict['dist']
    dist_str       = '%s_%s_dist_%s'       % (obs,species,mode)
    dist_err_str   = '%s_%s_disterr_%s'    % (obs,species,mode)
    energy_str     = '%s_%s_energy_%s'     % (obs,species,mode)
    error_flag_str = '%s_%s_errorflags_%s' % (obs,species,mode)
    phi_str        = '%s_%s_phi_%s'        % (obs,species,mode)
    theta_str      = '%s_%s_theta_%s'      % (obs,species,mode)

    uncorrected_dist = np.asarray(dist_cdf[dist_str][:,:,:,:])
    
    if corrections_on == 0:
        corrected_dist   = uncorrected_dist
    else:
        corrected_dist   = subtract_internal_photoelectrons(cdf_dict,uncorrected_dist,obs,mode,species,correction_override)
    
    dist_data = {'Dist'   : corrected_dist,
                 'Err'    : np.asarray(dist_cdf[dist_err_str][:,:,:,:]),
                 'Epoch'  : np.asarray(dist_cdf['Epoch'][:]),
                 'Flag'   : np.asarray(dist_cdf[error_flag_str][:])}
                   
    if ver == 'ver2':
        Energy = np.asarray(dist_cdf[energy_str][:])
    if ver == 'ver3':
        Energy = np.asarray(dist_cdf[energy_str][0,:])
                  
    params   = {'Erg'    : Energy,
                'Phi'    : np.asarray(dist_cdf[phi_str][:]),
                'Theta'  : np.asarray(dist_cdf[theta_str][:])}
    
    return dist_data, params
    
    
###############################################################################    
#
#  subtract_internal_photoelectrons takes a raw phase space distribution
#  and corrects it for internally-generated photoelectons 
#
###############################################################################    
def subtract_internal_photoelectrons(cdf_dict,raw_dist,obs,mode,species,correction_override):
    #allocate space for the corrected phase-space density data structure
    corrected_dist = np.zeros(raw_dist.shape)
    
    #construct startdelphi_count_str
    start_delphi_count_str = '%s_%s_startdelphi_count_fast' % (obs,species)
    
    #pull out n_photo - the structure of the following line is understood 
    #by noting that cdf_dict['debug'] gives the moments file and 
    #.attrs['Photoelectron_model_scaling_factor'][0] get the attributes 
    #(attrs) returned as a dictionary from which the key 
    #'Photoelectron_model_scaling_factor' gives the value in some 
    #spacepy/CDF way whose 0th component is a string which is then
    #case to float (whew!!!)
    if correction_override == 0:
        n_photo = float(cdf_dict['debug'].attrs['Photoelectron_model_scaling_factor'][0])
    else:
        n_photo = correction_override
    print cdf_dict.keys(), n_photo
    
    #get the photoelectron model
    photo_mode_str = 'mms_des_bgdist_%s' % mode
    f_photo = np.asarray(cdf_dict['photo'][photo_mode_str])
    
    #loop over all times
    for k in range(len(cdf_dict['dist']['Epoch'])):
        #find start_delphi_index
        startdelphi_index       = cdf_dict['dist'][start_delphi_count_str][k]

        #construct correction_index    
        correction_index        = int(np.floor(startdelphi_index/16.0))
        
        #subtract off the correction
        corrected_dist[k,:,:,:] = raw_dist[k,:,:,:] - n_photo*f_photo[correction_index,:,:,:]
        
    #floor the negative values at zero
    corrected_dist[corrected_dist < 0.0 ] = 0.0

    #return results    
    return corrected_dist
    

###############################################################################
###############################################################################
#  Magnetic Field
###############################################################################
###############################################################################
    
###############################################################################	
#
# fetch_magnetic_field takes a file handle for an FPI debug file and creates 
# a data structure bring back the magnetic field time series.
# The structure is:
#
# B_field:   np.array Nx4 (N-number of time samples)
#            Bx - X-component of the B unit vector in index 0
#            By - Y-component of the B unit vector in index 1
#            Bz - Z-component of the B unit vector in index 2
#            B  - magnetic field magnitude         in index 3
#
###############################################################################
def fetch_magnetic_field(cdf_dict,obs,species,source="DEBUG"):
    Bx_str = '%s_%s_bentPipeB_X_DSC'       % (obs,species)
    By_str = '%s_%s_bentPipeB_Y_DSC'       % (obs,species)
    Bz_str = '%s_%s_bentPipeB_Z_DSC'       % (obs,species)
    B_str  = '%s_%s_bentPipeB_Norm'        % (obs,species)  
    
    if source == "DEBUG":
        Bx     = np.asarray(cdf_dict['debug'][Bx_str])
        By     = np.asarray(cdf_dict['debug'][By_str])
        Bz     = np.asarray(cdf_dict['debug'][Bz_str])
        B      = np.asarray(cdf_dict['debug'][B_str])
        B_field = np.vstack((Bx,By,Bz,B)).transpose()
    if source == "AFG":
        B_field = fetch_AFG_B_field(cdf_dict,species)

    return B_field

###############################################################################	
#
# fetch_AFG_B_field takes the set of CDFs (cdf_dict) and interpolated the AFG
# magnetic field to the FPI points
#
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

###############################################################################	
#
# bisect_epochs finds a set of epochs in one signal that bracket a target
# epoch.  The indices for these bracketing points are returned.
#
###############################################################################
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

###############################################################################	
#
# detect_DMPA_str returns the appropriate AFG key for processing
# obviates the need for hardcoded key 
#
###############################################################################
def detect_DMPA_str(cdf_dict):
    for k in cdf_dict['AFG'].keys():
        if 'afg_srvy_l2pre_dmpa' in k:
            return k 


###############################################################################
###############################################################################
#  Particle Computations
###############################################################################
###############################################################################

###############################################################################	
#
# compute_incoming_particle_directions takes the polar and azimuthal angles
# of the instrument look direction and returns a 2-D array (32,16) of
# the particle velocity directions
#
###############################################################################
def compute_incoming_particle_directions(mode,params):

    deg       = sp.pi/180    
    
    #in fast mode, all particle directions are the same at the same
    #time so only one 32x16 array is needed
    if mode == 'fast':
        v_dirs    = np.zeros((num_az,num_polar,3))
        for i in range(num_az):
            for j in range(num_polar):
                phi         = params['Phi']  [i]*deg
                theta       = params['Theta'][j]*deg
                v_dirs[i,j] = [-sp.sin(theta)*sp.cos(phi),\
                               -sp.sin(theta)*sp.sin(phi),\
                               -sp.cos(theta)]
    #in brst mode, the despinning can't be done completely so the 
    #phi is record varying
    if mode == 'brst':
        num_time  = params['Phi'].shape[0]
        v_dirs    = np.zeros((num_time,num_az,num_polar,3))
        for k in range(num_time):
            for i in range(num_az):
                for j in range(num_polar):
                    phi           = params['Phi'][k][i]*deg
                    theta         = params['Theta'] [j]*deg
                    v_dirs[k,i,j] = [-sp.sin(theta)*sp.cos(phi),\
                                     -sp.sin(theta)*sp.sin(phi),\
                                     -sp.cos(theta)]
                        
    return v_dirs
	
###############################################################################    
#
# compute_pitch_angles takes the velocity directions and magnetic field 
# direction at a given time and returns a (32,16) array of pitch angle per 
# pixel
#
###############################################################################
def compute_pitch_angles(mode,v_dirs,bfield,time_label):
    #Construct the B field unit vector at chosen time
    Bx = bfield[time_label,0]
    By = bfield[time_label,1]
    Bz = bfield[time_label,2]
    pitch_angles = np.zeros((num_az,num_polar))
    
    if mode == 'fast':
        v = v_dirs[:,:]
    if mode == 'brst':
        v = v_dirs[time_label,:,:]
       
    for i in range(num_az):
        for j in range(num_polar):
            U                 = v[i,j]
            dot_prod          = U[0]*Bx + U[1]*By + U[2]*Bz
            if np.abs(dot_prod) > 1.0:
                if np.abs(dot_prod) > 1.001:
                    print "dot prod. 'tween vel and B was %.3e different from 1.0" % (dot_prod - np.sign(dot_prod)*1.0)
                dot_prod = 0.999*dot_prod
            pitch_angles[i,j] = sp.arccos(dot_prod)*180/sp.pi
    

    return pitch_angles

###############################################################################	
#
# compute_counts takes the full phase space distribution and returns the counts
#
###############################################################################
def compute_counts(dist):
    counts                   = (dist['Dist']/dist['Err'])**2
    counts[np.isnan(counts)] = 0
    counts                   = np.floor(counts)
    return counts

###############################################################################
#
# compute_number_flux the full phase space distribution and the energies and
# returns the number flux
#
###############################################################################
def compute_number_flux(dist,parms):
    m_e    = 9.10938356e-31  #mass of electron in Kg
    eV     = 1.60218e-19     #energy of 1 electron-volt in Joules
    E      = parms['Erg']
    coeff  = 2.0*(E*eV/m_e)**2/E *100.0**4 #100**2 for m -> cm
    jN     = np.zeros(dist['Dist'].shape)
    counter = 0
    for c in coeff:
        jN[:,:,:,counter] = c * dist['Dist'][:,:,:,counter]
        counter           = counter + 1
        
    #jN     = coeff*FS_dist['Dist'] - don't know why this ever worked
    return jN   

    
###############################################################################
#
# compute_number_flux_error computes the errors in the number flux given the
# errors in the distribution function
#
###############################################################################
def compute_number_flux_errors(dist,parms):
    m_e    = 9.10938356e-31  #mass of electron in Kg
    eV     = 1.60218e-19     #energy of 1 electron-volt in Joules
    E      = parms['Erg']
    coeff  = 2.0*(E*eV/m_e)**2/E *100.0**4 #100**2 for m -> cm
    jN_err = np.zeros(dist['Err'].shape)
    counter = 0
    for c in coeff:
        jN_err[:,:,:,counter] = c * dist['Err'][:,:,:,counter]
        counter               = counter + 1
    
    return jN_err
    
###############################################################################	
#
# compute_omni_number_flux takes the number flux and computes the 
# averaged flux over all directions and returns the results
#
###############################################################################
def compute_omni_number_flux(jN,time_label,energy_label):
    omni_jN = np.average(jN[time_label,:,:,energy_label],axis=0)/(4.0*np.pi)
    
    return omni_jN

###############################################################################
#
# compute_limited_PAD returns a truncated PAD limited by energy and pitch 
# angle from the full data
#
###############################################################################
def compute_limited_PAD(mode,time_label,minE,maxE,minPA,maxPA,core_data):
    #find the truncated energy range
    Edata                     = core_data['parms']['Erg'][minE:maxE]
    
    #get the pitch angles and then flatten to a 1-D array
    pitch_angles              = compute_pitch_angles(mode,core_data['v_dirs'],core_data['bfield'],time_label)
    flat_PA                   = np.ndarray.flatten(pitch_angles)
    
    sub_PAD                   = np.zeros(len(range(minE,maxE)))
    counter                   = 0
    for energy_label in range(minE,maxE):
        local_jN              = np.ndarray.flatten(core_data['jN'][time_label,:,:,energy_label])
        PA_table              = np.array(zip(flat_PA,local_jN))
        PA_range              = np.where( (flat_PA > minPA) & (flat_PA < maxPA))
        sub_PAD[counter]      = np.sum(PA_table[PA_range][:,1])
        counter              += 1        
    
    sterads    = spherical_cap_area(maxPA) - spherical_cap_area(minPA)        
    return Edata, sub_PAD/sterads

###############################################################################
#
# compute_limited_PAD_error returns a truncated PAD_errpr limited by energy and 
# pitch angle from the full data
#
###############################################################################
def compute_limited_PAD_error(time_label,minE,maxE,minPA,maxPA,core_data):
    #find the truncated energy range
    Edata                     = core_data['parms']['Erg'][minE:maxE]
    
    #get the pitch angles and then flatten to a 1-D array
    pitch_angles              = compute_pitch_angles(core_data['v_dirs'],core_data['bfield'],time_label)
    flat_PA                   = np.ndarray.flatten(pitch_angles)
    
    sub_PAD_err               = np.zeros(len(range(minE,maxE)))
    counter                   = 0
    for energy_label in range(minE,maxE):
        local_jN_err          = np.ndarray.flatten(core_data['jN_err'][time_label,:,:,energy_label])
        PA_table              = np.array(zip(flat_PA,local_jN_err))
        PA_range              = np.where( (flat_PA > minPA) & (flat_PA < maxPA))
        sub_PAD_err[counter]  = np.mean(PA_table[PA_range][:,1])
        counter              += 1        
    
    sterads    = spherical_cap_area(maxPA) - spherical_cap_area(minPA)        
    return Edata, sub_PAD_err/sterads

############################################################################### 
def spherical_cap_area(theta):
    #put theta into radians
    theta_loc = np.pi*theta/180.0
    
    #compute spherical cap area based on the formula in Wikipedia
    A = 2*np.pi*(1-np.cos(theta_loc))
    
    return A        
    
###############################################################################    
#
#  load_e_data reads a distribution file and debug file and returns a host of
#  numpy arrays
#
###############################################################################
def load_particle_data(cdf_dict,obs,mode,species,ver,corrections_on,correction_override = 0,source="DEBUG"):
    
    dist_data, parms    = unpack_dist_cdf(cdf_dict,obs,mode,species,ver,corrections_on,correction_override)
    bfield              = fetch_magnetic_field(cdf_dict,obs,species,source)
    v_dirs              = compute_incoming_particle_directions(mode,parms)
    counts              = compute_counts(dist_data)
    jN                  = compute_number_flux(dist_data,parms)
    jN_err              = compute_number_flux_errors(dist_data,parms)
    
    core_data                = {}
    core_data['dist_data']   = dist_data
    core_data['parms']       = parms
    core_data['bfield']      = bfield
    core_data['v_dirs']      = v_dirs
    core_data['counts']      = counts
    core_data['jN']          = jN
    core_data['jN_err']      = jN_err
    
    return core_data    