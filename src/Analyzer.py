# -*- coding: utf-8 -*-
import numpy                   as np
import scipy                   as sp
import Munger

num_az         = 32
num_polar      = 16
num_ergs       = 32
num_components = 3
deg            = sp.pi/180.0

############################################################################### 
def make_record_varying(munge,key):
    """Analysis function to turn a non-record-varying field to record varying.

       Arguments:
          munge: any munge whatsoever
          key:   key in the munge that is to be converted

       Returns:
           nothing per se - adds to the munge a time-varying (record-varying) 
           list of data where once there was only a non-record varying data
       
       Example use:  
           make_record_varying(edist_brst_mms1)
              
       Note:  None
       """    
    
    #determine the number of strides
    num_strides = len(munge)
    
    for N in range(num_strides):
        #determine the number of segments
        num_segs = munge[N]['num_segs']
        
        #determine the number of tknots
        num_tknots = len(munge[N]['epochs'])
        
        #determine the number of NRV values
        num_NRV_vals = len(munge[N][key])/num_segs
        
        #copy the NRV data
        temp    = np.zeros(num_NRV_vals)
        temp[:] = munge[N][key][:num_NRV_vals]
       
        #delete the entry in the dictionary
        munge[N].pop(key)
        
        #resize the NRV data to fit the size and restore 
        #it to the munge
        munge[N][key] = np.tile(temp,(num_tknots,1))        

############################################################################### 
def calculate_incoming_particle_directions(sdist_munge):
    """Analysis function to calculate the unit vectors for each pixel that
       corresponds to the incoming particle direction.

       Arguments:
          sdist_munge: a munge of either electron or ion distribution

       Returns:
           nothing per se - adds to the munge a time-varying list of 
           unit vectors
       
       Example use:  
           calculate_incoming_particle_directions(edist_brst_mms1)
              
       Note:  works fine with either species and in brst or fast mode
       """
       
    #determine the number of strides
    num_strides = len(sdist_munge)
    
    #define spectrometer dependent constants
    
    for N in range(num_strides):
        num_tknots = len(sdist_munge[N]['epochs'])
        v_dirs = np.nan*np.zeros((num_tknots,num_az,num_polar,num_components))
        for i in range(num_az):
            for j in range(num_polar):
                v_dirs[:,i,j,0] = -sp.sin(sdist_munge[N]['thetas'][:,j]*deg)\
                                  *sp.cos(sdist_munge[N]['phis'][:,i]*deg)
                v_dirs[:,i,j,1] = -sp.sin(sdist_munge[N]['thetas'][:,j]*deg)\
                                  *sp.sin(sdist_munge[N]['phis'][:,i]*deg)
                v_dirs[:,i,j,2] = -sp.cos(sdist_munge[N]['thetas'][:,j]*deg)
        sdist_munge[N]['v_dirs'] = v_dirs

###############################################################################         
def calculate_pitch_angles(fgm_munge,sdist_munge):
    """Analysis function to calculate the pitch angles of a distribution.

       Arguments:
          fgm_munge:   a magnetic field munge
          sdist_munge: a munge of either electron or ion distribution

       Returns:
           nothing per se - adds to the sdist_munge a field 'pitch_angs' 
           for pitch angles of the distribution
       
       Example use:  
           compute_pitch_angles(fgm_munge_mms1_brst,edist_munge_mms1_brst)
              
       Note:  None
       """    
    #determine the number of strides
    num_fgm_strides   = len(fgm_munge)
    num_sdist_strides = len(sdist_munge)
    if num_fgm_strides != num_sdist_strides:
        print 'Direct comparison of fgm_munge to sdist_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
        return False
    num_strides = num_fgm_strides
    
    #adapt the fgm_munge to the sdist_munge
    sfgm_munge = Munger.adapt_munge_to_munge(fgm_munge,sdist_munge)
    
    for N in range(num_strides):
        num_tknots = len(sdist_munge[N]['epochs'])
        dot_prods  = np.zeros((num_tknots,num_az,num_polar))
        for i in range(num_az):
            for j in range(num_polar):
                dot_prods[:,i,j] = np.sum(sfgm_munge[N]['Bbcs'][:,0:3]*sdist_munge[N]['v_dirs'][:,i,j,0:3],axis=1)/\
                                   sfgm_munge[N]['Bbcs'][:,3]
        mask = np.ma.getmask(sfgm_munge[N]['Bbcs'][:,0])
        dot_prods[mask,:,:] = np.nan
        sdist_munge[N]['pitch_angs'] = np.arccos(dot_prods)*180.0/np.pi

###############################################################################         
def calculate_differential_number_flux(sdist_munge,species):
    """Analysis function to calculate the differential number flux of a 
       distribution.

       Arguments:
          sdist_munge: a munge of either electron or ion distribution
          species:     either 'electrons' or 'ions'

       Returns:
           nothing per se - adds to the sdist_munge a field 'jn' for 
           differential number flux
       
       Example use:  
           calculate_differential_number_flux(sdist_munge,species)
              
       Note:  None
       """    
    #determine number of strides
    num_strides = len(sdist_munge)
    
    #set physical constants, including the mass of the species
    c = 29979245800 #cm/s
    if species == 'electrons':
        ms = 0.5109989461e6/c**2 #ev/c^2
    if species == 'ions':
        ms = 938.2720813e6/c**2 #ev/c^2
        
    #differential number flux: j_n = 1/2 f v^4/E
    #v^2 = 2 E / m -> j_n = 2(E/m)^2/E = 2E/m^2
    for N in range(num_strides):
        jn = np.zeros(sdist_munge[N]['dist'].shape)
        f  = sdist_munge[N]['dist']
        Es = sdist_munge[N]['ergs']
        for j in range(num_az):
            for k in range(num_polar):
                for m in range(num_ergs):
                    jn[:,j,k,m] = 2.0*f[:,j,k,m]*Es[:,m]/ms**2
                    #print 2.0*Es[:,m]/ms**2
        sdist_munge[N]['jn'] = jn    
        
###############################################################################        
def calculate_omni_number_flux(sdist_munge):
    """Analysis function to calculate the omnidirectional number flux of a 
       distribution.

       Arguments:
          sdist_munge: a munge of either electron or ion distribution

       Returns:
           nothing per se - adds to the sdist_munge a field 'omni_jn' for 
           differential number flux
       
       Example use:  
           calculate_omni_number_flux(sdist_munge)
              
       Note:  None
       """    
       
    #determine the number of strides
    num_strides = len(sdist_munge)
    
    for N in range(num_strides):
        omni_jn = np.average(np.average(sdist_munge[N]['jn'][:,:,:,:],axis=1),axis=1)/(4.0*np.pi)
        sdist_munge[N]['omni_jn'] = omni_jn

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
def spherical_cap_area(theta):
    #put theta into radians
    theta_loc = np.pi*theta/180.0
    
    #compute spherical cap area based on the formula in Wikipedia
    A = 2*np.pi*(1-np.cos(theta_loc))
    
    return A        
        