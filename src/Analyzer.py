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
    """Analysis function to compute the pitch angles of a distribution.

       Arguments:
          fgm_munge:   a magnetic field munge
          sdist_munge: a munge of either electron or ion distribution

       Returns:
           nothing per se - adds to the sdist_munge a field for pitch angles
       
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
    #determine number of strides
    num_strides = len(sdist_munge)
    
    #set physical constants, including the mass of the species
    c = 29979245800 #cm/s
    if species == 'emoms':
        ms = 0.5109989461e6/c**2 #ev/c^2
    if species == 'imoms':
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
                    print 2.0*Es[:,m]/ms**2
        sdist_munge[N]['jn'] = jn    
        
        
def compute_number_flux(sdist_munge):
    m_e    = 9.10938356e-31  #mass of electron in Kg
    eV     = 1.60218e-19     #energy of 1 electron-volt in Joules
    N      = 0
    E      = sdist_munge[N]['ergs'][0,:]
    coeff  = 2.0*(E*eV/m_e)**2/E *100.0**4 #100**2 for m -> cm
    jN     = np.zeros(sdist_munge[N]['dist'].shape)
    counter = 0
    for c in coeff:
        jN[:,:,:,counter] = c * sdist_munge[N]['dist'][:,:,:,counter]
        counter           = counter + 1
        
    #jN     = coeff*FS_dist['Dist'] - don't know why this ever worked
    return jN          