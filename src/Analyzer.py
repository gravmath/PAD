# -*- coding: utf-8 -*-
import matplotlib.dates        as mdates
import numpy                   as np
import scipy                   as sp
import scipy.interpolate       as interp
import Munger

num_az         = 32
num_polar      = 16
num_ergs       = 32
num_components = 3

num_pixels     = num_az*num_polar
num_channels   = num_az*num_polar*num_ergs
deg            = sp.pi/180.0

############################################################################### 
def calculate_mec_to_tknot(tknot,mec_munge,frame):
    """Helper function to interpolate MEC data to a given tknot.
       The primary use is to decorate the tick labels on a matplotlib plot

       Arguments:
          tknot:     specific time-knot
          mec_munge: munge of MEC data - usually ID
          frame:     'gse' or 'gsm'

       Returns:
          [xknot, yknot, zknot, mltknot] - 4-D array with the spececraft 
                                           position and magnetic local time. 
                                           
       Example use:  
           [x,y,z,M] = calc_mec_to_tknot(tknot,mec_munge,frame)
              
       Note:  None
       """    

    if frame == 'gse':
        mec_interp_x = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['gse_pos'][:,0])
        mec_interp_y = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['gse_pos'][:,1])
        mec_interp_z = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['gse_pos'][:,2])
    if frame == 'gsm':
        mec_interp_x = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['gsm_pos'][:,0])
        mec_interp_y = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['gsm_pos'][:,1])
        mec_interp_z = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['gsm_pos'][:,2])
    mec_interp_mlt = interp.interp1d(mdates.date2num(mec_munge[0]['epochs']),mec_munge[0]['mlt'][:])
        
    Re = 6378.14
    
    xknot   = float(mec_interp_x(tknot)/Re)
    yknot   = float(mec_interp_y(tknot)/Re)
    zknot   = float(mec_interp_z(tknot)/Re)
    mltknot = float(mec_interp_mlt(tknot))
    
    return [xknot,yknot,zknot,mltknot]

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
def calculate_index(data,data_bin):
    """Analysis function to calculate the indices that map a 1D array of values
       to a specfic bin.

       Arguments:
          data:       1D array of raw values to be mapped to a bin
          data_bin:   1D array with the bins

       Returns:
           index:     1D array of indices stating were the data belongs in 
                      the bin
           flag:      1D array of 0's or 1's indicating where the data belongs
                      in the bin (1) or falls outside (0)
       
       Example use:  
           index, flag = calculate_index(vals,my_bin)
              
       Note:  Assumes a 1D array
       """    
    delta   = data_bin[1] - data_bin[0]
    num_pts = len(data_bin) 
    
    index = np.around((data-data_bin[0])/delta).astype(int)
    flag  = np.ones(len(data))
    
    subzero_indices          = np.where(index < 0)
    over_9000_indices        = np.where(index > num_pts - 1)
    index[subzero_indices]   = 0
    flag[subzero_indices]    = 0.0
    index[over_9000_indices] = num_pts - 1
    flag[over_9000_indices]  = 0.0
    
    return index, flag  

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
def calculate_incoming_particle_velocities(sdist_munge,species):
    """Analysis function to calculate the particle velocities for each pixel 
       at each time step and add those data (time, num_channels,vector) to
       the munge.

       Arguments:
          sdist_munge: a munge of either electron or ion distribution
          species:     string of either 'electrons' or 'ions'

       Returns:
           nothing per se - adds to the munge a time-varying list of 
           particle velocity vectors
       
       Example use:  
           calculate_incoming_particle_velocities(edist_brst_mms1,'electrons')
              
       Note:  works fine with either species and in brst or fast mode
       """
    #determine the number of strides
    num_stride = len(sdist_munge)

    #set physical constants, including the mass of the species
    c    = 299792.458000 #km/s
    if species == 'electrons':
        ms = 0.5109989461e6/c**2 #ev/c^2
    if species == 'ions':
        ms = 938.2720813e6/c**2 #ev/c^2

    for N in range(num_stride):
        #determine the number of tknots
        num_tknots = len(sdist_munge[N]['epochs'])
        
        #tile and reshape the incoming velocity directions
        v_dirs     = np.tile(sdist_munge[N]['v_dirs'],num_ergs)
        v_dirs     = v_dirs.reshape(num_tknots,num_channels,num_components)
        
        #tile and reshape the energies
        ergs       = np.tile(sdist_munge[N]['ergs'],num_pixels)
        ergs       = ergs.reshape(num_tknots,num_channels)
        
        #compute the speeds and velocities
        speeds      = np.sqrt(2.0*ergs/ms)
        vels        = np.zeros((num_tknots,num_channels,num_components))
        vels[:,:,0] = speeds*v_dirs[:,:,0]
        vels[:,:,1] = speeds*v_dirs[:,:,1]
        vels[:,:,2] = speeds*v_dirs[:,:,2]        
        
        #add to the munge
        sdist_munge[N]['vels']       = vels
        
############################################################################### 
def calculate_B_local_transformation(sfgm_munge,smoms_munge):
    """Analysis function to calculate the transformation matrices from GSE to
       local magnetic coordinates (para, perp1, perp2) at each time step and 
       add those data to the fgm_munge.

       Arguments:
          sfgm_munge:  fgm_munge adapted to the smoms_munge
          smoms_munge: a munge of either electron or ion moments

       Returns:
           nothing per se - adds to the fgm_munge a time-varying list of 
           transformation matrices
       
       Example use:  
           calculate_B_local_transformation(efgm_munge,emoms_munge)
              
       Note:  works fine with either species and in brst or fast mode
       """
    #determine the number of strides
    num_strides = len(smoms_munge)
    
   
    for N in range(num_strides):
        #determine the number of tknots in the current time series
        num_tknots = len(sfgm_munge[N]['epochs'])
        Bgse       = sfgm_munge[N]['Bgse'][:,0:3]
        U          = smoms_munge[N]['bulk_vs']
        para       = Bgse/np.linalg.norm(Bgse,axis=1).reshape((num_tknots,1))
        perp2      = np.cross(Bgse,U)
        perp2      = perp2/np.linalg.norm(perp2,axis=1).reshape((num_tknots,1))
        perp1      = np.cross(perp2,para)
        T          = np.hstack((para,perp1,perp2)).reshape(num_tknots,3,3)      
        sfgm_munge[N]['T_localB'] = T        

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
def calculate_energy_PAD_at_time_and_energy(sdist_munge,timelabel,erg_index,num_angs):
    """Analysis function to calculate the pitch angle distribution
    

       Arguments:
          sdist_munge: a munge of either electron or ion distribution
          timelabel:   integer specifying the array index in the 'epochs' array
                       for the time in question
          erg_index:   integer specifying the array index in the 'ergs' array
                       for the energy in question
          num_angs:    integer specifying the number of angular bins for the PAD

       Returns:
          ang_bins:    angular bins for the PAD
          PAD:         the desired pitch angle distribution           
          
       Example use:  
           ang_bins, PAD = calculate_energy_PAD_at_time_and_energy(sdist_munge,timelabel,erg_index,num_angs)
              
       Note:  None
       """ 
    sliced_jn = np.ndarray.flatten(sdist_munge[0]['jn'][timelabel,:,:,erg_index])
    sliced_pa = np.ndarray.flatten(sdist_munge[0]['pitch_angs'][timelabel,:,:])

    PAD      = np.zeros(num_angs-1)
    angs     = np.linspace(0,180,num_angs)
    ang_bins = np.array([(angs[j]+angs[j-1])/2.0 for j in range(1,num_angs)])
    for N in range(1,num_angs):
        valid_pts  = np.logical_and(sliced_pa>=angs[N-1],sliced_pa<angs[N])
        PAD[N-1]   = np.average(sliced_jn[valid_pts])
    return ang_bins, PAD
        
        
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
def calculate_scalar_pressure(smoms_munge):
    """Analysis function to calculate the scalar pressure of a particle moments
       munge.

       Arguments:
          smoms_munge: a munge of either electron or ion moments

       Returns:
           nothing per se - adds to the sdist_munge a field 'scalar_p' for 
           differential number flux
       
       Example use:  
           calculate_scalar_pressure(smoms_munge)
              
       Note:  None
       """    
    #determine the number of strides
    num_strides = len(smoms_munge)
    
    for N in range(num_strides):
        scalar_p    = np.zeros(len(smoms_munge[N]['epochs']))
        scalar_p[:] = (smoms_munge[N]['pres_s'][:,0,0]+smoms_munge[N]['pres_s'][:,1,1]+smoms_munge[N]['pres_s'][:,2,2])/3.0
        
        smoms_munge[N]['scalar_p'] = scalar_p

###############################################################################          
def calculate_plasma_beta(fgm_munge,smoms_munge):
    """Analysis function to calculate the plasma beta for a particle moments
       munge.

       Arguments:
          fgm_munge:   a munge for fgm (unadapted)
          smoms_munge: a munge of either electron or ion moments

       Returns:
           nothing per se - adds to the sdist_munge a field 'beta' for 
           differential number flux
       
       Example use:  
           calculate_plasma_beta(fgm_munge,emoms_munge)
              
       Note:  None
       """ 
    #determine the number of strides
    num_strides = len(smoms_munge)
    
    #make a compatible fgm_munge
    sfgm_munge = Munger.adapt_munge_to_munge(fgm_munge,smoms_munge)
    
    #enter some physical constants
    mu_0 = 4*np.pi*1.0e-7
    
    #calculate plasma beta
    for N in range(num_strides):
        p_mag = (sfgm_munge[N]['Bbcs'][:,3]*1.0e-9)**2/(2*mu_0)
        beta  = smoms_munge[N]['scalar_p']*1.0e-9/p_mag
        
        smoms_munge[N]['beta'] = beta      

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
        