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


#set physical constants, including the mass of the species
c    = 299792458.000       #m/s
mu_0 = 4*np.pi*1.0e-7      #N/A^2
ep_0 = 5.52634940620899e7  #eV/(V^2 m)
m_e  = 0.5109989461e6/c**2 #ev/c^2
m_p  = 938.2720813e6/c**2  #ev/c^2
c_e  = 1.6021766208E-19    #electronic charge in Coulombs 

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
        
        #determine the number of NRV values (multiple copies made on 
        #ingestion from the CDF associated with a segment)
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

    if species == 'electrons':
        ms = m_e
    if species == 'ions':
        ms = m_p

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
        speeds      = np.sqrt(2.0*ergs/ms)/1000.0  #division by 1000.0 for km/s
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
def calculate_pitch_angles(sfgm_munge,sdist_munge):
    """Analysis function to calculate the pitch angles of a distribution.

       Arguments:
          sfgm_munge:  a magnetic field munge adapted to the sdist
          sdist_munge: a munge of either electron or ion distribution

       Returns:
           nothing per se - adds to the sdist_munge a field 'pitch_angs' 
           for pitch angles of the distribution
       
       Example use:  
           compute_pitch_angles(fgm_munge_mms1_brst,edist_munge_mms1_brst)
              
       Note:  None
       """    
    #determine the number of strides
    num_sfgm_strides   = len(sfgm_munge)
    num_sdist_strides = len(sdist_munge)
    if num_sfgm_strides != num_sdist_strides:
        print 'Direct comparison of sfgm_munge to sdist_munge not possible.'
        print 'Terminating with extreme prejudice!!!'
        return False
    num_strides = num_sfgm_strides
    
    #adapt the fgm_munge to the sdist_munge
    #sfgm_munge = Munger.adapt_munge_to_munge(fgm_munge,sdist_munge)
    
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
    
    if species == 'electrons':
        ms = m_e
    if species == 'ions':
        ms = m_p
        
    #differential number flux: j_n = 1/2 f v^4/E
    #v^2 = 2 E / m -> j_n = 2(E/m)^2/E = 2E/m^2
    for N in range(num_strides):
        jn = np.zeros(sdist_munge[N]['dist'].shape)
        f  = sdist_munge[N]['dist']
        Es = sdist_munge[N]['ergs']
        for j in range(num_az):
            for k in range(num_polar):
                for m in range(num_ergs):
                    jn[:,j,k,m] = 2.0*f[:,j,k,m]*Es[:,m]/ms**2*(100.0)**4 #(100.0)^4 for cm
                    #print 2.0*Es[:,m]/ms**2
        sdist_munge[N]['jn'] = jn         
        
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
def calculate_plasma_beta(sfgm_munge,smoms_munge):
    """Analysis function to calculate the plasma beta for a particle moments
       munge.

       Arguments:
          sfgm_munge:  a munge for fgm (adapted)
          smoms_munge: a munge of either electron or ion moments

       Returns:
           nothing per se - adds to the smoms_munge a field 'beta' for 
           plasma beta
       
       Example use:  
           calculate_plasma_beta(fgm_munge,emoms_munge)
              
       Note:  None
       """ 
    #determine the number of strides
    num_strides = len(smoms_munge)
    
    #make a compatible fgm_munge
    #sfgm_munge = Munger.adapt_munge_to_munge(fgm_munge,smoms_munge)
    
    #calculate plasma beta
    #note that the factors of 1.0e-9 are needed to get nT -> T and nPa -> Pa
    for N in range(num_strides):
        mag_pres = (sfgm_munge[N]['Bbcs'][:,3]*1.0e-9)**2/(2*mu_0)
        beta     = smoms_munge[N]['scalar_p']*1.0e-9/mag_pres
        
        smoms_munge[N]['beta'] = beta      
        
###############################################################################          
def calculate_plasma_frequency(smoms_munge,species):
    """Analysis function to calculate the plasma frequency for a particle moments
       munge.

       Arguments:
          smoms_munge: a munge of either electron or ion moments
          species:     either 'electrons' or 'ions'

       Returns:
           nothing per se - adds to the smoms_munge a field 'f_ps' for 
           plasma frequency
       
       Example use:  
           calculate_plasma_frequency(emoms_munge,'electrons')
              
       Note:  None
       """ 
       
    #determine the number of strides
    num_strides = len(smoms_munge)
    
    if species == 'electrons':
        ms = m_e
    if species == 'ions':
        ms = m_p
    
    #calculate plasma frequency
    for N in range(num_strides):
        #multiplication of 100.0**3 to get num_den in units of per m^3
        f_ps = np.sqrt( smoms_munge[N]['num_den']*100.0**3/ep_0/ms )/2.0/np.pi
       
        smoms_munge[N]['f_ps'] = f_ps
        
###############################################################################          
def calculate_debye_length(smoms_munge):
    """Analysis function to calculate the Debye length for a particle moments
       munge.

       Arguments:
          smoms_munge: a munge of either electron or ion moments

       Returns:
           nothing per se - adds to the smoms_munge a field 'lambda_D_par' for 
           Debye for T_par, 'lambda_D_perp' for T_perp
       
       Example use:  
           calculate_debye_length(emoms_munge)
              
       Note:  None
       """ 
       
    #determine the number of strides
    num_strides = len(smoms_munge)
    
    #calculate debye length
    for N in range(num_strides):
        #multiplication of 100.0**3 to get num_den in units of per m^3
        lambda_D_par  = np.sqrt( ep_0 * smoms_munge[N]['T_par']  / smoms_munge[N]['num_den']/100.0**3 )
        lambda_D_perp = np.sqrt( ep_0 * smoms_munge[N]['T_perp'] / smoms_munge[N]['num_den']/100.0**3 )
       
        smoms_munge[N]['lambda_D_par']  = lambda_D_par
        smoms_munge[N]['lambda_D_perp'] = lambda_D_perp
        

###############################################################################          
def calculate_cyclotron_frequency(sfgm_munge,smoms_munge,species):
    """Analysis function to calculate the cyclotron frequency for a particle moments
       munge.

       Arguments:
          fgm_munge:   a munge of the magnetic field
          smoms_munge: a munge of either electron or ion moments

       Returns:
           nothing per se - adds to the smoms_munge a field 'f_cs' for 
           cyclotron frequency
       
       Example use:  
           calculate_cyclotron_frequency(fgm_munge,emoms_munge)
              
       Note:  None
       """ 
       
    #determine the number of strides
    num_strides = len(smoms_munge)
    
    if species == 'electrons':
        ms = m_e
    if species == 'ions':
        ms = m_p

    #adapt the fgm_munge to the sdist_munge
    #sfgm_munge = Munger.adapt_munge_to_munge(fgm_munge,smoms_munge)
        
    #calculate plasma frequency
    for N in range(num_strides):
        f_cs = sfgm_munge[N]['Bgse'][:,3]*1e-9/ms/2.0/np.pi
        
        smoms_munge[N]['f_cs'] = f_cs
        
###############################################################################          
def calculate_current(emoms_munge,imoms_munge,species):
    """Analysis function to calculate the current from particle
       measurements

       Arguments:
          emoms_munge: a munge of the electron moments
          imoms_munge: a munge of the ion moments
          species:     string from 'electrons' or 'ions'

       Returns:
           nothing per se - adds the current to the smoms_munge as determined
           by the species string
       
       Example use:  
           calculate_current(emoms_munge,imoms_munge,'electrons')
              
       Note:  None
       """ 
   
    if species == 'electrons':
        #determine the number of strides
        num_strides = len(emoms_munge)
        amoms_munge = Munger.adapt_munge_to_munge(imoms_munge,emoms_munge)
    if species == 'ions':
        #determine the number of strides
        num_strides = len(imoms_munge)
        amoms_munge = Munger.adapt_munge_to_munge(emoms_munge,imoms_munge)

    #calculate current
    for N in range(num_strides):
        J = np.zeros(amoms_munge[N]['bulk_vs'].shape)
        if species == 'electrons':
             J[:,0] = amoms_munge[N]['num_den'][:]*(amoms_munge[N]['bulk_vs'][:,0] - emoms_munge[N]['bulk_vs'][:,0])
             J[:,1] = amoms_munge[N]['num_den'][:]*(amoms_munge[N]['bulk_vs'][:,1] - emoms_munge[N]['bulk_vs'][:,1])
             J[:,2] = amoms_munge[N]['num_den'][:]*(amoms_munge[N]['bulk_vs'][:,2] - emoms_munge[N]['bulk_vs'][:,2])
             emoms_munge[N]['current'] = c_e*J*1e15 #1e15 - to go from cm^-3 and km/s to microamps / m^2
        if species == 'ions':
             J[:,0] = imoms_munge[N]['num_den'][:]*(imoms_munge[N]['bulk_vs'][:,0] - amoms_munge[N]['bulk_vs'][:,0])
             J[:,1] = imoms_munge[N]['num_den'][:]*(imoms_munge[N]['bulk_vs'][:,1] - amoms_munge[N]['bulk_vs'][:,1])
             J[:,2] = imoms_munge[N]['num_den'][:]*(imoms_munge[N]['bulk_vs'][:,2] - amoms_munge[N]['bulk_vs'][:,2])
             imoms_munge[N]['current'] = c_e*J*1e15 #1e15 - to go from cm^-3 and km/s to microamps / m^2
 
###############################################################################  
def calculate_JdotE(adce_munge,smoms_munge):
    """Analysis function to calculate the power delivered (JdotE) 
       from particle measurements

       Arguments:
          adce_munge:  a munge of the electric fields adapted to the 
                       particle times
          smoms_munge: a munge of the species moments

       Returns:
           nothing per se - adds the JdotE to the smoms_munge
       
       Example use:  
           calculate_current(emoms_munge,imoms_munge,'electrons')
              
       Note:  None
       """ 
    #determine number of strides
    num_strides = len(smoms_munge)
    
    for N in range(num_strides):
        JdotE = np.sum(smoms_munge[N]['current']*adce_munge[N]['Egse'],axis=1)
        smoms_munge[N]['JdotE'] = JdotE
             
###############################################################################    
def subtract_internal_photoelectrons(sdist_munge,mode,n_photo,photo_f_file):
    """Analysis function to subtract the instrument internal photoelectrons
       from a electron distribution

       Arguments:
          sdist_munge:   a munge of the electron distributions
          mode:          either 'fast' or 'brst'
          n_photo:       value of n_photo from the moms file
          photo_f_file:  instrument photoelectron model file

       Returns:
           nothing per se - adds a corrected distribution to the sdist_munge
       
       Example use:  
           subtract_internal_photoelectrons(edist_munge,'fast',0.45,photo_f_file)
              
       Note:  None
    
    
    """
    #determine the number of strides
    num_strides = len(sdist_munge)
    
    #import the photoelectron phase space density
    photo_f_cdf = pycdf.CDF(photo_f_file)
    if mode == 'fast':
        photo_mode_str = 'mms_des_bgdist_fast'
        f_photo        = np.asarray(photo_f_cdf[photo_mode_str][:])
    if mode == 'brst':
        photo_mode_str0 = 'mms_des_bgdist_p0_brst'
        photo_mode_str1 = 'mms_des_bgdist_p1_brst'
        f_photo         = {}
        f_photo[0]      = np.asarray(photo_f_cdf[photo_mode_str0][:])
        f_photo[1]      = np.asarray(photo_f_cdf[photo_mode_str1][:])    
      
    for N in range(num_strides):
        #allocate space for the corrected phase-space density data structure
        corrected_dist = np.zeros(sdist_munge[N]['dist'].shape)
    
        #construct startdelphi_count_str
        start_delphi_count_str = 'mms_des_startdelphi_counts_%s' % (mode,)
                
        #loop over all times
        for k in range(len(sdist_munge[N]['epochs'])):
            
            #find start_delphi_index
            startdelphi_index       = sdist_munge[N]['start_dphi'][k]
            
            #construct correction_index    
            correction_index        = int(np.floor(startdelphi_index/16.0))
            
            #subtract off the correction
            if mode == 'fast':
                corrected_dist[k,:,:,:] = sdist_munge[N]['dist'][k,:,:,:] - n_photo*f_photo[correction_index,:,:,:]
            if mode == 'brst':
                parity = sdist_munge[N]['parity'][k]
                corrected_dist[k,:,:,:] = sdist_munge[N]['dist'][k,:,:,:] - n_photo*f_photo[parity][correction_index,:,:,:]
        
        #floor the negative values at zero
        corrected_dist[corrected_dist < 0.0 ] = 0.0
        sdist_munge[N]['cdist'] = corrected_dist                
        
        
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
    theta_rad = np.pi*theta/180.0
    
    #compute spherical cap area based on the formula in Wikipedia
    A = 2*np.pi*(1-np.cos(theta_rad))
    
    return A        
        