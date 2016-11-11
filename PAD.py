from ipy import *

############################################################################### 
#
# unpack_FS_dist_CDF takes a file handle for a fast survey (FS) distribution 
# file and creates a set of data structures to move the data around.  
# These structures are:
#
# FS_dist:   dictionary with the following keys
#            Dist  - phase space distribution
#            Err   - phase space errors
#            Epoch - name says it all
#            Flag  - goodness flag  
#
# FS_params: dictionary with the following keys
#            Erg   - the fast survey energies
#            Phi   - instrument look direction azimuths
#            Theta - instrument look direction polar angles
#
###############################################################################
def unpack_FS_dist_CDF(cdf_fh,MMS,species,ver,corrections_on,correction_override):
    dist           = cdf_fh['dist']
    dist_str       = '%s_%s_dist_fast'       % (MMS,species)
    dist_err_str   = '%s_%s_disterr_fast'    % (MMS,species)
    energy_str     = '%s_%s_energy_fast'     % (MMS,species)
    error_flag_str = '%s_%s_errorflags_fast' % (MMS,species)
    phi_str        = '%s_%s_phi_fast'        % (MMS,species)
    theta_str      = '%s_%s_theta_fast'      % (MMS,species)

    uncorrected_Dist = np.asarray(dist[dist_str][:,:,:,:])
    
    if corrections_on == 0:
        corrected_Dist   = uncorrected_Dist
    else:
        corrected_Dist   = subtract_internal_photoelectrons(cdf_fh,uncorrected_Dist,MMS,species,correction_override)
    
    FS_dist     = {'Dist'   : corrected_Dist,
                   'Err'    : np.asarray(dist[dist_err_str][:,:,:,:]),
                   'Epoch'  : np.asarray(dist['Epoch'][:]),
                   'Flag'   : np.asarray(dist[error_flag_str][:])}
                   
    if ver == 'ver2':
        Energy = np.asarray(dist[energy_str][:])
    if ver == 'ver3':
        Energy = np.asarray(dist[energy_str][0,:])
                  
    FS_params   = {'Erg'    : Energy,
                   'Phi'    : np.asarray(dist[phi_str][:]),
                   'Theta'  : np.asarray(dist[theta_str][:])}
    
    return FS_dist, FS_params

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
# fetch_AFG_B_field takes the set of CDFs (cdf_fh) and interpolated the AFG
# magnetic field to the FPI points
#
###############################################################################
def fetch_AFG_B_field(cdf_fh,species):
    B       = np.zeros( (len(cdf_fh['dist']['Epoch']),4) )
    B_str   = '%s_afg_srvy_l2pre_dmpa' % species

    counter = 0
    for e in cdf_fh['dist']['Epoch']:
        low_index, high_index = bisect_epochs(e,cdf_fh['AFG']['Epoch'])
        epoch_low    = cdf_fh['AFG']['Epoch'][low_index]
        epoch_high   = cdf_fh['AFG']['Epoch'][high_index]        
        h            = (epoch_high - epoch_low).total_seconds()
        dt           = (epoch_high - e).total_seconds()
        B_low        = cdf_fh['AFG'][B_str][low_index]
        B_high       = cdf_fh['AFG'][B_str][high_index]
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
def fetch_magnetic_field(cdf_fh,MMS,species,source="DEBUG"):
    Bx_str = '%s_%s_bentPipeB_X_DSC'       % (MMS,species)
    By_str = '%s_%s_bentPipeB_Y_DSC'       % (MMS,species)
    Bz_str = '%s_%s_bentPipeB_Z_DSC'       % (MMS,species)
    B_str  = '%s_%s_bentPipeB_Norm'        % (MMS,species)  
    
    #B_field = {'Bx' : np.asarray(cdf_fh[Bx_str][:]),
    #           'By' : np.asarray(cdf_fh[By_str][:]),
    #           'Bz' : np.asarray(cdf_fh[Bz_str][:]),
    #           'B'  : np.asarray(cdf_fh[B_str][:]),}
    if source == "DEBUG":
        Bx     = np.asarray(cdf_fh[Bx_str])
        By     = np.asarray(cdf_fh[By_str])
        Bz     = np.asarray(cdf_fh[Bz_str])
        B      = np.asarray(cdf_fh[B_str])
        B_field = np.vstack((Bx,By,Bz,B)).transpose()
    if source == "AFG":
        B_field = fetch_AFG_B_field(cdf_fh,species)

    return B_field

###############################################################################	
#
# calculate_incoming_particle_directions takes the polar and azimuthal angles
# of the instrument look direction and returns a 2-D array (32,16) of
# the particle velocity directions
#
###############################################################################
def calculate_incoming_particle_directions(parms):
    num_az    = parms['Phi'].shape[0]
    num_polar = parms['Theta'].shape[0]
    deg       = sp.pi/180    

    V = np.zeros((num_az,num_polar,3))
    for i in range(num_az):
        for j in range(num_polar):
            phi   = parms['Phi']  [i]*deg
            theta = parms['Theta'][j]*deg
            V[i,j] = [-sp.sin(theta)*sp.cos(phi),-sp.sin(theta)*sp.sin(phi),-sp.cos(theta)]
            
    return V
	
###############################################################################    
#
# calculate_pitch_angles takes the velocity directions and magnetic field 
# direction at a given time and returns a (32,16) array of pitch angle per 
# pixel
#
###############################################################################
def calculate_pitch_angles(v,bfield,time_label):
    #Construct the B field unit vector at chosen time
    Bx = bfield[time_label,0]
    By = bfield[time_label,1]
    Bz = bfield[time_label,2]
    
    pitch_angles = np.zeros(v.shape[0:2])
    
    for i in range(0,v.shape[0]):
        for j in range(0,v.shape[1]):
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
def compute_counts(FS_dist):
    counts                   = (FS_dist['Dist']/FS_dist['Err'])**2
    counts[np.isnan(counts)] = 0
    counts                   = np.floor(counts)
    return counts

###############################################################################
#
# compute_number_flux the full phase space distribution and the energies and
# returns the number flux
#
###############################################################################
def compute_number_flux(FS_dist,FS_parms):
    m_e    = 9.10938356e-31  #mass of electron in Kg
    eV     = 1.60218e-19     #energy of 1 electron-volt in Joules
    E      = FS_parms['Erg']
    coeff  = 2.0*(E*eV/m_e)**2/E *100.0**4 #100**2 for m -> cm
    jN     = np.zeros(FS_dist['Dist'].shape)
    counter = 0
    for c in coeff:
        jN[:,:,:,counter] = c * FS_dist['Dist'][:,:,:,counter]
        counter           = counter + 1
        
    #jN     = coeff*FS_dist['Dist'] - don't know why this worked in the first place
    return jN   

###############################################################################	
#
# compute_ave_number_flux takes the number flux and gives the azimuthally-
# averaged flux and returns the results
#
###############################################################################
def compute_ave_number_flux(jN,time_label,energy_label):
    ave_jN = np.average(jN[time_label,:,:,energy_label],axis=0)
    
    return ave_jN

###############################################################################    
#
#  load_e_data reads a distribution file and debug file and returns a host of
#  numpy arrays
#
###############################################################################
def load_e_data(cdf_fh,MMS,species,ver,corrections_on,correction_override = 0):
    B                   = cdf_fh['bfield']
    
    edist, parms        = unpack_FS_dist_CDF(cdf_fh,MMS,species,ver,corrections_on,correction_override)
    bfield              = fetch_magnetic_field(B,MMS,species)
    v_dirs              = calculate_incoming_particle_directions(parms)
    counts              = compute_counts(edist)
    jN                  = compute_number_flux(edist,parms)
    
    core_data           = {}
    core_data['edist']  = edist
    core_data['parms']  = parms
    core_data['bfield'] = bfield
    core_data['v_dirs'] = v_dirs
    core_data['counts'] = counts
    core_data['jN']     = jN    
    
    return core_data    
    
###############################################################################    
#
#  subtract_internal_photoelectrons takes a raw phase space distribution
#  and corrects it for internally-generated photoelectons 
#
###############################################################################    
def subtract_internal_photoelectrons(cdf_fh,raw_Dist,MMS,species,correction_override):
    #allocate space for the corrected phase-space density data structure
    corrected_Dist = np.zeros(raw_Dist.shape)
    
    #construct startdelphi_count_str
    start_delphi_count_str = '%s_%s_startdelphi_count_fast' % (MMS,species)
    
    #pull out n_photo - the structure of the following line is understood 
    #by noting that cdf_fh['moms'] gives the moments file and 
    #.attrs['Photoelectron_model_scaling_factor'][0] get the attributes 
    #(attrs) returned as a dictionary from which the key 
    #'Photoelectron_model_scaling_factor' gives the value in some 
    #spacepy/CDF way whose 0th component is a string which is then
    #case to float (whew!!!)
    if correction_override == 0:
        n_photo = float(cdf_fh['moms'].attrs['Photoelectron_model_scaling_factor'][0])
    else:
        n_photo = correction_override
    print cdf_fh.keys(), n_photo
    
    #get the photoelectron model
    f_photo = np.asarray(cdf_fh['photo']['mms_des_bgdist_fast'])
    
    #loop over all times
    for k in range(len(cdf_fh['dist']['Epoch'])):
        #find start_delphi_index
        startdelphi_index       = cdf_fh['dist'][start_delphi_count_str][k]

        #construct correction_index    
        correction_index        = int(np.floor(startdelphi_index/16.0))
        
        #subtract off the correction
        corrected_Dist[k,:,:,:] = raw_Dist[k,:,:,:] - n_photo*f_photo[correction_index,:,:,:]
        
    #return results
    return corrected_Dist
    
    
###############################################################################
#
# compute_limited_PAD returns a truncated PAD limited by energy and pitch 
# angle from the full data
#
###############################################################################
def compute_limited_PAD(time_label,minE,maxE,minPA,maxPA,core_data):
    #find the truncated energy range
    Edata                     = core_data['parms']['Erg'][minE:maxE]
    
    #get the pitch angles and then flatten to a 1-D array
    pitch_angles              = calculate_pitch_angles(core_data['v_dirs'],core_data['bfield'],time_label)
    flat_PA                   = np.ndarray.flatten(pitch_angles)
    
    FAC                       = np.zeros(len(range(minE,maxE)))
    counter                   = 0
    for energy_label in range(minE,maxE):
        local_jN              = np.ndarray.flatten(core_data['jN'][time_label,:,:,energy_label])
        PA_table              = np.array(zip(flat_PA,local_jN))
        PA_range              = np.where( (flat_PA > minPA) & (flat_PA < maxPA))
        FAC[counter]          = np.sum(PA_table[PA_range][:,1])
        counter              += 1        
        
    return Edata, FAC


    
    
###############################################################################    
def create_flux_survey_spectrum(time_label,time_range_str,filepath,core_data):
    minE       = 0
    maxE       = 32
    jN         = core_data['jN']

    omni       = np.zeros(len(range(minE,maxE)))
    counter    = 0
    for energy_label in range(minE,maxE):
        local_jN              = np.ndarray.flatten(jN[time_label,:,:,energy_label])
        omni[counter]         = np.sum(local_jN)/4.0/np.pi
        counter              += 1
    
    fig1,axes  = plt.subplots(nrows=3,ncols=3,figsize=(9,9),sharex=True,sharey=True) 
    plt.subplots_adjust(hspace=0.0, wspace=0.0)
    for i in range(0,9):
        minPA      = i*20.0
        maxPA      = (i+1)*20.0
        sterads    = spherical_cap_area(maxPA) - spherical_cap_area(minPA)
        row        = int(np.floor(i/3))
        col        = i%3
        Edata, FAC = compute_limited_PAD(time_label,minE,maxE,minPA,maxPA,core_data)
        axes[row][col].loglog(Edata,FAC/sterads,'r-',linewidth = 3.0)
        axes[row][col].loglog(Edata,omni,'b-')
        axes[row][col].set_xlim([1e1,1e5])
        axes[row][col].set_ylim([1e2,1e8])
        axes[row][col].xaxis.set_ticks([1e1,1e2,1e3,1e4,1e5])      
        axes[row][col].yaxis.set_ticks([1e2,1e3,1e4,1e5,1e6,1e7,1e8])
        annotation_string = r'$%2.0f\degree \geq \alpha \geq %2.0f\degree$' % (minPA,maxPA)
        axes[row][col].annotate(annotation_string,xy=(2e3,1e7))
    big_ax = fig1.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_color('none')
    big_ax.spines['bottom'].set_color('none')
    big_ax.spines['left'].set_color('none')
    big_ax.spines['right'].set_color('none')
    big_ax.set_ylabel('Differential Number flux (s^-1 cm^-2 eV^-1 sr^-1)',fontsize=14)
    big_ax.set_xlabel('Energy (eV)', fontsize=14)
    title_string = 'Binned PADs for %s\n (red -binned, blue-omnidirectional)' % core_data['edist']['Epoch'][time_label]
    big_ax.set_title(title_string, fontsize=14)
    fig1.tight_layout()    
    filename = filepath+'flux_specturm_survey_'+time_range_str+'_%i.png' % time_label
    fig1.savefig(filename,dpi=300)
    fig1.clf()
    plt.close()

###############################################################################    
def create_smooth_survey_PAD_plot(time_label,time_range_str,filepath,core_data):
    fig1         = plt.figure(figsize=(10,20))
    ax           = fig1.add_subplot(111)
    pitch_angles = calculate_pitch_angles(core_data['v_dirs'],core_data['bfield'],time_label)
    cols         = list(reversed((cm.rainbow(np.linspace(0,1,32)))))
    poly_deg     = 2
    jN_LM        = calculate_flux_LM(time_label,core_data)
    ELow         = 0
    EHigh        = 32
    for e in range(ELow,EHigh):
        A        = jN_LM[:,:,e]
        my_poly  = np.polyfit(A[np.argsort(A[:,0])][:,0],A[np.argsort(A[:,0])][:,2],poly_deg)
        Q        = np.polyval(my_poly,A[np.argsort(A[:,0])][:,0])
        plt.semilogy(A[np.argsort(A[:,0])][:,0],Q,color=cols[e],marker='None',linewidth=1.0,label = "%2.2f eV" % (core_data['parms']['Erg'][e]) )    
    plt.xlabel('Pitch angle (deg)',fontsize=14)
    plt.ylabel('Fitted Differential Number Flux (cm^-2 s^-1 eV^-1)',fontsize=14)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='lower left',bbox_to_anchor=(1, 0.5))
    plt.title(core_data['edist']['Epoch'][time_label])
    plt.grid(b=True, which='major', color='gray', linestyle='-')
    filename = filepath+'Smooth_PAD_survey_'+time_range_str+'_%i.png' % time_label
    fig1.savefig(filename,dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')    
    fig1.clf()
    plt.close()

###############################################################################
def create_raw_survey_PAD_plot(time_label,Elow,Ehigh,time_range_str,filepath,core_data):
    fig1                   = plt.figure(figsize=(10,20))
    ax                     = fig1.add_subplot(111)
    pitch_angles           = calculate_pitch_angles(core_data['v_dirs'],core_data['bfield'],time_label)
    flat_pitch_angles      = np.ndarray.flatten(pitch_angles)
    n_energies             = Ehigh - Elow + 1
    cols                   = list(reversed((cm.rainbow(np.linspace(0,1,n_energies)))))
    for energy_label in range(Elow,Ehigh):
        flat_particle_flux = np.ndarray.flatten(core_data['jN'][time_label,:,:,energy_label])  
        plt.semilogy(flat_pitch_angles,flat_particle_flux,color=cols[energy_label],marker='.',linewidth=0.0,label = "%2.2f eV" % (core_data['parms']['Erg'][energy_label]) )
    plt.xlabel('Pitch angle (deg)',fontsize=14)
    plt.ylabel('PSD s^3/cm^6',fontsize=14)
    plt.ylabel('Differential Number Flux (cm^-2 s^-1 eV^-1)',fontsize=14)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='lower left',bbox_to_anchor=(1, 0.5))
    plt.title(core_data['edist']['Epoch'][time_label])
    plt.grid(b=True, which='major', color='gray', linestyle='-')
    plt.grid(b=True, which='minor', color='gray', linestyle='--')
    filename = filepath+'Raw_PAD_survey_'+time_range_str+'_%i.png' % time_label
    fig1.savefig(filename,dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')    
    fig1.clf()
    plt.close()

###############################################################################
def convert_to_LM(B):
    B_norm = np.sqrt(B.dot(B))
    z_LM   =  B/B_norm
    z_sc   = np.array([0,0,1])
    x      = np.cross(z_sc,B)
    x_norm = np.sqrt(x.dot(x))
    x_LM   = x/x_norm
    y_LM   = np.cross(z_LM,x_LM)
    
    return np.vstack((x_LM,y_LM,z_LM))   

###############################################################################
def calculate_flux_LM(time_label,core_data):
    B          = core_data['bfield'][time_label,0:3]
    T          = convert_to_LM(B)
    num_pixels = len(core_data['parms']['Phi'])*len(core_data['parms']['Theta'])
    num_ergs   = len(core_data['parms']['Erg'])
    
    jN_LM      = np.zeros((num_pixels,3,num_ergs))    
    
    counter = 0
    for i in range(32):
        for j in range(16):
            v_dir_LM           = T.dot(core_data['v_dirs'][i][j])
            jN_LM[counter,0,:] = np.arccos(v_dir_LM[2])*180/np.pi
            jN_LM[counter,1,:] = np.arctan2(v_dir_LM[1],v_dir_LM[0])*180/np.pi
            jN_LM[counter,2,:] = core_data['jN'][time_label,i,j,:]
            counter             += 1
            
    return jN_LM    

###############################################################################    
def create_raw_survey_PAD_plot_LM(time_label,Elow,Ehigh,time_range_str,filepath,core_data):
    jN_LM = calculate_flux_LM(time_label,core_data)
    fig1  = plt.figure(figsize=(10,20))
    ax    = fig1.add_subplot(111) 
    cols                   = list(reversed((cm.rainbow(np.linspace(0,1,32)))))
    for e in range(Elow,Ehigh):
        plt.semilogy(jN_LM[:,0,e],jN_LM[:,2,e],color=cols[e],marker='.',linewidth=0.0,label = "%2.2f eV" % (core_data['parms']['Erg'][e]) )    
    plt.xlabel('Pitch angle (deg)',fontsize=14)
    plt.ylabel('Differential Number Flux (cm^-2 s^-1 eV^-1)',fontsize=14)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='lower left',bbox_to_anchor=(1, 0.5))
    plt.title(core_data['edist']['Epoch'][time_label])
    plt.grid(b=True, which='major', color='gray', linestyle='-')
    plt.grid(b=True, which='minor', color='gray', linestyle='--')   
    filename = filepath+'Raw_PAD_survey_LM_'+time_range_str+'_%i.png' % time_label
    fig1.savefig(filename,dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')
    fig1.clf()
    plt.close()

###############################################################################    
def visualize_FPI_pixels_in_LM(time_label,time_range_str,filepath,core_data):
    jN_LM = calculate_flux_LM(time_label,core_data)

    fig1 = plt.figure(figsize=(10,10))
    plt.plot(jN_LM[:,1,0],jN_LM[:,0,0],'b.')
    plt.grid('on')
    plt.xlabel('LM Azimuth (deg)',fontsize = 14)
    plt.ylabel('LM Polar or Pitch Angle (deg)', fontsize = 14)
    bin20  = np.where(   (0.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<= 20.0) )[0]
    bin40  = np.where(  (20.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<= 40.0) )[0]
    bin60  = np.where(  (40.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<= 60.0) )[0]
    bin80  = np.where(  (60.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<= 80.0) )[0]
    bin100 = np.where(  (80.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<=100.0) )[0]
    bin120 = np.where( (100.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<=120.0) )[0]
    bin140 = np.where( (120.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<=140.0) )[0]
    bin160 = np.where( (140.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<=160.0) )[0]
    bin180 = np.where( (160.0<=jN_LM[:,0,0]) & (jN_LM[:,0,0]<=180.0) )[0]
    plt.annotate('%d pts' % len(bin20), xy=(-190,170))
    plt.annotate('%d pts' % len(bin40), xy=(-190,150))
    plt.annotate('%d pts' % len(bin60), xy=(-190,130))
    plt.annotate('%d pts' % len(bin80), xy=(-190,110))
    plt.annotate('%d pts' % len(bin100),xy=(-190,90))
    plt.annotate('%d pts' % len(bin120),xy=(-190,70))
    plt.annotate('%d pts' % len(bin140),xy=(-190,50))
    plt.annotate('%d pts' % len(bin160),xy=(-190,30))
    plt.annotate('%d pts' % len(bin180),xy=(-190,10))

    plt.annotate('0.38 sr',xy=(170,170))
    plt.annotate('1.09 sr',xy=(170,150))
    plt.annotate('1.67 sr',xy=(170,130))
    plt.annotate('2.05 sr',xy=(170,110))
    plt.annotate('2.18 sr',xy=(170,90))
    plt.annotate('2.05 sr',xy=(170,70))
    plt.annotate('1.67 sr',xy=(170,50))
    plt.annotate('1.09 sr',xy=(170,30))
    plt.annotate('0.38 sr',xy=(170,10))
    plt.title('FPI look directions in LM Coordinates at %s' % core_data['edist']['Epoch'][time_label])
    filename = filepath+'FPI_look_angles_LM_'+time_range_str+'_%i.png' % time_label
    fig1.savefig(filename,dpi=300)
    fig1.clf()
    plt.close()
    

############################################################################### 
def spherical_cap_area(theta):
    #put theta into radians
    theta_loc = np.pi*theta/180.0
    
    #calculate spherical cap area based on the formula in Wikipedia
    A = 2*np.pi*(1-np.cos(theta_loc))
    
    return A    
    
###############################################################################	
def unpack_FS_moms_CDF(cdf_fh,MMS,species):
    vx_str  = '%s_%s_bulkx_dbcs_fast'     % (MMS,species)
    vy_str  = '%s_%s_bulky_dbcs_fast'     % (MMS,species)
    vz_str  = '%s_%s_bulkz_dbcs_fast'     % (MMS,species)    
    v_str   = '%s_%s_bulkspeed_dbcs_fast' % (MMS,species)    
    
    FS_moms = {'Vx' : np.asarray(cdf_fh[vx_str]),
               'Vy' : np.asarray(cdf_fh[vy_str]),
               'Vz' : np.asarray(cdf_fh[vz_str]),
               'V'  : np.asarray(cdf_fh[v_str])}
    
    return FS_moms

###############################################################################	
def calculate_instrument_average_pitch_angles(pitch_angles):
    num_az           = pitch_angles.shape[0]
    ave_pitch_angles = np.average(pitch_angles,axis=0) #average over azimuth
    std_pitch_angles = np.std    (pitch_angles,axis=0) #std over azimuth
    
    return ave_pitch_angles, std_pitch_angles 


###############################################################################	
def fetch_PSD(FS_dist,time_label,energy_label):
    
    return FS_dist['Dist'][time_label,:,:,energy_label]

###############################################################################	
def fetch_ave_PSD(FS_dist,time_label,energy_label):
    f     = fetch_PSD(FS_dist,time_label,energy_label)
    ave_f = np.average(f,axis=0)
    
    return ave_f
    