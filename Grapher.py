    
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
    filename = filepath+'flux_specturm_survey_'+time_range_str+'_%i.pdf' % time_label
    fig1.savefig(filename,format='pdf',dpi=1200)
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
    filename = filepath+'Smooth_PAD_survey_'+time_range_str+'_%i.pdf' % time_label
    fig1.savefig(filename,format='pdf',dpi=1200,bbox_extra_artists=(lgd,),bbox_inches='tight')    
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
    filename = filepath+'Raw_PAD_survey_'+time_range_str+'_%i.pdf' % time_label
    fig1.savefig(filename,format='pdf',dpi=1200,bbox_extra_artists=(lgd,),bbox_inches='tight')    
    fig1.clf()
    plt.close()


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
    filename = filepath+'Raw_PAD_survey_LM_'+time_range_str+'_%i.pdf' % time_label
    fig1.savefig(filename,format='pdf',dpi=1200,bbox_extra_artists=(lgd,),bbox_inches='tight')
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
    



    