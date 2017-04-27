import datetime           as dt
import matplotlib.cm      as cmap
import matplotlib.dates   as mdates
import matplotlib.patches as patches
import matplotlib.pyplot  as plt
import matplotlib.ticker  as ticker
import numpy              as np

    
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
    

#####################################################################################
def compute_fill_val(array):
    
    if np.max(array) > 0.0:
        fill_val = np.max(array)
    else:
        fill_val = 1.0
        
    return fill_val

#####################################################################################
def find_sign_intervals(times,signal):
    curr_index = 0
    curr_sign  = np.sign(signal[0])
    counter    = 0
    intervals  = [[curr_index,curr_sign]]
    for s in signal:
        s_index = counter 
        s_sign  = np.sign(s)
        counter = counter + 1
        if s_sign != curr_sign:
            intervals.append([s_index,s_sign])
            curr_index = s_index
            curr_sign  = s_sign
    #finally clean up the end
    intervals.append([s_index,s_sign])
    
    Rects = []
    for i in range(len(intervals)-1):
        start_index = intervals[i][0]
        stop_index  = intervals[i+1][0]
        start_time  = mdates.date2num(times[start_index])
        stop_time   = mdates.date2num(times[stop_index])
        time_delta  = stop_time - start_time
        if intervals[i][1] == -1.0:
            color = 'k'
        elif intervals[i][1] == +1.0:
            color = 'g'
        else:
            color = 'y'
        Rects.append([start_time,time_delta,color])
    
    return Rects

#####################################################################################
#  Summary Plots
#####################################################################################
def make_brst_summary_plot(obs,curr_debug):
    #unpack the CDF
    times           = np.asarray(curr_debug['Epoch'])
    energies        = np.asarray(curr_debug['%s_des_energy_brst' % obs])[0,:]
    sc_pot          = np.asarray(curr_debug['%s_des_scpot_mean_brst' % obs])
    bx              = np.asarray(curr_debug['%s_des_bentPipeB_X_DSC' % obs])
    by              = np.asarray(curr_debug['%s_des_bentPipeB_Y_DSC' % obs])
    bz              = np.asarray(curr_debug['%s_des_bentPipeB_Z_DSC' % obs])
    bnorm           = np.asarray(curr_debug['%s_des_bentPipeB_Norm' % obs])
    par             = np.asarray(curr_debug['%s_des_energyspectr_par_brst' % obs]).T
    anti            = np.asarray(curr_debug['%s_des_energyspectr_anti_brst' % obs]).T
    temp            = np.asarray(curr_debug['%s_des_energyspectr_omni_brst' % obs])
    omni            = np.ma.masked_invalid(np.log10(temp).T)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_lowen_brst' % obs])
    low             = np.ma.masked_invalid(np.log10(temp).T)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_miden_brst' % obs])
    mid             = np.ma.masked_invalid(np.log10(temp).T)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_highen_brst' % obs])
    high            = np.ma.masked_invalid(np.log10(temp).T)
    angles          = np.linspace(0,180,30)
    mean_par        = np.mean(par)
    mean_anti       = np.mean(anti)
    ratio           = np.ma.masked_invalid(np.divide(par,anti))
    
    BRST_time = dt.datetime.strftime(times[0],'%Y-%m-%d_%H%M')
    
    ###########################################################################
    #create the figure and axes
    fig = plt.figure(figsize=(16,20))
    #fig.autofmt_xdate()
    
    cbar_off = 0.01
    cbar_wth = 0.01
    cmap.jet.set_bad('k',alpha=1.0)    
    cmap.bwr.set_bad('k',alpha=1.0)    
    
    ###########################################################################
    #0th pane - spectrogram
    min_exp = 4
    max_exp = 8
    ax0     = fig.add_subplot(611)
    
    #plot the data
    spec_data = ax0.pcolormesh(times,energies,omni,cmap=cmap.jet,vmin=max_exp,vmax=min_exp)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax0.xaxis.set_major_locator(minutes)
    ax0.xaxis.set_major_formatter(time_format)
    ax0.set_xlabel('Time')
    
    #deal with the y-axis
    ax0.set_yscale('log')
    ax0.set_ylim([energies[0],energies[31]])
    ax0.set_ylabel('Energy (eV)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax0        = fig.add_axes(cbar_position(ax0,cbar_off,cbar_wth))
    spec_cbar     = fig.colorbar(spec_data,cax=cb_ax0,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10.0**log_cbar_span]
    #cb_ax0.set_yticklabels(cbar_span)
    
    #put on a title
    ax0.set_title('%s on %s' % (obs.upper(),BRST_time))
    
    ax0.plot(times,sc_pot,'k-')
    
    ###########################################################################
    #1st pane - low energy PAD
    min_exp = int(np.floor(np.min(low)))
    max_exp = int(np.ceil(np.max(low)))
    ax1 = fig.add_subplot(612)
    
    #plot the data
    PAD_low_dat = ax1.pcolormesh(times,angles,low,cmap=cmap.bwr,vmin=min_exp,vmax=max_exp)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax1.xaxis.set_major_locator(minutes)
    ax1.xaxis.set_major_formatter(time_format)
    ax1.set_xlabel('Time')
    
    #deal with the y-axis
    ax1.set_yscale('linear')
    ax1.set_ylim([0,180])
    ax1.set_ylabel('Low Pitch Ang (deg)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax1        = fig.add_axes(cbar_position(ax1,cbar_off,cbar_wth))
    PAD_low_cbar  = fig.colorbar(PAD_low_dat,cax=cb_ax1,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10**log_cbar_span]
    #cb_ax1.set_yticklabels(cbar_span)

     
    ###########################################################################
    #2nd pane - mid energy PAD
    min_exp = int(np.floor(np.min(mid)))
    max_exp = int(np.ceil(np.max(mid))) 
    ax2 = fig.add_subplot(613)
    
    #plot the data
    PAD_mid_dat = ax2.pcolormesh(times,angles,mid,cmap=cmap.bwr)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax2.xaxis.set_major_locator(minutes)
    ax2.xaxis.set_major_formatter(time_format)
    ax2.set_xlabel('Time')
    
    #deal with the y-axis
    ax2.set_yscale('linear')
    ax2.set_ylim([0,180])
    ax2.set_ylabel('Mid Pitch Ang (deg)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax2        = fig.add_axes(cbar_position(ax2,cbar_off,cbar_wth))
    PAD_mid_cbar  = fig.colorbar(PAD_mid_dat,cax=cb_ax2,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10**log_cbar_span]
    #cb_ax2.set_yticklabels(cbar_span)
    

    ###########################################################################
    #3rd pane - high energy PAD
    min_exp = int(np.floor(np.min(high)))
    max_exp = int(np.ceil(np.max(high))) 
    ax3 = fig.add_subplot(614)
    
    #plot the data
    PAD_high_dat = ax3.pcolormesh(times,angles,high,cmap=cmap.bwr)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax3.xaxis.set_major_locator(minutes)
    ax3.xaxis.set_major_formatter(time_format)
    ax3.set_xlabel('Time')
    
    #deal with the y-axis
    ax3.set_yscale('linear')
    ax3.set_ylim([0,180])
    ax3.set_ylabel('Hi Pitch Ang (deg)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax3        = fig.add_axes(cbar_position(ax3,cbar_off,cbar_wth))
    PAD_high_cbar = fig.colorbar(PAD_high_dat,cax=cb_ax3,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10**log_cbar_span]
    #cb_ax3.set_yticklabels(cbar_span)

    ###########################################################################
    #4th pane - ratio of PADs
    ax4 = fig.add_subplot(615)
    
    #plot the data
    asym_data = ax4.pcolormesh(times,energies,ratio,cmap=cmap.bwr,vmin=0.6,vmax=1.4)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax4.xaxis.set_major_locator(minutes)
    ax4.xaxis.set_major_formatter(time_format)
    ax4.set_xlabel('Time')
    
    #deal with the y-axis
    ax4.set_yscale('log')
    ax4.set_ylim([energies[0],energies[31]])
    ax4.set_ylabel('Energy (eV)')
    
    #create the colorbar
    cb_ax4          = fig.add_axes(cbar_position(ax4,cbar_off,cbar_wth))
    asym_cbar       = fig.colorbar(asym_data,cax=cb_ax4)

    ###########################################################################
    #5th pane - mag field
    ax5 = fig.add_subplot(616)
    
    #plot the data
    ax5.plot(times,bx*bnorm,label='Bx')
    ax5.plot(times,by*bnorm,label='By')
    ax5.plot(times,bz*bnorm,label='Bz')
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax5.xaxis.set_major_locator(minutes)
    ax5.xaxis.set_major_formatter(time_format)
    ax5.set_xlabel('Time')
    
    #deal with the y-axis
    ax5.set_yscale('linear')
    ax5.set_ylabel('Magnetic Field (nT)')
    
    ax5.legend(bbox_to_anchor=[1.12,1.0])
    ###########################################################################
    #save the file
    
    #fig.tight_layout()
    filename = 'x:/fpishare/Conrad/PAD_plots/BRST_summary_%s_%s.png' % (obs,BRST_time)
    fig.savefig(filename)
    fig.clf()

#####################################################################################
def make_fast_summary_plot(obs,curr_debug):
    #unpack the CDF
    times           = np.asarray(curr_debug['Epoch'])
    energies        = np.asarray(curr_debug['%s_des_energy_fast' % obs])[0,:]
    sc_pot          = np.asarray(curr_debug['%s_des_scpot_mean_fast' % obs])
    bx              = np.asarray(curr_debug['%s_des_bentPipeB_X_DSC' % obs])
    by              = np.asarray(curr_debug['%s_des_bentPipeB_Y_DSC' % obs])
    bz              = np.asarray(curr_debug['%s_des_bentPipeB_Z_DSC' % obs])
    bnorm           = np.asarray(curr_debug['%s_des_bentPipeB_Norm' % obs])
    par             = np.asarray(curr_debug['%s_des_energyspectr_par_fast' % obs]).T
    anti            = np.asarray(curr_debug['%s_des_energyspectr_anti_fast' % obs]).T
    temp            = np.asarray(curr_debug['%s_des_energyspectr_omni_fast' % obs])
    omni            = np.ma.masked_invalid(np.log10(temp).T)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_lowen_fast' % obs])
    low             = np.ma.masked_invalid(np.log10(temp).T)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_miden_fast' % obs])
    mid             = np.ma.masked_invalid(np.log10(temp).T)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_highen_fast' % obs])
    high            = np.ma.masked_invalid(np.log10(temp).T)
    angles          = np.linspace(0,180,30)
    mean_par        = np.mean(par)
    mean_anti       = np.mean(anti)
    ratio           = np.ma.masked_invalid(np.divide(par,anti))
    
    FAST_time = dt.datetime.strftime(times[0],'%Y-%m-%d_%H%M')
    
    ###########################################################################
    #create the figure and axes
    fig = plt.figure(figsize=(16,20))
    #fig.autofmt_xdate()
    
    cbar_off = 0.01
    cbar_wth = 0.01
    cmap.jet.set_bad('k',alpha=1.0)    
    cmap.bwr.set_bad('k',alpha=1.0)    
    
    ###########################################################################
    #0th pane - spectrogram
    min_exp = 4
    max_exp = 8
    ax0     = fig.add_subplot(611)
    
    #plot the data
    spec_data = ax0.pcolormesh(times,energies,omni,cmap=cmap.jet,vmin=max_exp,vmax=min_exp)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax0.xaxis.set_major_locator(minutes)
    ax0.xaxis.set_major_formatter(time_format)
    ax0.set_xlabel('Time')
    
    #deal with the y-axis
    ax0.set_yscale('log')
    ax0.set_ylim([energies[0],energies[31]])
    ax0.set_ylabel('Energy (eV)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax0        = fig.add_axes(cbar_position(ax0,cbar_off,cbar_wth))
    spec_cbar     = fig.colorbar(spec_data,cax=cb_ax0,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10.0**log_cbar_span]
    #cb_ax0.set_yticklabels(cbar_span)
    
    #put on a title
    ax0.set_title('%s on %s' % (obs.upper(),FAST_time))
    
    ax0.plot(times,sc_pot,'k-')
    
    ###########################################################################
    #1st pane - low energy PAD
    min_exp = int(np.floor(np.min(low)))
    max_exp = int(np.ceil(np.max(low)))
    ax1 = fig.add_subplot(612)
    
    #plot the data
    PAD_low_dat = ax1.pcolormesh(times,angles,low,cmap=cmap.bwr,vmin=min_exp,vmax=max_exp)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax1.xaxis.set_major_locator(minutes)
    ax1.xaxis.set_major_formatter(time_format)
    ax1.set_xlabel('Time')
    
    #deal with the y-axis
    ax1.set_yscale('linear')
    ax1.set_ylim([0,180])
    ax1.set_ylabel('Low Pitch Ang (deg)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax1        = fig.add_axes(cbar_position(ax1,cbar_off,cbar_wth))
    PAD_low_cbar  = fig.colorbar(PAD_low_dat,cax=cb_ax1,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10**log_cbar_span]
    #cb_ax1.set_yticklabels(cbar_span)

     
    ###########################################################################
    #2nd pane - mid energy PAD
    min_exp = int(np.floor(np.min(mid)))
    max_exp = int(np.ceil(np.max(mid))) 
    ax2 = fig.add_subplot(613)
    
    #plot the data
    PAD_mid_dat = ax2.pcolormesh(times,angles,mid,cmap=cmap.bwr)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax2.xaxis.set_major_locator(minutes)
    ax2.xaxis.set_major_formatter(time_format)
    ax2.set_xlabel('Time')
    
    #deal with the y-axis
    ax2.set_yscale('linear')
    ax2.set_ylim([0,180])
    ax2.set_ylabel('Mid Pitch Ang (deg)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax2        = fig.add_axes(cbar_position(ax2,cbar_off,cbar_wth))
    PAD_mid_cbar  = fig.colorbar(PAD_mid_dat,cax=cb_ax2,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10**log_cbar_span]
    #cb_ax2.set_yticklabels(cbar_span)
    

    ###########################################################################
    #3rd pane - high energy PAD
    min_exp = int(np.floor(np.min(high)))
    max_exp = int(np.ceil(np.max(high))) 
    ax3 = fig.add_subplot(614)
    
    #plot the data
    PAD_high_dat = ax3.pcolormesh(times,angles,high,cmap=cmap.bwr)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax3.xaxis.set_major_locator(minutes)
    ax3.xaxis.set_major_formatter(time_format)
    ax3.set_xlabel('Time')
    
    #deal with the y-axis
    ax3.set_yscale('linear')
    ax3.set_ylim([0,180])
    ax3.set_ylabel('Hi Pitch Ang (deg)')
    
    #create the colorbar
    log_cbar_span = np.array(range(min_exp,max_exp+1))
    cb_ax3        = fig.add_axes(cbar_position(ax3,cbar_off,cbar_wth))
    PAD_high_cbar = fig.colorbar(PAD_high_dat,cax=cb_ax3,ticks=log_cbar_span,format=ticker.FormatStrFormatter('$10^{%d}$'))
    #cbar_span     = ['%1.0e' % val for val in 10**log_cbar_span]
    #cb_ax3.set_yticklabels(cbar_span)

    ###########################################################################
    #4th pane - ratio of PADs
    ax4 = fig.add_subplot(615)
    
    #plot the data
    asym_data = ax4.pcolormesh(times,energies,ratio,cmap=cmap.bwr,vmin=0.6,vmax=1.4)
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax4.xaxis.set_major_locator(minutes)
    ax4.xaxis.set_major_formatter(time_format)
    ax4.set_xlabel('Time')
    
    #deal with the y-axis
    ax4.set_yscale('log')
    ax4.set_ylim([energies[0],energies[31]])
    ax4.set_ylabel('Energy (eV)')
    
    #create the colorbar
    cb_ax4          = fig.add_axes(cbar_position(ax4,cbar_off,cbar_wth))
    asym_cbar       = fig.colorbar(asym_data,cax=cb_ax4)

    ###########################################################################
    #5th pane - mag field
    ax5 = fig.add_subplot(616)
    
    #plot the data
    ax5.plot(times,bx*bnorm,label='Bx')
    ax5.plot(times,by*bnorm,label='By')
    ax5.plot(times,bz*bnorm,label='Bz')
    
    #deal with epochs on the x-axis
    time_format = mdates.DateFormatter('%H:%M')
    minutes     = mdates.MinuteLocator(range(0,59),interval = 10,tz=None)
    ax5.xaxis.set_major_locator(minutes)
    ax5.xaxis.set_major_formatter(time_format)
    ax5.set_xlabel('Time')
    
    #deal with the y-axis
    ax5.set_yscale('linear')
    ax5.set_ylabel('Magnetic Field (nT)')
    
    ax5.legend(bbox_to_anchor=[1.12,1.0])
    ###########################################################################
    #save the file
    
    #fig.tight_layout()
    filename = 'x:/fpishare/Conrad/PAD_plots/FAST_summary_%s_%s.png' % (obs,FAST_time)
    fig.savefig(filename)
    fig.clf()
    
#####################################################################################   
def cbar_position(current_ax,offset,cbar_width):
    #get the tuple describing the lower x,y corner position
    #and the width and height
    x_pos, y_pos, width, height = current_ax.get_position().bounds
    return [x_pos + width + offset, y_pos, cbar_width, height]    