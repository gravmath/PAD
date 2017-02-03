import datetime           as dt
import matplotlib.cm      as cmap
import matplotlib.dates   as mdates
import matplotlib.patches as patches
import matplotlib.pyplot  as plt
import matplotlib.ticker  as ticker
import numpy              as np
import os
import re
import scipy              as sp
import spacepy.pycdf      as pycdf
import sys

#####################################################################################
def compute_fill_val(array):
    
    if np.max(array) > 0.0:
        fill_val = np.max(temp)
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
def make_FS_summary_plot(obs,curr_debug):
    #unpack the CDF
    times           = np.asarray(curr_debug['Epoch'])
    energies        = np.asarray(curr_debug['%s_des_energy_fast' % obs])[0,:]
    sc_pot          = np.asarray(curr_debug['%s_des_scpot_mean_fast' % obs])
    bx              = np.asarray(curr_debug['%s_des_bentPipeB_X_DSC' % obs])
    by              = np.asarray(curr_debug['%s_des_bentPipeB_Y_DSC' % obs])
    bz              = np.asarray(curr_debug['%s_des_bentPipeB_Z_DSC' % obs])
    bnorm           = np.asarray(curr_debug['%s_des_bentPipeB_Norm' % obs])
    par             = np.asarray(curr_debug['%s_des_energyspectr_par_fast' % obs])
    anti            = np.asarray(curr_debug['%s_des_energyspectr_anti_fast' % obs])
    temp            = np.asarray(curr_debug['%s_des_energyspectr_omni_fast' % obs])
    temp[temp == 0] = 1.0
    omni            = np.log10(temp)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_lowen_fast' % obs])
    temp[temp == 0] = compute_fill_val(temp)
    low             = np.log10(temp)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_miden_fast' % obs])
    temp[temp == 0] = compute_fill_val(temp)
    mid             = np.log10(temp)
    temp            = np.asarray(curr_debug['%s_des_pitchangdist_highen_fast' % obs])
    temp[temp == 0] = compute_fill_val(temp)
    high            = np.log10(temp)
    angles          = np.linspace(0,180,30)
    mean_par  = np.mean(par)
    mean_anti = np.mean(anti)
    #if mean_anti > mean_par:
    #    ratio = anti/par
    #else:
    ratio                  = np.divide(par,anti)
    ratio[np.isnan(ratio)] = 1.0
    ratio[np.isinf(ratio)] = 1.0
    ratio[ratio > 1.4]     = 1.4
    ratio[ratio < 0.6]     = 0.6
    
    FS_time = dt.datetime.strftime(times[0],'%Y-%m-%d_%H%M')
    Rects   = find_sign_intervals(times,bx)
    ###########################################################################
    #create the figure and axes
    fig = plt.figure(figsize=(16,20))
    #fig.autofmt_xdate()
    
    cbar_xloc = 0.925
    ###########################################################################
    #0th pane - spectrogram
    ax0 = fig.add_subplot(611)
    
    #plot the data
    spec_data = ax0.pcolormesh(times,energies,omni.transpose())#,cmap=cmap.bwr)
    
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
    max_exp       = np.ceil(np.max(omni))
    min_exp       = np.floor(np.min(omni))
    log_cbar_span = np.array(range(min_exp,max_exp))
    cb_ax0        = fig.add_axes([cbar_xloc, 0.785714, 0.01, 0.114286])
    spec_cbar     = fig.colorbar(spec_data,ticks=log_cbar_span,cax=cb_ax0)
    cbar_span     = ['%1.1e' % val for val in 10**log_cbar_span]
    cb_ax0.set_yticklabels(cbar_span)
    
    #put on a title
    ax0.set_title('MMS1 on %s' % FS_time)
    
    #add the spacecraft potential
    ax6 = ax0.twinx()
    ax6.plot(times,sc_pot,'k-')
    ax6.set_ylim([0,20])
    #ax6.set_ylabel('SC Potential (V)')
    
    ###########################################################################
    #1st pane - low energy PAD
    ax1 = fig.add_subplot(612)
    
    #plot the data
    PAD_low_ax = ax1.pcolormesh(times,angles,low.transpose(),cmap=cmap.bwr)
    
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
    max_exp       = np.ceil(np.max(low))
    min_exp       = np.floor(np.min(low))
    log_cbar_span = np.array(range(min_exp,max_exp))
    cb_ax1        = fig.add_axes([cbar_xloc,0.648571,0.01,0.114286])
    PAD_low_cbar  = fig.colorbar(PAD_low_ax,ticks=log_cbar_span,cax=cb_ax1)
    cbar_span     = ['%1.1e' % val for val in 10**log_cbar_span]
    cb_ax1.set_yticklabels(cbar_span)

     
    ###########################################################################
    #2nd pane - mid energy PAD
    ax2 = fig.add_subplot(613)
    
    #plot the data
    PAD_mid_ax = ax2.pcolormesh(times,angles,mid.transpose(),cmap=cmap.bwr)
    
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
    max_exp       = np.ceil(np.max(mid))
    min_exp       = np.floor(np.min(mid))
    log_cbar_span = np.array(range(min_exp,max_exp))
    cb_ax2        = fig.add_axes([cbar_xloc,0.511429,0.01,0.114286])
    PAD_mid_cbar  = fig.colorbar(PAD_mid_ax,ticks=log_cbar_span,cax=cb_ax2)
    cbar_span     = ['%1.1e' % val for val in 10**log_cbar_span]
    cb_ax2.set_yticklabels(cbar_span)
    

    ###########################################################################
    #3rd pane - high energy PAD
    ax3 = fig.add_subplot(614)
    
    #plot the data
    PAD_high_ax = ax3.pcolormesh(times,angles,high.transpose(),cmap=cmap.bwr)
    
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
    max_exp       = np.ceil(np.max(high))
    min_exp       = np.floor(np.min(high))
    log_cbar_span = np.array(range(min_exp,max_exp))
    cb_ax3        = fig.add_axes([cbar_xloc,0.374286,0.01,0.114286])
    PAD_high_cbar = fig.colorbar(PAD_high_ax,cax=cb_ax3,ticks=log_cbar_span)
    cbar_span     = ['%1.1e' % val for val in 10**log_cbar_span]
    cb_ax3.set_yticklabels(cbar_span)

    ###########################################################################
    #4th pane - ratio of PADs
    ax4 = fig.add_subplot(615)
    
    #plot the data
    asym_data = ax4.pcolormesh(times,energies,ratio.transpose(),cmap=cmap.bwr)
    #ax4.axhline(50.0,c='k')
    for r in Rects:
        start_x   = r[0]
        start_y   = 50
        width     = r[1]
        thickness = 100
        ax4.add_patch(patches.Rectangle((start_x,start_y),width,thickness,color=r[2],alpha=0.5))
    
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
    max_ratio       = 2#int(np.ceil( np.max(ratio)))
    min_ratio       = 0#int(np.floor(np.min(ratio)))
    cbar_span       = np.asarray(range(min_ratio,max_ratio+1))
    cb_ax4          = fig.add_axes([cbar_xloc,0.237143,0.01,0.114286])
    asym_cbar       = fig.colorbar(asym_data,cax=cb_ax4)#ticks=cbar_span
    #cb_ax4.set_yticklabels(cbar_span)

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
    filename = 'X:/fpishare/Conrad/PAD_plots/FS_summary_%s_%s.png' % (obs,FS_time)
    fig.savefig(filename)
    fig.clf()