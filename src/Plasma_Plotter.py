# -*- coding: utf-8 -*-
import datetime           as dt
import matplotlib.cm      as cmap
import matplotlib.dates   as mdates
import matplotlib.pyplot  as plt
import matplotlib.ticker  as ticker
import numpy              as np
import scipy.interpolate  as interp

import Convert
import Grapher
import PAD

###############################################################################
def add_info_box(fig,date,geometry):
    dax = fig.add_axes(geometry)
    Grapher.quiet_axis(dax)
    string = 'hh:mm\nX-GSM (Re)\nY-GSM (Re)\nZ-GSM (Re)\n%s' % date
    dax.annotate(string,xy=(0.5,0.5))    

###############################################################################   
def pos_on_time_axis(fig,ax,cursor,fpi_prd1,obs,mode,descriptor,year,month,day):
    x_interp, y_interp, z_interp = Grapher.Convert.construct_interpolants\
                                   (cursor,fpi_prd1,obs,mode,\
                                    descriptor,year,month,day)
    fig.canvas.draw_idle()
    xticks_labs = ax.get_xticklabels()
    for xt in xticks_labs:
        xt_text              = xt.get_text()
        hour, minute, second = [int(t) for t in xt_text.split(':')]
        curr_epoch           = dt.datetime(year,month,day,hour,minute,second)
        curr_time            = mdates.date2num(curr_epoch)
        if x_interp != False:
            x                = interp.splev(curr_time,x_interp)
        else:
            x                = 0
        if y_interp != False:
            y                = interp.splev(curr_time,y_interp)
        else:
            y                = 0
        if z_interp != False:
            z                = interp.splev(curr_time,z_interp)        
        else:
            z                = 0
        xt.set_text(xt_text+'\n%2.1f\n%2.1f\n%2.1f' % (x,y,z))  
    ax.set_xticklabels(xticks_labs)
           
###############################################################################
def make_density_panel(ax,obs,emoms_munge,imoms_munge):
    #determine the number of segments 
    num_estrides = len(emoms_munge)
    num_istrides = len(imoms_munge)
    if num_estrides != num_istrides:
        print "Burst segments for electrons and ions don't match!"
        print "Terminating with extreme prejudice!!!"
        return False
    else:
        num_strides = num_estrides
       
    #graph and label the first segment
    te1 = emoms_munge[0]['epochs']
    ne1 = emoms_munge[0]['num_den']
    ti1 = imoms_munge[0]['epochs']
    ni1 = imoms_munge[0]['num_den']
    n_trace = Grapher.curves(ax,te1,ne1)
    n_trace.add_line(ti1,ni1)
    n_trace.customize_li(0,{'color':'black','label':'Ne'})
    n_trace.customize_li(1,{'color':'green','label':'Ni'})
    
    for j in range(1,num_strides):
        tej = emoms_munge[j]['epochs']
        nej = emoms_munge[j]['num_den']
        tij = imoms_munge[j]['epochs']
        nij = imoms_munge[j]['num_den']
        n_trace.add_line(tej,nej)
        n_trace.add_line(tij,nij)
        n_trace.customize_li(2*j,{'color':'black'})
        n_trace.customize_li(2*j+1,{'color':'green'})
        
    n_trace.customize_ax({'ylabel':'%s\\ndensity\\n[$cm^{-3}$]'%obs})
    n_trace.show_legend()
    return n_trace
    
###############################################################################
def make_beta_panel(ax,obs,emoms_munge,imoms_munge):
    #determine the number of segments 
    num_estrides = len(emoms_munge)
    num_istrides = len(imoms_munge)
    if num_estrides != num_istrides:
        print "Burst segments for electrons and ions don't match!"
        print "Terminating with extreme prejudice!!!"
        return False
    else:
        num_strides = num_estrides
       
    #graph and label the first segment
    te1 = emoms_munge[0]['epochs']
    be1 = emoms_munge[0]['beta']
    ti1 = imoms_munge[0]['epochs']
    bi1 = imoms_munge[0]['beta']
    b_trace = Grapher.curves(ax,te1,be1)
    b_trace.add_line(ti1,bi1)
    b_trace.customize_li(0,{'color':'black','label':'Beta_e'})
    b_trace.customize_li(1,{'color':'green','label':'Beta_i'})
    
    for j in range(1,num_strides):
        tej = emoms_munge[j]['epochs']
        bej = emoms_munge[j]['beta']
        tij = imoms_munge[j]['epochs']
        bij = imoms_munge[j]['beta']
        b_trace.add_line(tej,bej)
        b_trace.add_line(tij,bij)
        b_trace.customize_li(2*j,{'color':'black'})
        b_trace.customize_li(2*j+1,{'color':'green'})
        
    b_trace.customize_ax({'ylabel':'%s\\nBeta'%obs,'yscale':'log'})
    b_trace.show_legend()
    return b_trace    

###############################################################################
def make_sdensity_panel(ax,obs,smoms_munge,species):
    #determine the number of strides 
    num_strides = len(smoms_munge)
       
    #graph and label the first segment
    ts1 = smoms_munge[0]['epochs']
    ns1 = smoms_munge[0]['num_den']
    n_trace = Grapher.curves(ax,ts1,ns1)
    if species == 'emoms':
        n_trace.customize_li(0,{'color':'black','label':'Ne'})
    if species == 'imoms':
        n_trace.customize_li(1,{'color':'green','label':'Ni'})
    
    for j in range(1,num_strides):
        tsj = emoms_munge[j]['epochs']
        nsj = emoms_munge[j]['num_den']
        n_trace.add_line(tsj,nsj)
        if species == 'emoms':
            n_trace.customize_li(2*j,{'color':'black'})
        if species == 'imoms':
            n_trace.customize_li(2*j+1,{'color':'green'})
    
    if species == 'emoms':
        n_trace.customize_ax({'ylabel':'%s\\n$N_e$\\n[$cm^{-3}$]'%obs})
    if species == 'imoms':
        n_trace.customize_ax({'ylabel':'%s\\n$N_i$\\n[$cm^{-3}$]'%obs})
    #n_trace.show_legend()
    return n_trace    
    
###############################################################################
def make_temperature_panel(ax,obs,emoms_munge,imoms_munge):
    #determine the number of segments
    num_estrides = len(emoms_munge)
    num_istrides = len(imoms_munge)
    if num_estrides != num_istrides:
        print "Burst segments for electrons and ions don't match!"
        print "Terminating with extreme prejudice!!!"
        return False
    else:
        num_strides = num_estrides
       
    #graph and label the first segment
    te1      = emoms_munge[0]['epochs']
    Te_perp1 = emoms_munge[0]['T_perp']
    Te_par1  = emoms_munge[0]['T_par']
    ti1      = imoms_munge[0]['epochs']
    Ti_perp1 = imoms_munge[0]['T_perp']
    Ti_par1  = imoms_munge[0]['T_par']
    T_trace = Grapher.curves(ax,te1,Te_perp1)
    T_trace.add_line(te1,Te_par1)
    T_trace.add_line(ti1,Ti_perp1)
    T_trace.add_line(ti1,Ti_par1)
    T_trace.customize_li(0,{'color':'red','label':'$T_{e\perp}$'})
    T_trace.customize_li(1,{'color':'blue','label':'$T_{e\parallel}$'})
    T_trace.customize_li(2,{'color':'black','label':'$T_{i\perp}$'})
    T_trace.customize_li(3,{'color':'green','label':'$T_{i\parallel}$'})
    
    for j in range(1,num_strides):
        tej      = emoms_munge[j]['epochs']
        Te_perpj = emoms_munge[j]['T_perp']
        Te_parj  = emoms_munge[j]['T_par']
        tij      = imoms_munge[j]['epochs']
        Ti_perpj = imoms_munge[j]['T_perp']
        Ti_parj  = imoms_munge[j]['T_par']
        T_trace.add_line(tej,Te_perpj)
        T_trace.add_line(tej,Te_parj)
        T_trace.add_line(tij,Ti_perpj)
        T_trace.add_line(tij,Ti_parj)    
        T_trace.customize_li(4*j,{'color':'red'})
        T_trace.customize_li(4*j+1,{'color':'blue'})
        T_trace.customize_li(4*j+2,{'color':'black'})
        T_trace.customize_li(4*j+3,{'color':'green'})
        
    T_trace.customize_ax({'ylabel':'%s\\nTemp\\n[eV]'%obs,'ylim':[1,1e4],'yscale':'log'})
    T_trace.show_legend()                          
    return T_trace                                           

###############################################################################
def make_stemperature_panel(ax,obs,smoms_munge,species):
       
    #determine the number of strides
    num_strides = len(smoms_munge)    
       
    #graph and label the first segment
    ts1      = smoms_munge[0]['epochs']
    Ts_perp1 = smoms_munge[0]['T_perp']
    Ts_par1  = smoms_munge[0]['T_par']
    T_trace = Grapher.curves(ax,ts1,Ts_perp1)
    T_trace.add_line(ts1,Ts_par1)
    if species == 'emoms':
        T_trace.customize_li(0,{'color':'red', 'label':'$T_{e\perp}$'})
        T_trace.customize_li(1,{'color':'blue','label':'$T_{e\parallel}$'})
    if species == 'imoms':
        T_trace.customize_li(2,{'color':'black','label':'$T_{i\perp}$'})
        T_trace.customize_li(3,{'color':'green','label':'$T_{i\parallel}$'})
    
    for j in range(1,num_strides):
        tsj      = smoms_munge[j]['epochs']
        Ts_perpj = smoms_munge[j]['T_perp']
        Ts_parj  = smoms_munge[j]['T_par']
        T_trace.add_line(tsj,Ts_perpj)
        T_trace.add_line(tsj,Ts_parj)
        if species == 'emoms':
            T_trace.customize_li(4*j,{'color':'red'})
            T_trace.customize_li(4*j+1,{'color':'blue'})
        if species == 'imoms':
            T_trace.customize_li(4*j+2,{'color':'black'})
            T_trace.customize_li(4*j+3,{'color':'green'})
            
    if species == 'emoms':
        T_trace.customize_ax({'ylabel':'%s\\n$T_e$\\n[eV]'%obs,'ylim':[1,1e4],'yscale':'log'})
    if species == 'imoms':
        T_trace.customize_ax({'ylabel':'%s\\n$T_i$\\n[eV]'%obs,'ylim':[1,1e4],'yscale':'log'})

    T_trace.show_legend()                          
    return T_trace                                           
    
    
###############################################################################
def make_Vvector_panel(ax,obs,emoms_munge,imoms_munge):
    #determine the number of segments
    num_estrides = len(emoms_munge)
    num_istrides = len(imoms_munge)
    if num_estrides != num_istrides:
        print "Burst segments for electrons and ions don't match!"
        print "Terminating with extreme prejudice!!!"
        return False
    else:
        num_strides = num_estrides
       
    #graph and label the first segment
    te1      = emoms_munge[0]['epochs']
    Vxe1     = emoms_munge[0]['bulk_vs'][:,0]
    Vye1     = emoms_munge[0]['bulk_vs'][:,1]
    Vze1     = emoms_munge[0]['bulk_vs'][:,2]    
    ti1      = imoms_munge[0]['epochs']
    Vxi1     = imoms_munge[0]['bulk_vs'][:,0]
    Vyi1     = imoms_munge[0]['bulk_vs'][:,1]
    Vzi1     = imoms_munge[0]['bulk_vs'][:,2]
    V_trace = Grapher.curves(ax,te1,Vxe1)
    V_trace.add_line(te1,Vye1)
    V_trace.add_line(te1,Vze1)
    V_trace.add_line(ti1,Vyi1)
    V_trace.add_line(ti1,Vzi1)
    V_trace.add_line(ti1,Vyi1)
    V_trace.customize_li(0,{'color':'blue',  'label':'Vxe_gse','linestyle':'-'})
    V_trace.customize_li(1,{'color':'green', 'label':'Vye_gse','linestyle':'-'})
    V_trace.customize_li(2,{'color':'red','label':'Vze_gse','linestyle':'-'})
    V_trace.customize_li(3,{'color':'blue',  'label':'Vxi_gse','linestyle':'--'})
    V_trace.customize_li(4,{'color':'green', 'label':'Vyi_gse','linestyle':'--'})
    V_trace.customize_li(5,{'color':'red','label':'Vzi_gse','linestyle':'--'})

    for j in range(1,num_strides):
        tej      = emoms_munge[j]['epochs']
        Vxej     = emoms_munge[j]['bulk_vs'][:,0]
        Vyej     = emoms_munge[j]['bulk_vs'][:,1]
        Vzej     = emoms_munge[j]['bulk_vs'][:,2]    
        tij      = imoms_munge[j]['epochs']
        Vxij     = imoms_munge[j]['bulk_vs'][:,0]
        Vyij     = imoms_munge[j]['bulk_vs'][:,1]
        Vzij     = imoms_munge[j]['bulk_vs'][:,2]
        V_trace.add_line(tej,Vxej)
        V_trace.add_line(tej,Vyej)
        V_trace.add_line(tej,Vzej)
        V_trace.add_line(tij,Vyij)
        V_trace.add_line(tij,Vzij)
        V_trace.add_line(tij,Vyij)
        V_trace.customize_li(6*j, {'color':'blue',  'linestyle':'-'})
        V_trace.customize_li(6*j+1, {'color':'green', 'linestyle':'-'})
        V_trace.customize_li(6*j+2, {'color':'red','linestyle':'-'})
        V_trace.customize_li(6*j+3, {'color':'blue',  'linestyle':'--'})
        V_trace.customize_li(6*j+4, {'color':'green', 'linestyle':'--'})
        V_trace.customize_li(6*j+5,{'color':'red','linestyle':'--'})
    
    
    V_trace.customize_ax({'ylabel':'%s\\nVelocity\\n[km/s]'%obs,\
                                        'ylim':[-1000,1000]})    
    V_trace.show_legend()
    return V_trace
    
###############################################################################
def make_sVvector_panel(ax,obs,smoms_munge,species):

    #determine the number of segments
    num_strides = len(smoms_munge)
       
    #graph and label the first segment
    ts1      = smoms_munge[0]['epochs']
    Vxs1     = smoms_munge[0]['bulk_vs'][:,0]
    Vys1     = smoms_munge[0]['bulk_vs'][:,1]
    Vzs1     = smoms_munge[0]['bulk_vs'][:,2]    
    V_trace  = Grapher.curves(ax,ts1,Vxs1)
    V_trace.add_line(ts1,Vys1)
    V_trace.add_line(ts1,Vzs1)
    V_trace.customize_li(0,{'color':'blue',  'label':'$V_x (GSE)$','linestyle':'-'})
    V_trace.customize_li(1,{'color':'green', 'label':'$V_y (GSE)$','linestyle':'-'})
    V_trace.customize_li(2,{'color':'red','label':'$V_z (GSE)$','linestyle':'-'})

    for j in range(1,num_strides):
        tsj      = smoms_munge[j]['epochs']
        Vxsj     = smoms_munge[j]['bulk_vs'][:,0]
        Vysj     = smoms_munge[j]['bulk_vs'][:,1]
        Vzsj     = smoms_munge[j]['bulk_vs'][:,2]    
        V_trace.add_line(tsj,Vxsj)
        V_trace.add_line(tsj,Vysj)
        V_trace.add_line(tsj,Vzsj)
        V_trace.customize_li(3*j, {'color':'blue',  'linestyle':'-'})
        V_trace.customize_li(3*j+1, {'color':'green', 'linestyle':'-'})
        V_trace.customize_li(3*j+2, {'color':'red','linestyle':'-'})
  
    if species == 'emoms':
        V_trace.customize_ax({'ylabel':'%s\\n$V_e$\\n[km/s]'%(obs)})
    if species == 'imoms':
        V_trace.customize_ax({'ylabel':'%s\\n$V_i$\\n[km/s]'%(obs)})    
    V_trace.show_legend()
    return V_trace
    
###############################################################################
def make_Jvector_panel(ax,obs,smoms_munge):
    #determine the number of segments
    num_strides = len(smoms_munge)
       
    #graph and label the first segment
    t1      = smoms_munge[0]['epochs']
    Jx1     = smoms_munge[0]['current'][:,0]
    Jy1     = smoms_munge[0]['current'][:,1]
    Jz1     = smoms_munge[0]['current'][:,2]    
    J_trace = Grapher.curves(ax,t1,Jx1)
    J_trace.add_line(t1,Jy1)
    J_trace.add_line(t1,Jz1)
    J_trace.customize_li(0,{'color':'blue',  'label':'$J_x (GSE)$','linestyle':'-'})
    J_trace.customize_li(1,{'color':'green', 'label':'$J_y (GSE)$','linestyle':'-'})
    J_trace.customize_li(2,{'color':'red',   'label':'$J_z (GSE)$','linestyle':'-'})

    for j in range(1,num_strides):
        tj      = smoms_munge[j]['epochs']
        Jxj     = smoms_munge[j]['current'][:,0]
        Jyj     = smoms_munge[j]['current'][:,1]
        Jzj     = smoms_munge[j]['current'][:,2]    
        J_trace.add_line(tj,Jxj)
        J_trace.add_line(tj,Jyj)
        J_trace.add_line(tj,Jzj)
        J_trace.customize_li(3*j,   {'color':'blue',  'linestyle':'-'})
        J_trace.customize_li(3*j+1, {'color':'green', 'linestyle':'-'})
        J_trace.customize_li(3*j+2, {'color':'red','linestyle':'-'})

    
    
    J_trace.customize_ax({'ylabel':'%s\\nCurrent\\n[$\mu A$/$m^2$]'%obs,\
                                        'ylim':[-10,10]})    
    J_trace.show_legend()
    return J_trace    
    
###############################################################################    
def make_Et_panel(fig,ax,obs,smoms_munge,species,scpot,min_val,max_val):
    #determine the number of segments
    num_strides = len(smoms_munge)

    s1 = np.ma.masked_invalid(np.log10(smoms_munge[0]['omnis'])).T
    t1 = smoms_munge[0]['epochs']
    E1 = smoms_munge[0]['ergs'][0,:]
    Et_spec = Grapher.patch(ax,t1,E1,s1,min_val,max_val)

    for j in range(1,num_strides):
        sj = np.ma.masked_invalid(np.log10(smoms_munge[j]['omnis'])).T
        tj = smoms_munge[j]['epochs']
        Ej = smoms_munge[j]['ergs'][0,:]
        ax.pcolormesh(tj,Ej,sj,vmin=min_val,vmax=max_val,cmap=cmap.jet)

    if scpot != 0:
        tpot1 = scpot[0]['epochs']
        vpot1 = scpot[0]['scpot']
        ax.plot(tpot1,vpot1,'k-',linewidth = 3)
        
    cmap.jet.set_bad('k',alpha=1.0) 
    Et_spec.set_colormap(cmap.jet)
    Et_spec.add_colorbar(fig)
    ax.set_yscale('log')
    
    ax.set_ylabel('%s\n%s\nEnergy\n[eV]'%(obs,species))
    Et_spec.cbar.set_label('$keV/(cm^2 s sr keV)$')

    return Et_spec

###############################################################################    
def make_counterstream_panel(fig,ax,obs,smoms_munge,species,min_val,max_val):
    #determine the number of segments
    num_strides = len(smoms_munge)

    c1 = np.ma.masked_invalid(smoms_munge[0]['par']/smoms_munge[0]['anti']).T
    t1 = smoms_munge[0]['epochs']
    E1 = smoms_munge[0]['ergs'][0,:]
    Et_c_spec = Grapher.patch(ax,t1,E1,c1,min_val,max_val)
    
    for j in range(1,num_strides):
        cj = np.ma.masked_invalid(smoms_munge[j]['par']/smoms_munge[j]['anti']).T
        tj = smoms_munge[j]['epochs']
        Ej = smoms_munge[j]['ergs'][0,:]
        ax.pcolormesh(tj,Ej,cj,vmin=min_val,vmax=max_val,cmap=cmap.bwr)
    
    cmap.bwr.set_bad('k',alpha=1.0) 
    cax  = fig.add_axes(Grapher.cbar_position(ax,0.01,0.01))
    cbar = fig.colorbar(Et_c_spec.patch,cax=cax)  
    Et_c_spec.set_colormap(cmap.bwr)
    ax.set_yscale('log')
    ax.set_ylabel('%s\n%s\nEnergy\n[eV]'%(obs,species))
    cbar.set_label('Counter-streaming\nFraction')

    return Et_c_spec
    
###############################################################################    
def make_trap_fraction_panel(fig,ax,obs,smoms_munge,species,min_val,max_val):
    #determine the number of segments
    num_strides = len(smoms_munge)

    f1 = np.ma.masked_invalid(2.0*smoms_munge[0]['perp']/(smoms_munge[0]['anti'] + smoms_munge[0]['par'])).T
    t1 = smoms_munge[0]['epochs']
    E1 = smoms_munge[0]['ergs'][0,:]
    Et_f_spec = Grapher.patch(ax,t1,E1,f1,min_val,max_val)
    
    for j in range(1,num_strides):
        fj = np.ma.masked_invalid(2.0*smoms_munge[j]['perp']/(smoms_munge[j]['anti'] + smoms_munge[j]['par'])).T
        tj = smoms_munge[j]['epochs']
        Ej = smoms_munge[j]['ergs'][0,:]
        ax.pcolormesh(tj,Ej,fj,vmin=min_val,vmax=max_val,cmap=cmap.bwr)
    
    cmap.bwr.set_bad('k',alpha=1.0) 
    cax  = fig.add_axes(Grapher.cbar_position(ax,0.01,0.01))
    cbar = fig.colorbar(Et_f_spec.patch,cax=cax)  
    Et_f_spec.set_colormap(cmap.bwr)
    ax.set_yscale('log')
    ax.set_ylabel('%s\n%s\nEnergy\n[eV]'%(obs,species))
    cbar.set_label('Trapped\nFraction')

    return Et_f_spec   


    
###############################################################################
def make_Bvector_panel(ax,obs,Bmunge):
    #determine the number of segments
    num_strides = len(Bmunge)
    
    #graph and label the first segment
    t0  = Bmunge[0]['epochs']
    Bx0 = Bmunge[0]['Bgsm'][:,0]
    By0 = Bmunge[0]['Bgsm'][:,1]
    Bz0 = Bmunge[0]['Bgsm'][:,2]
    
    Bvector = Grapher.curves(ax,t0,Bx0)
    Bvector.add_line(t0,By0)
    Bvector.add_line(t0,Bz0)
    Bvector.customize_li(0,{'color':'blue',  'label':'$B_x (GSE)$'})
    Bvector.customize_li(1,{'color':'green', 'label':'$B_y (GSE)$'})
    Bvector.customize_li(2,{'color':'red',   'label':'$B_z (GSE)$'})   
    
    for i in range(1,num_strides):
        ti  = Bmunge[i]['epochs']
        Bxi = Bmunge[i]['Bgse'][:,0]
        Byi = Bmunge[i]['Bgse'][:,1]
        Bzi = Bmunge[i]['Bgse'][:,2]        
        Bvector.add_line(ti,Bxi)
        Bvector.add_line(ti,Byi)
        Bvector.add_line(ti,Bzi)
        Bvector.customize_li(3*i,{'color':'blue'})
        Bvector.customize_li(3*i+1,{'color':'green'})
        Bvector.customize_li(3*i+2,{'color':'red'})

    ylab = '%s\\nFGM\\n[nT]' % (obs,)
    Bvector.customize_ax({'ylabel':ylab})
    Bvector.show_legend()    
    return Bvector
    
###############################################################################    
def make_psd_panel(fig,ax,obs,psd_munge,field,smoms_munge,type):
    #determine the number of segments
    num_strides = len(psd_munge)

    value_min = -11
    value_max = -4
    s1 = np.ma.masked_invalid(np.log10(psd_munge[0][field])).T
    t1 = psd_munge[0]['epochs']
    f1 = psd_munge[0]['freqs']
    psdt_spec = Grapher.patch(ax,t1,f1,s1,value_min,value_max)
    if type == 'plasma':
        fs1 = smoms_munge[0]['f_ps']
        ts1 = smoms_munge[0]['epochs']
        ax.plot(ts1,fs1,'w-')
    if type == 'cyclotron':
        fs1 = smoms_munge[0]['f_cs']
        ts1 = smoms_munge[0]['epochs']
        ax.plot(ts1,fs1,'w-')
    if type == 'half_cyclotron':
        fs1 = smoms_munge[0]['f_cs']
        ts1 = smoms_munge[0]['epochs']
        ax.plot(ts1,fs1,'w-')
        ax.plot(ts1,0.5*fs1,'w-')
    if type == 'chorus':
        fs1 = smoms_munge[0]['f_cs']
        ts1 = smoms_munge[0]['epochs']
        ax.plot(ts1,fs1,'w-')
        ax.plot(ts1,0.5*fs1,'w-')
        ax.plot(ts1,0.1*fs1,'w-')

        
    for j in range(1,num_strides):
        sj = np.ma.masked_invalid(np.log10(psd_munge[j][field])).T
        tj = psd_munge[j]['epochs']
        fj = psd_munge[j]['freqs']
        ax.pcolormesh(tj,fj,sj,vmin=value_min,vmax=value_max,cmap=cmap.jet)
        if type == 'plasma':
            fsj = smoms_munge[j]['f_ps']
            ts1 = smoms_munge[j]['epochs']
            ax.plot(tsj,fsj,'w-')
        if type == 'cyclotron':
            fsj = smoms_munge[j]['f_cs']
            tsj = smoms_munge[j]['epochs']
            ax.plot(tsj,fsj,'w-')
        if type == 'half_cyclotron':
            fs1 = smoms_munge[0]['f_cs']
            ts1 = smoms_munge[0]['epochs']            
            ax.plot(tsj,fsj,'w-')
            ax.plot(tsj,0.5*fsj,'w-')
        if type == 'chorus':
            fsj = smoms_munge[0]['f_cs']
            tsj = smoms_munge[0]['epochs']
            ax.plot(tsj,fsj,'w-')
            ax.plot(tsj,0.5*fsj,'w-')
            ax.plot(tsj,0.1*fsj,'w-')            
    
    cmap.jet.set_bad('k',alpha=1.0) 
    psdt_spec.set_colormap(cmap.jet)
    psdt_spec.add_colorbar(fig)
    ax.set_yscale('log')
    ax.set_ylabel('%s\nFreq\n[Hz]'%(obs))
    if field == 'bpsd':
        psdt_spec.cbar.set_label('$T^2/Hz$')
    if field == 'epsd':
        psdt_spec.cbar.set_label('$(V/m)^2/Hz$')
    return psdt_spec    


###############################################################################
def make_JdotE_panel(ax,obs,smoms_munge):
    #determine the number of strides 
    num_strides = len(smoms_munge)
       
    #graph and label the first segment
    ts1 = smoms_munge[0]['epochs']
    ns1 = smoms_munge[0]['JdotE']
    n_trace = Grapher.curves(ax,ts1,ns1)
   
    for j in range(1,num_strides):
        tsj = smoms_munge[j]['epochs']
        nsj = smoms_munge[j]['JdotE']
        n_trace.add_line(tsj,nsj)
        n_trace.customize_li(j,{'color':'black'})
    
    n_trace.customize_ax({'ylabel':'%s\\n$J\vec \cdot E\vec$\\n[$TBS$]'%obs})
    return n_trace     
    
    
    #'$10^{%d}$'
    