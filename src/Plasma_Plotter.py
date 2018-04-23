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
        n_trace.customize_li(j+1,{'color':'black'})
        n_trace.customize_li(j+2,{'color':'green'})
        
    n_trace.customize_ax({'ylabel':'%s\\ndensity\\n[cm^-3]'%obs})
    n_trace.show_legend()
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
    T_trace.customize_li(0,{'color':'red','label':'Te_perp'})
    T_trace.customize_li(1,{'color':'blue','label':'Te_par'})
    T_trace.customize_li(2,{'color':'black','label':'Ti_perp'})
    T_trace.customize_li(3,{'color':'green','label':'Ti_par'})
    
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
        T_trace.customize_li(j+3,{'color':'red'})
        T_trace.customize_li(j+4,{'color':'blue'})
        T_trace.customize_li(j+5,{'color':'black'})
        T_trace.customize_li(j+6,{'color':'green'})
        
    T_trace.customize_ax({'ylabel':'%s\\nTemp\\n[eV]'%obs,'ylim':[1,1e4],'yscale':'log'})
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
    V_trace.customize_li(0,{'color':'red',  'label':'Vxe_gse','linestyle':'-'})
    V_trace.customize_li(1,{'color':'blue', 'label':'Vye_gse','linestyle':'-'})
    V_trace.customize_li(2,{'color':'green','label':'Vze_gse','linestyle':'-'})
    V_trace.customize_li(3,{'color':'red',  'label':'Vxi_gse','linestyle':'--'})
    V_trace.customize_li(4,{'color':'blue', 'label':'Vyi_gse','linestyle':'--'})
    V_trace.customize_li(5,{'color':'green','label':'Vzi_gse','linestyle':'--'})

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
        V_trace.customize_li(j+5, {'color':'red',  'linestyle':'-'})
        V_trace.customize_li(j+6, {'color':'blue', 'linestyle':'-'})
        V_trace.customize_li(j+7, {'color':'green','linestyle':'-'})
        V_trace.customize_li(j+8, {'color':'red',  'linestyle':'--'})
        V_trace.customize_li(j+9, {'color':'blue', 'linestyle':'--'})
        V_trace.customize_li(j+10,{'color':'green','linestyle':'--'})
    
    
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
    V_trace.customize_li(0,{'color':'red',  'label':'Vxe_gse','linestyle':'-'})
    V_trace.customize_li(1,{'color':'blue', 'label':'Vye_gse','linestyle':'-'})
    V_trace.customize_li(2,{'color':'green','label':'Vze_gse','linestyle':'-'})

    for j in range(1,num_strides):
        tsj      = smoms_munge[j]['epochs']
        Vxsj     = smoms_munge[j]['bulk_vs'][:,0]
        Vysj     = smoms_munge[j]['bulk_vs'][:,1]
        Vzsj     = smoms_munge[j]['bulk_vs'][:,2]    
        V_trace.add_line(tsj,Vxsj)
        V_trace.add_line(tsj,Vysj)
        V_trace.add_line(tsj,Vzsj)
        V_trace.customize_li(j+2, {'color':'red',  'linestyle':'-'})
        V_trace.customize_li(j+3, {'color':'blue', 'linestyle':'-'})
        V_trace.customize_li(j+4, {'color':'green','linestyle':'-'})
  
    
    V_trace.customize_ax({'ylabel':'%s\\n%s Velocity\\n[km/s]'%(obs,species)})    
    V_trace.show_legend()
    return V_trace
    
###############################################################################    
def make_Et_panel(fig,ax,obs,smoms_munge,species,sc_pot):
    #determine the number of segments
    num_strides = len(smoms_munge)

    s1 = np.ma.masked_invalid(np.log10(smoms_munge[0]['omnis'])).T
    t1 = smoms_munge[0]['epochs']
    E1 = smoms_munge[0]['ergs'][0,:]
    Et_spec = Grapher.patch(ax,t1,E1,s1,4,8)
    
    for j in range(1,num_strides):
        sj = np.ma.masked_invalid(np.log10(smoms_munge[j]['omnis'])).T
        tj = smoms_munge[j]['epochs']
        Ej = smoms_munge[j]['ergs'][0,:]
        ax.pcolormesh(tj,Ej,sj,vmin=4,vmax=8,cmap=cmap.jet)
    
    cmap.jet.set_bad('k',alpha=1.0) 
    Et_spec.set_colormap(cmap.jet)
    Et_spec.add_colorbar(fig)
    ax.set_yscale('log')
    #ax.plot(t,sc_pot,'w-',linewidth = 3)
    ax.set_ylabel('%s\n%s\nEnergy\n[eV]'%(obs,species))
    Et_spec.cbar.set_label('keV/(cm^2 s sr keV)')

    return Et_spec
    
###############################################################################
def make_Bvector_panel(ax,obs,Bmunge):
    #determine the number of segments
    num_strides = len(Bmunge)
    
    #graph and label the first segment
    t0  = Bmunge[0]['epochs']
    Bx0 = Bmunge[0]['Bgse'][:,0]
    By0 = Bmunge[0]['Bgse'][:,1]
    Bz0 = Bmunge[0]['Bgse'][:,2]
    
    Bvector = Grapher.curves(ax,t0,Bx0)
    Bvector.add_line(t0,By0)
    Bvector.add_line(t0,Bz0)
    Bvector.customize_li(0,{'color':'blue',  'label':'Bx_GSE'})
    Bvector.customize_li(1,{'color':'green', 'label':'By_GSE'})
    Bvector.customize_li(2,{'color':'red',   'label':'Bz_GSE'})   
    
    for i in range(1,num_strides):
        ti  = Bmunge[i]['epochs']
        Bxi = Bmunge[i]['Bgse'][:,0]
        Byi = Bmunge[i]['Bgse'][:,1]
        Bzi = Bmunge[i]['Bgse'][:,2]        
        Bvector.add_line(ti,Bxi)
        Bvector.add_line(ti,Byi)
        Bvector.add_line(ti,Bzi)
        Bvector.customize_li(i+2,{'color':'blue'})
        Bvector.customize_li(i+3,{'color':'green'})
        Bvector.customize_li(i+4,{'color':'red'})

    ylab = '%s\\nFGM\\n[nT]' % (obs,)
    Bvector.customize_ax({'ylabel':ylab})
    Bvector.show_legend()    
    return Bvector
    