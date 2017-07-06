# -*- coding: utf-8 -*-
import datetime           as dt
import matplotlib.cm      as cmap
import matplotlib.dates   as mdates
import matplotlib.pyplot  as plt
import matplotlib.ticker  as ticker
import numpy              as np

import Convert
import Grapher
import PAD

###############################################################################
def add_info_box(fig,date,geometry):
    dax = fig.add_axes(geometry)
    quiet_axis(dax)
    string = 'hh:mm\nX-GSM (Re)\nY-GSM (Re)\nZ-GSM (Re)\n%s' % date
    dax.annotate(string,xy=(0.5,0.5))    

           
###############################################################################
def make_density_panel(ax,obs,e_t,e_n,i_t,i_n):
    n_trace = Grapher.traces(ax,e_t,e_n)
    n_trace.add_line(i_t,i_n)
    n_trace.customize_line(0,{'color':'black','label':'Ne'})
    n_trace.customize_line(1,{'color':'green','label':'Ni'})
    n_trace.customize_ax({'loc':'best','xlabel':'','xscale':'','xlim':[],'ylabel':'%s\nDen\n[cm^-3]'%obs,'ylim':'','yscale':''})

###############################################################################
def make_Tperp_panel(ax,obs,e_t,e_Tperp,i_t,i_Tperp):
     Tperp_trace = Grapher.traces(ax,e_t,e_Tperp)
     Tperp_trace.add_line(i_t,i_Tperp)
     Tperp_trace.customize_line(0,{'color':'black','label':'Te_perp'})
     Tperp_trace.customize_line(1,{'color':'green','label':'Ti_perp'})
     Tperp_trace.customize_ax({'loc':'best','xlabel':'','xscale':'','xlim':[],'ylabel':'%s\nTemp\n[eV]'%obs,'ylim':[1e2,1e4],'yscale':'log'})

###############################################################################
def make_eVvector_panel(ax,obs,e_t,e_V):
    eVvector = Grapher.traces(ax,e_t,e_V[:,0])
    eVvector.add_line(e_t,e_V[:,1])
    eVvector.add_line(e_t,e_V[:,2])
    eVvector.customize_line(0,{'color':'blue',  'label':'Vx_GSE'})
    eVvector.customize_line(1,{'color':'green', 'label':'Vy_GSE'})
    eVvector.customize_line(2,{'color':'red',   'label':'Vz_GSE'})
    eVvector.customize_ax({'loc':'best','xlabel':'','xscale':'','xlim':[],'ylabel':'%s\nDES Velocity\n[km/s]'%obs,'ylim':'','yscale':''})
    
###############################################################################
def make_iVvector_panel(ax,obs,i_t,i_V):
    iVvector = Grapher.traces(ax,i_t,i_V[:,0])
    iVvector.add_line(i_t,i_V[:,1])
    iVvector.add_line(i_t,i_V[:,2])
    iVvector.customize_line(0,{'color':'blue',  'label':'Vx_GSE'})
    iVvector.customize_line(1,{'color':'green', 'label':'Vy_GSE'})
    iVvector.customize_line(2,{'color':'red',   'label':'Vz_GSE'})
    iVvector.customize_ax({'loc':'best','xlabel':'','xscale':'','xlim':[],'ylabel':'%s\nDIS Velocity\n[km/s]'%obs,'ylim':'','yscale':''})

###############################################################################    
def make_Et_panel(fig,ax,obs,t,E,s,sc_pot):
    Et_spec = Grapher.patch(ax,t,E,np.ma.masked_invalid(np.log10(s)).T,4,8)
    cmap.jet.set_bad('k',alpha=1.0) 
    Et_spec.set_colormap(cmap.jet)
    Et_spec.add_colorbar(fig)
    ax.set_yscale('log')
    ax.plot(t,sc_pot,'w-',linewidth = 3)
    ax.set_ylabel('Energy\[eV]')
    Et_spec.cbar.set_label('Stuff')
    
###############################################################################
def make_Bvector_panel(ax,obs,t,B):
    Bvector = Grapher.traces(ax,t,B[:,0])
    Bvector.add_line(t,B[:,1])
    Bvector.add_line(t,B[:,2])
    Bvector.customize_line(0,{'color':'blue',  'label':'Bx_GSE'})
    Bvector.customize_line(1,{'color':'green', 'label':'By_GSE'})
    Bvector.customize_line(2,{'color':'red',   'label':'Bz_GSE'})
    Bvector.customize_ax({'loc':'best','xlabel':'','xscale':'','xlim':[],'ylabel':'%s\nFGM\n[nT]'%obs,'ylim':[],'yscale':''})
    Bvector.format_ax_time(t,'brst')
    ax.set_xlabel('Epoch')
