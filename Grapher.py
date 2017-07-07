# -*- coding: utf-8 -*-
import datetime           as dt
import matplotlib.cm      as cmap
import matplotlib.dates   as mdates
import matplotlib.pyplot  as plt
import matplotlib.ticker  as ticker
import numpy              as np

import Convert
import PAD

brst_parms = {'time_format'  :mdates.DateFormatter('%H:%M:%S'),
              'time_location':mdates.SecondLocator([0,20,40])}
              
brst_delta = dt.timedelta(seconds=20)

###############################################################################
def add_info_box(fig,date,geometry):
    dax = fig.add_axes(geometry)
    quiet_axis(dax)
    string = 'hh:mm:ss\nX-GSM (Re)\nY-GSM (Re)\nZ-GSM (Re)\n%s' % date
    dax.annotate(string,xy=(0.5,0.5))    

###############################################################################
def cbar_position(current_ax,offset,cbar_width):
    #get the tuple describing the lower x,y corner position
    #and the width and height
    x_pos, y_pos, width, height = current_ax.get_position().bounds
    return [x_pos + width + offset, y_pos, cbar_width, height] 

###############################################################################
def quiet_axis(ax):
    ax.set_axis_bgcolor('none')
    ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    
###############################################################################
class line():
    def __init__(self,ax,x,y):
        self.ax                  = ax
        self.line                = self.ax.plot(x,y)[0]
        self.line.set_color('blue')
        self.line.set_label('')
        self.line.set_linestyle('-')
        self.line.set_linewidth(2)
        self.line.set_marker('None')
        self.line.set_markersize(10)
    def color(self,col):
        self.line.set_color(col)
    def label(self,lab):
        self.line.set_label(lab)
    def linestyle(self,lstyle):
        self.line.set_linestyle(lstyle)
    def linewidth(self,lwidth):
        self.line.set_linewidth(lwidth)
    def marker(self,mark):
        self.line.set_marker(mark)
    def markersize(self,marksize):
        self.line.set_markersize(marksize)
        
###############################################################################
class curves():
    def __init__(self,ax,x,y):
        self.ax    = ax
        self.lines = [line(ax,x,y)]
    def add_line(self,x,y):
        self.lines.append(line(self.ax,x,y))
    def customize_line(self,index,line_parms):
        self.lines[index].color(line_parms['color'])
        self.lines[index].label(line_parms['label'])
    def customize_ax(self,ax_parms):
        self.ax.set_xlabel(ax_parms['xlabel'])
        if len(ax_parms['xlim']) == 2:
            self.ax.set_xlim  (ax_parms['xlim'])
        if ax_parms['xscale'] != '':
            self.ax.set_xscale(ax_parms['xscale'])
        self.ax.set_ylabel(ax_parms['ylabel'])
        if len(ax_parms['ylim']) == 2:        
            self.ax.set_ylim  (ax_parms['ylim'])
        if ax_parms['yscale'] != '':            
            self.ax.set_yscale(ax_parms['yscale'])
        if ax_parms['loc'] != '':
            self.ax.legend(loc=ax_parms['loc'])   
            
###############################################################################
class traces():
    def __init__(self,ax,x,y):
        self.ax    = ax
        self.lines = [line(ax,x,y)]
    def add_line(self,x,y):
        self.lines.append(line(self.ax,x,y))
    def customize_line(self,index,line_parms):
        self.lines[index].color(line_parms['color'])
        self.lines[index].label(line_parms['label'])
    def format_ax_time(self,t,tformat):
        if tformat == 'brst':
            self.ax.xaxis.set_major_locator(brst_parms['time_location'])
            self.ax.xaxis.set_major_formatter(brst_parms['time_format'])
            t0 = Convert.round_time(t[0], date_delta=brst_delta,to='down')
            tf = Convert.round_time(t[-1],date_delta=brst_delta,to='up')
            self.ax.set_xlim([t0,tf])
    def customize_ax(self,ax_parms):
        self.ax.set_xlabel(ax_parms['xlabel'])
        if len(ax_parms['xlim']) == 2:
            self.ax.set_xlim  (ax_parms['xlim'])
        self.ax.set_xscale('linear')
        self.ax.set_ylabel(ax_parms['ylabel'])
        if len(ax_parms['ylim']) == 2:        
            self.ax.set_ylim  (ax_parms['ylim'])
        if ax_parms['yscale'] != '':            
            self.ax.set_yscale(ax_parms['yscale'])
        if ax_parms['loc'] != '':
            self.ax.legend(loc=ax_parms['loc'])  

###############################################################################
class patch():
    def __init__(self,ax,x,y,z,val_min,val_max):
        self.x           = x
        self.y           = y
        self.z           = z
        self.ax          = ax
        self.val_min     = val_min
        self.val_max     = val_max
        self.patch       = self.ax.pcolormesh(x,y,z,vmin=self.val_min,vmax=self.val_max)
        self.cax         = 'null'
        self.cbar        = 'null'
        self.cbar_off    = 0.01
        self.cbar_wth    = 0.01   
        self.cbar_format = ticker.FormatStrFormatter('$10^{%d}$')
    def set_colormap(self,cmap):
        self.patch.set_cmap(cmap=cmap)
    def add_colorbar(self,fig):
        self.cbar_span = np.array(range(self.val_min,self.val_max+1))
        self.cax       = fig.add_axes(cbar_position(self.ax,self.cbar_off,self.cbar_wth))
        self.cbar      = fig.colorbar(self.patch,cax=self.cax,ticks=self.cbar_span,format=self.cbar_format)  
        
#################################################################################