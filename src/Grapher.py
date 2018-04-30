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

fast_parms = {'time_format'  :mdates.DateFormatter('%H:%M:%S'),
              'time_location':mdates.MinuteLocator([0,10,20,30,40,50])}

fast_delta = dt.timedelta(minutes=10)

tick_length = 0.02

###############################################################################
def add_info_box(fig,geometry,annotation_string,border='off'):
    """A simple helper function to add an axis that bears an annotation
    string - typically something that has the date and the position of the 
    spacecraft in some coordinate system.

       Arguments:
          fig:                 figure object that owns the info box axis
          geometry:            tuple containing (x,y,width,height)
          annotation_string:   string to be placed in the info box
          border (optional):   defaults 'off', setting it 'on' 
                               turns on the facecolor to allow for sizing 
                               and placement

       Returns:
           nothing
       
       Example use:  
           add_info_box(fig,(0.05,0.085,0.01,0.01),'Seconds\nMinutes\nPosition')
       
    """    
    dax = fig.add_axes(geometry)
    quiet_axis(dax)
    if border == 'on':
        #dax.tick_params(top='on',bottom='on',left='on',right='on')
        dax.set_facecolor('#eeddcc')
    dax.annotate(annotation_string,xy=(0.0,0.0))  
    return dax

###############################################################################
def cbar_position(current_ax,offset,cbar_width):
    """A simple helper function that returns the position of a colorbar axis
       (or any new axis) relative to a desired axis.
     
       Arguments:
          current_ax:          axis object relative to which the new axis will
                               be positioned  
          offset:              the spacing in x between the new axis
                               and the current axis (It is a horizontal offset)
          cbar_width:          number, in relative coordinates, between 
                               string to be placed in the info box

       Returns:
           list with the geometry for the new axis
       
       Example use:  
           self.cax       = fig.add_axes(cbar_position(self.ax,self.cbar_off,self.cbar_wth))    
.  
    
    Notes:  The cbar_width determines the width of the new axis, also in figure
    relative coordinates.
    
    The final ingredient is the height of the new axis, which is set to be the 
    same as current_axis.  The height of current_axis is queried by making
    the call x_pos, y_pos, width, height = current_ax.get_position().bounds, 
    which returns a tuple with the lower x,y corner position of current_ax and
    its width and height."""
    
    x_pos, y_pos, width, height = current_ax.get_position().bounds
    return [x_pos + width + offset, y_pos, cbar_width, height]

    
###############################################################################
def customize_axis(ax,ax_parms):
    """A small function to 'set' all the common 'setters' for an axis all at 
       once. 

       Arguments:
          ax:          axis object to customize
          ax_parms:    a dictionary of parameters defined as:
          
           The common set of parameters in ax_parms are:
       
               xlabel
               xlim
               xscale
               xticks
               xticklabels
               ylabel
               ylim
               yscale
               yticks
               yticklabels
       
            These parameters should be keys in a dictionary and the 
            corresponding value what would ordinarily be put into the 
            command.  For example, for an axis ax, the conventional way 
            of setting the xlabel would be:
            
                ax.set_xlabel('foo')
            
            and the new way would be
            
                ax_parms = {'xlabel':'foo'}
                
            or 
                 
                ax_parms = {'xlim',[-1,1]}
            
            then run
            
                customize_ax(ax,ax_parms)
            
            In this way, multiple plots can all share the same style, 
            without redundant coding.

       Returns:
           none
       
       Example use:  
           customize_ax(ax,ax_parms) 
       """
    for k in ax_parms.keys():
        if type(ax_parms[k]) == str:
            cmd = "ax."+"set_"+k+"('%s')" % (ax_parms[k])
        else:
            cmd = "ax."+"set_"+k+"(%s)" % (ax_parms[k])
        #print cmd
        eval(cmd)
    
    
###############################################################################
def customize_line(line,line_parms):
    """A small function to 'set' all the common 'setters' for a line all at 
       once.  The common set of parameters are:
 
        Arguments:
          line:        line object to customize
          line_parms:  a dictionary of parameters defined as:
          
          The common set of parameters in line_parms are:

              color
              label
              linestyle
              linewidth
              marker
              markersize 

          These parameters should be keys in a dictionary and the corresponding value
          what would ordinarily be put into the command.  For example, for an line
          li, the conventional way of setting the color would be:
           
              li.set_color('#aaff99')
           
          and the new way would be
           
              li_parms = {'color':'#aaff99'}
               
          then run
           
              customize_line(li,li_parms)
           
          In this way, multiple plots can all share the same style, without redundant
          coding.
           
           
          Returns:
              none
           
          Example use:  
              customize_line(li,li_parms)
        """
    for k in line_parms.keys():
        if type(line_parms[k]) == str:
            cmd = "line."+"set_"+k+"('%s')" % (line_parms[k])
        else:
            cmd = "line."+"set_"+k+"(%s)" % (line_parms[k])
        #print cmd
        eval(cmd)        


###############################################################################   
def decorate_xmajorticks(fig,ax,descriptors):
    """
    Helper function that decorates the x-ticks with extra labels.  To make
    this work, the fig first must be 'drawn' with a draw_idle() and then the
    labels must be pulled.
    
    Assumes that the x-axis is in epoch format
    
        Arguments:
          fig:            figure that owns the axes
          ax:             axis which owns the ticks
          descriptors:    additional text to be given to the ticks
           
        Returns:
          none
           
        Example use:  
          decorate_xmajorticks(my_fig,my_ax,list-of-strings)    
    
    
    """
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
def info_box_position(ax,info_width,info_height):
    """A simple helper function to determine the geometry of an info box
    relative to an axis.

        Arguments:
            ax:          axis that is relative to the info_box
            info_width:  relative width of the info_box
            info_height: relative height of the info_box
           
        Returns:
            geometry:  usual geometry of x,y position and width, height
           
        Example use:  
            info_box_geometry = info_box_position(ax)

        Notes:
            Currently assumes lower left for the info_box location         
"""
    x_pos, y_pos, width, height = ax.get_position().bounds
    x_iax = x_pos - tick_length - info_width
    y_iax = y_pos - tick_length - info_height
    return [x_iax, y_iax, info_width, info_height]
    
  
    
###############################################################################
def quiet_axis(ax):
    """A simple helper function to turn all default options on an axis to 
    off.  The idea here is to exploit the trick of adding a big axis to a
    pre-existing set of axes so that they all share a common xlabel, ylabel, 
    and title"""
    
    #ax.set_axis_bgcolor('none')
    ax.set_facecolor('none')
    ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    
    
###############################################################################
class line(object):
    def __init__(self,ax,x,y):
        self.ax                  = ax
        self.line                = self.ax.plot(x,y)[0]
        self.line.set_color('blue')
        self.line.set_label('')
        self.line.set_linestyle('-')
        self.line.set_linewidth(2)
        self.line.set_marker('None')
        self.line.set_markersize(10)
    def set_color(self,col):
        self.line.set_color(col)
    def set_label(self,lab):
        self.line.set_label(lab)
    def set_linestyle(self,lstyle):
        self.line.set_linestyle(lstyle)
    def set_linewidth(self,lwidth):
        self.line.set_linewidth(lwidth)
    def set_marker(self,mark):
        self.line.set_marker(mark)
    def set_markersize(self,marksize):
        self.line.set_markersize(marksize)

        
###############################################################################
class curves(object):
    def __init__(self,ax,x,y):
        self.ax    = ax
        self.lines = [line(ax,x,y)]
    def add_line(self,x,y):
        self.lines.append(line(self.ax,x,y))
    def customize_ax(self,ax_parms):
        customize_axis(self.ax,ax_parms)
    def customize_li(self,index,line_parms):
        customize_line(self.lines[index],line_parms)
    def format_ax_time(self,t,major_tick_parms,minor_tick_parms,tstyle='custom'):
        if tstyle == 'custom':
            maj_locator   = major_tick_parms['loc'](major_tick_parms['by'],interval=major_tick_parms['int'])
            maj_formatter = mdates.DateFormatter(major_tick_parms['form'])
            min_locator   = minor_tick_parms['loc'](minor_tick_parms['by'],interval=minor_tick_parms['int'])
            min_formatter = mdates.DateFormatter(minor_tick_parms['form'])
        if 'seconds' in tstyle:
            if tstyle == 'seconds5':
                min_locator   = mdates.SecondLocator(bysecond=range(0,60,5),interval=1)
            if tstyle == 'seconds10':
                min_locator   = mdates.SecondLocator(bysecond=[0,10,20,30,40,50],interval=1)
            if tstyle == 'seconds15':
                min_locator   = mdates.SecondLocator(bysecond=[0,15,30,45],interval=1)
            if tstyle == 'seconds20':
                min_locator   = mdates.SecondLocator(bysecond=[0,20,40],interval=1)
            if tstyle == 'seconds30':
                min_locator   = mdates.SecondLocator(bysecond=[0,30],interval=1)
            maj_locator   = mdates.MinuteLocator()
            maj_formatter = mdates.DateFormatter('%M')
            min_formatter = mdates.DateFormatter('\n%S')
        if 'minutes' in tstyle:
            if tstyle == 'minutes5':
                min_locator   = mdates.MinuteLocator(byminute=range(0,60,5),interval=1)
            if tstyle == 'minutes10':
                min_locator   = mdates.MinuteLocator(byminute=[0,10,20,30,40,50],interval=1)
            if tstyle == 'minutes15':
                min_locator   = mdates.MinuteLocator(byminute=[0,15,30,45],interval=1)
            if tstyle == 'minutes20':
                min_locator   = mdates.MinuteLocator(byminute=[0,20,40],interval=1)
            if tstyle == 'minutes30':
                min_locator   = mdates.MinuteLocator(byminute=[0,30],interval=1)
            maj_locator = mdates.HourLocator()
            maj_formatter = mdates.DateFormatter('%H')
            min_formatter = mdates.DateFormatter('\n%M')
        if 'hours' in tstyle:
            if tstyle == 'hours':
                min_locator   = mdates.HourLocator(byhour=range(24),interval=1)
            if tstyle == 'hours2':
                min_locator   = mdates.HourLocator(byhour=range(24,2),interval=1)
            if tstyle == 'hours3':
                min_locator   = mdates.HourLocator(byhour=range(24,3),interval=1)
            if tstyle == 'hours4':
                min_locator   = mdates.HourLocator(byhour=range(24,4),interval=1)
            if tstyle == 'hours6':
                min_locator   = mdates.HourLocator(byhour=range(24,6),interval=1)
            if tstyle == 'hours8':
                min_locator   = mdates.HourLocator(byhour=range(24,8),interval=1)
            if tstyle == 'hours12':
                min_locator   = mdates.HourLocator(byhour=range(24,12),interval=1)
            maj_locator = mdates.DayLocator()
            maj_formatter = mdates.DateFormatter('%d')
            min_formatter = mdates.DateFormatter('\n%H')
        self.ax.xaxis.set_major_locator(maj_locator)
        self.ax.xaxis.set_major_formatter(maj_formatter)
        self.ax.xaxis.set_minor_locator(min_locator)
        self.ax.xaxis.set_minor_formatter(min_formatter)
    def rotate_xticks(self,angle):
        for tick in self.ax.get_xticklabels():
            tick.set_rotation(angle)
    def show_legend(self,location='best'):
        self.ax.legend(loc=location)
    def show_grid(self):
        self.ax.grid('on')    

        
###############################################################################
class patch(object):
    def __init__(self,ax,x,y,z,val_min,val_max):
        self.x           = x
        self.y           = y
        self.z           = z
        self.ax          = ax
        self.val_min     = val_min
        self.val_max     = val_max
        self.patch       = self.ax.pcolormesh(x,y,z,vmin=self.val_min,vmax=self.val_max)
        self.lines       = []
        self.cax         = 'null'
        self.cbar        = 'null'
        self.cbar_off    = 0.01
        self.cbar_wth    = 0.01   
        self.cbar_format = ticker.FormatStrFormatter('$10^{%d}$')
    def add_line(self,x,y):
        self.lines.append(line(self.ax,x,y))        
    def set_colormap(self,cmap):
        self.patch.set_cmap(cmap=cmap)
    def add_colorbar(self,fig):
        self.cbar_span = np.array(range(self.val_min,self.val_max+1))
        self.cax       = fig.add_axes(cbar_position(self.ax,self.cbar_off,self.cbar_wth))
        self.cbar      = fig.colorbar(self.patch,cax=self.cax,ticks=self.cbar_span,format=self.cbar_format)  
    def format_ax_time(self,t,major_tick_parms,minor_tick_parms,tstyle='custom'):
        if tstyle == 'custom':
            maj_locator   = major_tick_parms['loc'](major_tick_parms['by'],interval=major_tick_parms['int'])
            maj_formatter = mdates.DateFormatter(major_tick_parms['form'])
            min_locator   = minor_tick_parms['loc'](minor_tick_parms['by'],interval=minor_tick_parms['int'])
            min_formatter = mdates.DateFormatter(minor_tick_parms['form'])
        if 'seconds' in tstyle:
            if tstyle == 'seconds5':
                min_locator   = mdates.SecondLocator(bysecond=range(0,60,5),interval=1)
            if tstyle == 'seconds10':
                min_locator   = mdates.SecondLocator(bysecond=[0,10,20,30,40,50],interval=1)
            if tstyle == 'seconds15':
                min_locator   = mdates.SecondLocator(bysecond=[0,15,30,45],interval=1)
            if tstyle == 'seconds20':
                min_locator   = mdates.SecondLocator(bysecond=[0,20,40],interval=1)
            if tstyle == 'seconds30':
                min_locator   = mdates.SecondLocator(bysecond=[0,30],interval=1)
            maj_locator   = mdates.MinuteLocator()
            maj_formatter = mdates.DateFormatter('%M')
            min_formatter = mdates.DateFormatter('\n%S')
        if 'minutes' in tstyle:
            if tstyle == 'minutes5':
                min_locator   = mdates.MinuteLocator(byminute=range(0,60,5),interval=1)
            if tstyle == 'minutes10':
                min_locator   = mdates.MinuteLocator(byminute=[0,10,20,30,40,50],interval=1)
            if tstyle == 'minutes15':
                min_locator   = mdates.MinuteLocator(byminute=[0,15,30,45],interval=1)
            if tstyle == 'minutes20':
                min_locator   = mdates.MinuteLocator(byminute=[0,20,40],interval=1)
            if tstyle == 'minutes30':
                min_locator   = mdates.MinuteLocator(byminute=[0,30],interval=1)
            maj_locator = mdates.HourLocator()
            maj_formatter = mdates.DateFormatter('%H')
            min_formatter = mdates.DateFormatter('\n%M')
        if 'hours' in tstyle:
            if tstyle == 'hours':
                min_locator   = mdates.HourLocator(byhour=range(24),interval=1)
            if tstyle == 'hours2':
                min_locator   = mdates.HourLocator(byhour=range(24,2),interval=1)
            if tstyle == 'hours3':
                min_locator   = mdates.HourLocator(byhour=range(24,3),interval=1)
            if tstyle == 'hours4':
                min_locator   = mdates.HourLocator(byhour=range(24,4),interval=1)
            if tstyle == 'hours6':
                min_locator   = mdates.HourLocator(byhour=range(24,6),interval=1)
            if tstyle == 'hours8':
                min_locator   = mdates.HourLocator(byhour=range(24,8),interval=1)
            if tstyle == 'hours12':
                min_locator   = mdates.HourLocator(byhour=range(24,12),interval=1)
            maj_locator = mdates.DayLocator()
            maj_formatter = mdates.DateFormatter('%d')
            min_formatter = mdates.DateFormatter('\n%H')
        self.ax.xaxis.set_major_locator(maj_locator)
        self.ax.xaxis.set_major_formatter(maj_formatter)
        self.ax.xaxis.set_minor_locator(min_locator)
        self.ax.xaxis.set_minor_formatter(min_formatter)
    def rotate_xticks(self,angle):
        for tick in self.ax.get_xticklabels():
            tick.set_rotation(angle)
    def show_legend(self,location='best'):
        self.ax.legend(loc=location)
    def show_grid(self):
        self.ax.grid('on')    

        
