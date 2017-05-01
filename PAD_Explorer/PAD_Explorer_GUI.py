# -*- coding: utf-8 -*-
import Tkinter as Tk
import sys
sys.path.append('c:/users/cschiff/Documents/GitHub/PAD')
import Grapher

###############################################################################
# Widgets go here
###############################################################################
#create the root window
root = Tk.Tk()
Debug_filename           = Tk.StringVar()
Dist_filename            = Tk.StringVar()
Photocorrection_filename = Tk.StringVar()
Observatory              = Tk.StringVar()
Species_type             = Tk.StringVar()
CDF_Ver                  = Tk.StringVar()
Corrections_on_flag      = Tk.StringVar()
Corrections_override_val = Tk.StringVar()

#specify some of the root characteristics
root.title("PAD Explorer")
root.geometry("1050x260")


########################################
#  Frame the variables
########################################
#Now define a frame to hold most everything else
all_frame = Tk.Frame(root)
all_frame.configure(width = 980,bd=1,relief='raised')
all_frame.grid(row = 1, column = 0)

#### 0th row
####

#Debug Filename
debug_fname_lab = Tk.Label(all_frame,text = 'Debug Filename')
debug_fname_lab.grid(row = 1, column = 0, sticky = 'W')
debug_fname_entry = Tk.Entry(all_frame,textvariable=Debug_filename)
debug_fname_entry.configure(state='normal',width = 80,justify='left')
debug_fname_entry.grid(row = 1, column = 1, sticky = 'W')

root.mainloop() 