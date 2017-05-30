# -*- coding: utf-8 -*-
import spacepy.pycdf as pycdf
import sys
sys.path.append('c:/users/cschiff/Documents/GitHub/PAD/')
#import Grapher
import PAD

debug_filename = 'C:\Yuggoth\BBF\mms1_fpi_brst_l2_des-debug_20160809092044_v3.1.1.cdf'.replace('\\','/')
dist_filename  = 'C:\Yuggoth\BBF\mms1_fpi_brst_l2_des-dist_20160809092044_v3.1.1.cdf'.replace('\\','/')
photo_filename = 'C:\Yuggoth\Photoelectron Model\mms_fpi_brst_l2_des-bgdist_v1.1.0_p0-2.cdf'.replace('\\','/')

cdf_dict       = {'debug' : pycdf.CDF(debug_filename),
                  'dist'  : pycdf.CDF(dist_filename),
                  'photo' : pycdf.CDF(photo_filename)}
obs            = 'mms1'
mode           = 'brst'
species        = 'des'
ver            = 'ver3'
corrections_on = 0
core_data      = PAD.load_particle_data(cdf_dict,obs,mode,species,ver,corrections_on)