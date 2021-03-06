import spacepy.pycdf as pycdf
import sys
sys.path.append('c:/users/cschiff/Documents/GitHub/PAD/')

import PAD
import Grapher

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

num_time_steps = len(core_data['dist_data']['Epoch'])

time_interval_str = 'Burst Interval 08_09_2017 09_20_44'
Elow              = 0
Ehigh             = 31
file_path         = 'c:/Users/cschiff/Documents/GitHub/PAD/'

for time_label in range(33,num_time_steps):
    Grapher.create_raw_survey_PAD_plot(mode,time_label,Elow,Ehigh,time_interval_str,file_path,core_data)