import sys
import pdb

sys.path.append('c:/Users/cschiff/Documents/GitHub/')
import Fetch_Debugs

base_dir = 'x:/data/ftp/%s/fpi/fast/l2/des-debug/'
obs_list = ['mms1']
obs   = obs_list[0]
ver   = r'3.[0123].[012]'
year  = 2016
months = [1]
for month in months:
    dude = Fetch_Debugs.find_debug_CDFs(obs,ver,base_dir,year,month)

print len(dude)

#pdb.set_trace()
culled_list = Fetch_Debugs.cull_list(dude)

print len(culled_list)